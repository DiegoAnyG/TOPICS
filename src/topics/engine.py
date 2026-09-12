"""A transparent baseline: place rigid binding partners on independently sampled PROTACs."""

import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem
from threadpoolctl import threadpool_limits

from . import __version__
from .chemistry import descriptors, generate, read_ccd
from .geometry import Distances, fit, interface, rmsd
from .structures import Cloud, read_cloud, write_complex


@dataclass
class Inputs:
    poi: Cloud
    e3: Cloud
    mol: object
    atom_names: list[str]
    poi_indices: list[int]
    e3_indices: list[int]
    poi_head: np.ndarray
    e3_head: np.ndarray
    title: str
    provenance: dict


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs(config_path):
    config_path = Path(config_path)
    data = json.loads(config_path.read_text())
    if data.get("schema_version") != 1:
        raise ValueError("Unsupported input schema; expected schema_version=1.")
    base = config_path.parent
    mol, names = read_ccd(base / data["ligand_ccd"])
    clouds, indices, heads = {}, {}, {}
    hashes = {"config": digest(config_path), "ligand_ccd": digest(base / data["ligand_ccd"])}
    for role in ("poi", "e3"):
        spec = data[role]
        clouds[role] = read_cloud(base / spec["structure"])
        hashes[role] = digest(base / spec["structure"])
        atom_names = spec["head_atom_names"]
        if len(atom_names) != len(set(atom_names)) or len(atom_names) < 3:
            raise ValueError(f"{role}: at least three unique attachment atoms are required.")
        try:
            indices[role] = [names.index(name) for name in atom_names]
        except ValueError as exc:
            raise ValueError(f"{role}: an attachment atom is absent from the chemical dictionary.") from exc
        heads[role] = np.asarray(spec["head_coordinates"], dtype=float)
        if heads[role].shape != (len(atom_names), 3):
            raise ValueError(f"{role}: attachment coordinates do not match atom names.")
        fit(heads[role], heads[role])
    if set(indices["poi"]) & set(indices["e3"]):
        raise ValueError("The two binding heads must have disjoint atom mappings.")
    return Inputs(clouds["poi"], clouds["e3"], mol, names, indices["poi"], indices["e3"],
                  heads["poi"], heads["e3"], str(data.get("title", "TOPICS assembly")),
                  {"hashes": hashes, "source": data.get("source", {}), "input_schema": 1})


def assemble(inputs, conformer):
    """Only local binary poses enter this function; no reference ternary is accepted."""
    r1, t1, error1 = fit(conformer[inputs.poi_indices], inputs.poi_head)
    ligand = conformer @ r1 + t1
    r2, t2, error2 = fit(inputs.e3_head, ligand[inputs.e3_indices])
    return inputs.e3.moved(r2, t2), ligand, max(error1, error2)


def ligand_cloud(inputs, xyz):
    code = inputs.mol.GetProp("ccd_id") if inputs.mol.HasProp("ccd_id") else "LIG"
    return Cloud(xyz, [("L", 1, "", code, n) for n in inputs.atom_names],
                 [atom.GetSymbol() for atom in inputs.mol.GetAtoms()])


def rank_candidates(rows):
    """Lexicographic feasibility ordering; never inspect reference-based metrics."""
    rows.sort(key=lambda r: (r["clash_pairs"], r["head_rmsd_A"], -r["ppi_contacts"], r["energy_kcal_mol"], r["candidate"]))
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank


def run(config_path, output, count=16, seeds=(42,), threads=1, device="cpu", pause=0.05, progress=print):
    if not 1 <= count <= 1000 or not 1 <= threads <= min(8, os.cpu_count() or 1):
        raise ValueError("Use 1–1000 conformers and 1–8 threads, within the available CPU count.")
    if not seeds or len(seeds) > 20 or len(set(seeds)) != len(seeds) or any(not 0 <= s < 2**31 for s in seeds):
        raise ValueError("Supply 1–20 unique integer seeds between 0 and 2147483647.")
    if not np.isfinite(pause) or not 0 <= pause <= 60:
        raise ValueError("Pause must be between 0 and 60 seconds.")
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory is not empty; choose a new run directory.")
    inputs = load_inputs(config_path)
    estimated_bytes = count * len(seeds) * (len(inputs.e3.xyz) + len(inputs.atom_names)) * 3 * 8 * 3
    if estimated_bytes > 256 * 1024**2:
        raise ValueError("Requested ensemble exceeds the 256 MiB coordinate budget. Reduce conformers or seeds.")
    backend = Distances(device)
    output.mkdir(parents=True, exist_ok=True)
    (output / "manifest.json").write_text(json.dumps({"schema_version": 1, "status": "running",
        "provenance": inputs.provenance, "parameters": {"seeds": list(seeds), "conformers_per_seed": count}}))
    start = time.perf_counter()
    rows, poses, stats = [], {}, []
    progress(f"{inputs.title}: {backend.reason}; {threads} CPU thread(s).")
    with threadpool_limits(limits=threads):
        for seed in seeds:
            progress(f"Generating {count} conformers, seed {seed}...")
            ensemble, seed_stats = generate(inputs.mol, count, seed, threads)
            stats.append({"seed": seed, **seed_stats})
            for cid, xyz, energy in ensemble:
                e3, ligand, error = assemble(inputs, xyz)
                ppi, pp_clashes = interface(inputs.poi, e3, backend)
                lig = ligand_cloud(inputs, ligand)
                pl, pl_clashes = interface(inputs.poi, lig, backend)
                el, el_clashes = interface(e3, lig, backend)
                candidate = f"s{seed}_c{cid}"
                rows.append({"candidate": candidate, "seed": seed, "conformer": cid,
                             "clash_pairs": pp_clashes + pl_clashes + el_clashes,
                             "ppi_contacts": len(ppi), "poi_ligand_contacts": len(pl),
                             "e3_ligand_contacts": len(el), "head_rmsd_A": error, "energy_kcal_mol": energy})
                poses[candidate] = (e3, ligand)
                if pause:
                    time.sleep(pause)
    rank_candidates(rows)
    elapsed = time.perf_counter() - start
    versions = {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "rdkit", "gemmi", "plotly", "matplotlib")}
    if backend.device == "cuda":
        versions["cupy-cuda12x"] = importlib.metadata.version("cupy-cuda12x")
    source_hash = hashlib.sha256()
    for source_file in sorted(Path(__file__).parent.glob("*.py")):
        source_hash.update(source_file.name.encode() + b"\0" + source_file.read_bytes())
    manifest = {"schema_version": 1, "topics_version": __version__, "implementation_sha256": source_hash.hexdigest(),
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "title": inputs.title, "provenance": inputs.provenance, "versions": versions,
                "platform": {"system": platform.system(), "machine": platform.machine(), "python": platform.python_version()},
                "parameters": {"conformers_per_seed": count, "seeds": list(seeds), "threads": threads, "device_requested": device,
                               "device_used": backend.device, "pause_seconds": pause, "clash_cutoff_A": 2.0, "contact_cutoff_A": 5.0,
                               "embedding": "ETKDGv3, random coordinates, enforced chirality", "forcefield": "MMFF94s, 1500 iterations"},
                "backend_note": backend.reason, "conformers": stats, "compute_seconds": elapsed,
                "descriptors": descriptors(inputs.mol), "status": "rendering",
                "method": "Rigid bound-component conformer-driven assembly; geometry-only lexicographic ranking",
                "limitations": ["Bound head coordinates and protein conformations are supplied as input.",
                                "No protein flexibility, solvent, PPI energy, ubiquitination or cooperativity prediction.",
                                "Fixed geometric cutoffs are heuristics, not physical energy or calibrated probabilities.",
                                "Low-power defaults limit work; they cannot guarantee a hardware temperature."]}
    with (output / "candidates.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    arrays = {"poi": inputs.poi.xyz, "e3": np.array([poses[r["candidate"]][0].xyz for r in rows]),
              "ligand": np.array([poses[r["candidate"]][1] for r in rows])}
    np.savez_compressed(output / "ensemble.npz", **arrays)
    manifest["atom_keys"] = {"poi": inputs.poi.keys, "e3": inputs.e3.keys, "ligand": inputs.atom_names}
    best_e3, best_ligand = poses[rows[0]["candidate"]]
    write_complex(output / "best.cif", {"POI": inputs.poi, "E3": best_e3, "PROTAC": ligand_cloud(inputs, best_ligand)})
    molecule = Chem.Mol(inputs.mol)
    with Chem.SDWriter(str(output / "ensemble.sdf")) as writer:
        for row in rows:
            molecule.RemoveAllConformers()
            conformer = Chem.Conformer(molecule.GetNumAtoms())
            for i, xyz in enumerate(poses[row["candidate"]][1]):
                conformer.SetAtomPosition(i, xyz)
            molecule.AddConformer(conformer)
            molecule.SetProp("_Name", row["candidate"])
            writer.write(molecule)
    (output / "view.pml").write_text("load best.cif, predicted\nhide everything\nshow cartoon, polymer\nshow sticks, chain PROTAC_L\ncolor marine, chain POI_*\ncolor orange, chain E3_*\nzoom\n")
    manifest["output_hashes"] = {name: digest(output / name) for name in ("candidates.csv", "ensemble.npz", "ensemble.sdf", "best.cif")}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False))
    from .report import render
    render(output)
    manifest["status"] = "complete"
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False))
    progress(f"Saved {len(rows)} candidates in {elapsed:.2f} s. Open {output / 'report.html'}")
    return output
