"""Retrospective bound-component controls, with reference withheld from assembly and ranking."""

import csv
import gzip
import json
import re
import shutil
import time
import urllib.request
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial.distance import cdist

from .chemistry import read_ccd
from .engine import assemble, digest, load_inputs
from .geometry import Distances, fit, interface, rmsd
from .structures import Cloud, read_cloud, write_complex


def fetch(url, target):
    """One explicit request at a time; atomic download, no silent empty response."""
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if target.stat().st_size == 0:
            raise ValueError(f"Empty cached file: {target.name}")
        return target
    time.sleep(0.25)
    request = urllib.request.Request(url, headers={"User-Agent": "TOPICS/0.1 structural-research"})
    with urllib.request.urlopen(request, timeout=60) as response:
        content = response.read(30_000_001)
    if not content or len(content) > 30_000_000:
        raise ValueError("Unexpected download size; retrieve the structure manually.")
    temporary = target.with_suffix(target.suffix + ".part")
    temporary.write_bytes(content)
    temporary.replace(target)
    return target


def prepare(spec_path, output, data_dir=None):
    spec = json.loads(Path(spec_path).read_text())
    code, ligand_code = spec["pdb_id"].upper(), spec["ligand"].upper()
    if not re.fullmatch(r"[0-9][A-Z0-9]{3}", code) or not re.fullmatch(r"[A-Z0-9]{1,5}", ligand_code):
        raise ValueError("Invalid PDB or chemical component identifier.")
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Benchmark preparation requires a new directory.")
    output.mkdir(parents=True, exist_ok=True)
    if data_dir:
        data_dir = Path(data_dir)
        cif = data_dir / f"{code}.cif"
        if not cif.exists():
            archive = data_dir / f"pdb_0000{code.lower()}.cif.gz"
            if not archive.exists():
                raise ValueError(f"Missing local structure {code}.cif or its compressed archive.")
            cif = output / "source.cif"
            cif.write_bytes(gzip.decompress(archive.read_bytes()))
        ccd = data_dir / f"{ligand_code}.cif"
    else:
        cif = fetch(f"https://files.rcsb.org/download/{code}.cif", output / "source.cif")
        ccd = fetch(f"https://files.rcsb.org/ligands/download/{ligand_code}.cif", output / "ligand.cif")
    molecule, names = read_ccd(ccd)
    poi = read_cloud(cif, spec["poi_chains"])
    e3 = read_cloud(cif, spec["e3_chains"])
    ligand = read_cloud(cif, [spec["ligand_chain"]], ligand_code, spec["ligand_residue"])
    coordinates = {key[-1]: xyz for key, xyz in zip(ligand.keys, ligand.xyz)}
    if set(coordinates) != set(names):
        raise ValueError("Crystal and CCD heavy atoms do not match exactly; curate missing or renamed atoms first.")
    native_ligand = np.array([coordinates[name] for name in names])
    poi_indices = [names.index(n) for n in spec["poi_head_atoms"]]
    e3_indices = [names.index(n) for n in spec["e3_head_atoms"]]
    if set(poi_indices) & set(e3_indices):
        raise ValueError("Head atoms must not overlap.")
    # Each binary component has its own origin. The native relative translation is discarded.
    center_poi = native_ligand[poi_indices].mean(0)
    center_e3 = native_ligand[e3_indices].mean(0)
    write_complex(output / "poi.cif", {"POI": poi.moved(np.eye(3), -center_poi)})
    write_complex(output / "e3.cif", {"E3": e3.moved(np.eye(3), -center_e3)})
    if ccd.resolve() != (output / "ligand.cif").resolve():
        shutil.copyfile(ccd, output / "ligand.cif")
    config = {"schema_version": 1, "title": spec["title"], "ligand_ccd": "ligand.cif",
              "source": {"pdb_id": code, "url": f"https://www.rcsb.org/structure/{code}",
                         "sha256": digest(cif), "resolution_A": gemmi.read_structure(str(cif)).resolution,
                         "benchmark_type": "Retrospective bound-component assembly; not an independent blind benchmark",
                         "head_definition": "Explicit atom selections curated from the known ligand binding heads",
                         "selection": spec, "excluded_components": spec.get("excluded_components", [])},
              "poi": {"structure": "poi.cif", "head_atom_names": spec["poi_head_atoms"],
                      "head_coordinates": (native_ligand[poi_indices] - center_poi).tolist()},
              "e3": {"structure": "e3.cif", "head_atom_names": spec["e3_head_atoms"],
                     "head_coordinates": (native_ligand[e3_indices] - center_e3).tolist()}}
    config_path = output / "input.json"
    config_path.write_text(json.dumps(config, indent=2))
    inputs = load_inputs(config_path)
    np.savez_compressed(output / "reference.npz", poi=poi.xyz, e3=e3.xyz, ligand=native_ligand,
                        poi_keys=np.asarray(inputs.poi.keys, dtype=str), e3_keys=np.asarray(inputs.e3.keys, dtype=str),
                        atom_names=np.asarray(names))
    # Positive control exercises the same rigid assembly using the native linker geometry.
    recovered_e3, recovered_ligand, head_error = assemble(inputs, native_ligand)
    r, t, _ = fit(inputs.poi.xyz, poi.xyz)
    positive = {"e3_rmsd_A": rmsd(recovered_e3.xyz @ r + t, e3.xyz),
                "ligand_rmsd_A": rmsd(recovered_ligand @ r + t, native_ligand), "head_rmsd_A": head_error,
                "meaning": "Geometry/mapping positive control using native ligand geometry; not predictive accuracy"}
    (output / "positive_control.json").write_text(json.dumps(positive, indent=2))
    return config_path


def evaluate(run_dir, reference_path):
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "manifest.json").read_text())
    if manifest.get("status") != "complete":
        raise ValueError("Only completed runs can be evaluated.")
    reference = np.load(reference_path, allow_pickle=False)
    arrays = np.load(run_dir / "ensemble.npz", allow_pickle=False)
    for role in ("poi", "e3"):
        if reference[f"{role}_keys"].tolist() != [[str(v) for v in k] for k in manifest["atom_keys"][role]]:
            raise ValueError(f"Reference {role} atom identities/order do not match this run.")
    if reference["atom_names"].tolist() != manifest["atom_keys"]["ligand"]:
        raise ValueError("Reference ligand atom mapping differs.")
    for role in ("poi", "e3", "ligand"):
        n = len(manifest["atom_keys"][role])
        if reference[role].shape != (n, 3) or not np.isfinite(reference[role]).all():
            raise ValueError(f"Invalid reference {role} coordinates.")
    r, t, _ = fit(arrays["poi"], reference["poi"])
    keys = {role: [tuple(k) for k in manifest["atom_keys"][role]] for role in ("poi", "e3")}
    ref_poi = Cloud(reference["poi"], keys["poi"], ["C"] * len(keys["poi"]))
    ref_e3 = Cloud(reference["e3"], keys["e3"], ["C"] * len(keys["e3"]))
    backend = Distances("cpu")
    native_contacts, _ = interface(ref_poi, ref_e3, backend)
    ca = np.array([k[-1] == "CA" for k in keys["e3"]])
    if not ca.any() or not native_contacts:
        raise ValueError("Evaluation needs E3 alpha carbons and a non-empty native PPI interface.")
    rows = list(csv.DictReader((run_dir / "candidates.csv").open()))
    evaluated = []
    for index, row in enumerate(rows):
        transformed = arrays["e3"][index] @ r + t
        predicted = Cloud(transformed, keys["e3"], ref_e3.elements)
        contacts, _ = interface(ref_poi, predicted, backend)
        shared = len(contacts & native_contacts)
        evaluated.append({"candidate": row["candidate"], "rank": int(row["rank"]), "seed": int(row["seed"]),
                          "e3_ca_rmsd_A": rmsd(transformed[ca], reference["e3"][ca]),
                          "ligand_rmsd_A": rmsd(arrays["ligand"][index] @ r + t, reference["ligand"]),
                          "native_contact_recall": shared / len(native_contacts),
                          "contact_precision": shared / len(contacts) if contacts else 0.0})
    errors = [row["e3_ca_rmsd_A"] for row in evaluated]
    summary = {"type": "Retrospective bound-component assembly",
               "reference_sha256": digest(reference_path), "reference_used_for_ranking": False,
               "top_ranked": evaluated[0], "best_sampled": min(evaluated, key=lambda x: x["e3_ca_rmsd_A"]),
               "median_e3_ca_rmsd_A": float(np.median(errors)), "candidate_count": len(evaluated),
               "native_ppi_contacts": len(native_contacts), "candidates": evaluated,
               "per_seed": [{"seed": seed, "top_ranked": next(x for x in evaluated if x["seed"] == seed),
                             "best_sampled_e3_ca_rmsd_A": min(x["e3_ca_rmsd_A"] for x in evaluated if x["seed"] == seed)}
                            for seed in sorted({x["seed"] for x in evaluated})],
               "limitations": "Bound binary shapes and curated head atom identities come from this crystal. Two controls do not establish general accuracy. RMSD uses exact atom names without symmetry permutations; contacts are heavy-atom residue pairs within 5 A. E3 CA RMSD is measured after fitting POI heavy atoms; it is not DockQ."}
    (run_dir / "evaluation.json").write_text(json.dumps(summary, indent=2, allow_nan=False))
    with (run_dir / "evaluation.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(evaluated[0])); writer.writeheader(); writer.writerows(evaluated)
    # Bring the reference into the input POI frame for the report viewer.
    np.savez_compressed(run_dir / "reference_view.npz", poi=(reference["poi"] - t) @ r.T,
                        e3=(reference["e3"] - t) @ r.T, ligand=(reference["ligand"] - t) @ r.T)
    from .report import render
    render(run_dir)
    return summary
