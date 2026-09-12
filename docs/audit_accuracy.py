"""Audit saved bound-component controls without changing predictions or ranking.

Run with the TOPICS environment. Optional --dockq-python names an isolated
Python environment containing DockQ; DockQ is not a production dependency.
"""

import argparse
from collections import Counter
import csv
import importlib.metadata
import json
from pathlib import Path
import subprocess
import tempfile

import gemmi
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

from topics.engine import digest, ligand_cloud, load_inputs
from topics.geometry import Distances, fit, interface, rmsd
from topics.structures import Cloud, write_complex


def write_pair(path, poi, e3):
    """Export one chain per role; fail instead of silently merging chains."""
    assert len({k[0] for k in poi.keys}) == len({k[0] for k in e3.keys}) == 1
    write_complex(path, {"POI": poi, "E3": e3})
    structure = gemmi.read_structure(str(path))
    for chain, name in zip(structure[0], "AB"):
        chain.name = name
    structure.make_mmcif_document().write_file(str(path))


def audit(case, dockq_python):
    inputs = load_inputs(case / "inputs/input.json")
    run = case / "run"
    manifest = json.loads((run / "manifest.json").read_text())
    evaluation = json.loads((run / "evaluation.json").read_text())
    assert manifest["status"] == "complete"
    assert evaluation["reference_used_for_ranking"] is False
    arrays = np.load(run / "ensemble.npz", allow_pickle=False)
    reference_path = case / "inputs/reference.npz"
    reference = np.load(reference_path, allow_pickle=False)
    assert digest(reference_path) == evaluation["reference_sha256"]
    assert reference["atom_names"].tolist() == inputs.atom_names
    for role in ("poi", "e3"):
        assert reference[f"{role}_keys"].tolist() == [list(map(str, key)) for key in getattr(inputs, role).keys]
    with (run / "candidates.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert [r["candidate"] for r in rows] == [r["candidate"] for r in evaluation["candidates"]]
    assert len(rows) == len(arrays["e3"]) == len(arrays["ligand"])
    native = {role: Cloud(reference[role], getattr(inputs, role).keys, getattr(inputs, role).elements)
              for role in ("poi", "e3")}
    backend = Distances("cpu")
    ligand_distances = {role: cdist(reference["ligand"], reference[role]).min(axis=1)
                        for role in native}
    # Reference-derived diagnostic sets only: these are NOT chemical head definitions
    # and must never become inputs to a blind prediction.
    masks = {role: (ligand_distances[role] <= 4.0) & (ligand_distances[other] > 5.0)
             for role, other in (("poi", "e3"), ("e3", "poi"))}
    assert all(mask.any() for mask in masks.values())
    selected = {"top_ranked": 0,
                "best_e3_rmsd": next(i for i, r in enumerate(rows)
                                     if r["candidate"] == evaluation["best_sampled"]["candidate"])}
    details = {}
    with tempfile.TemporaryDirectory(prefix="topics-audit-") as folder:
        folder = Path(folder)
        native_file = folder / "native.cif"
        write_pair(native_file, native["poi"], native["e3"])
        for label, index in {"native": None, **selected}.items():
            proteins = native if index is None else {
                "poi": inputs.poi,
                "e3": Cloud(arrays["e3"][index], inputs.e3.keys, inputs.e3.elements)}
            ligand = ligand_cloud(inputs, reference["ligand"] if index is None else arrays["ligand"][index])
            pair_counts = {}
            for name, a, b in (("ppi", proteins["poi"], proteins["e3"]),
                               ("poi_ligand", proteins["poi"], ligand),
                               ("e3_ligand", proteins["e3"], ligand)):
                contacts, clashes = interface(a, b, backend)
                pair_counts[name] = {"contacts": len(contacts), "clashes": clashes}
            pocket_errors = {}
            for role, mask in masks.items():
                rotation, translation, protein_error = fit(proteins[role].xyz, native[role].xyz)
                assert protein_error < 1e-5, "Audit requires rigid bound-component controls"
                pocket_errors[role] = rmsd((ligand.xyz @ rotation + translation)[mask], reference["ligand"][mask])
            result = {"candidate": "native" if index is None else rows[index]["candidate"],
                      "pairs": pair_counts, "local_contact_atom_rmsd_A": pocket_errors}
            if index is not None:
                assert sum(p["clashes"] for p in pair_counts.values()) == int(rows[index]["clash_pairs"])
                result["rank"] = int(rows[index]["rank"])
            if dockq_python:
                model_file = folder / "model.cif"
                write_pair(model_file, proteins["poi"], proteins["e3"])
                json_file = folder / "dockq.json"
                subprocess.run([dockq_python, "-m", "DockQ.DockQ", str(model_file), str(native_file),
                                "--mapping", "AB:AB", "--no_align", "--n_cpu", "1",
                                "--json", str(json_file), "--short"],
                               check=True, capture_output=True, text=True, timeout=120)
                data = json.loads(json_file.read_text())
                result["dockq"] = data["best_result"]
                if index is None:
                    assert data["GlobalDockQ"] > 0.99, "Native self-comparison must pass"
            details[label] = result
    clashes = [int(row["clash_pairs"]) for row in rows]
    counts = Counter(clashes)
    errors = [row["e3_ca_rmsd_A"] for row in evaluation["candidates"]]
    return {
        "input_hashes": inputs.provenance["hashes"],
        "file_hashes": {name: digest(run / name) for name in ("ensemble.npz", "candidates.csv", "evaluation.json", "manifest.json")},
        "reference_sha256": digest(reference_path), "candidate_count": len(rows),
        "zero_clash_candidates": sum(c == 0 for c in clashes),
        "candidates_with_unique_clash_count": sum(counts[c] == 1 for c in clashes),
        "correlations_to_e3_error_descriptive_only": {
            key: float(spearmanr([float(row[key]) for row in rows], errors).statistic)
            for key in ("rank", "clash_pairs", "head_rmsd_A", "ppi_contacts", "energy_kcal_mol")},
        "contact_atom_selection": "Native ligand atoms <=4 A from own protein and >5 A from the other; evaluation only",
        "contact_atom_names": {role: list(np.asarray(inputs.atom_names)[mask]) for role, mask in masks.items()},
        "selected_candidates": details,
        "limitation": "DockQ is evaluated only for native, rank 1 and the best existing E3 RMSD candidate, not all candidates."
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-root", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dockq-python")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use a new output file.")
    if len({p.name for p in args.case_root}) != len(args.case_root):
        parser.error("Case directory names must be unique.")
    versions = {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "gemmi", "rdkit")}
    if args.dockq_python:
        versions["DockQ"] = subprocess.check_output(
            [args.dockq_python, "-c", "from importlib.metadata import version; print(version('DockQ'))"],
            text=True, timeout=30).strip()
    with threadpool_limits(limits=1):
        result = {"schema_version": 1, "versions": versions, "script_sha256": digest(__file__),
                  "cases": {case.name: audit(case, args.dockq_python) for case in args.case_root}}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
