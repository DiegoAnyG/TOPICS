"""Read-only, reference-separated assessment of a saved candidate ensemble."""

import csv
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time

import gemmi
import numpy as np
from rdkit import Chem
from threadpoolctl import threadpool_limits

from .engine import digest, ligand_cloud, load_inputs
from .geometry import Distances, fit, interface, rmsd
from .structures import Cloud, read_cloud, write_complex


def identity_order(expected, observed):
    """Map complete identities, allowing permutations but never positional guesses."""
    expected, observed = list(map(tuple, expected)), list(map(tuple, observed))
    if len(set(expected)) != len(expected) or len(set(observed)) != len(observed):
        raise ValueError("Duplicate atom identities in prediction or reference.")
    if set(expected) != set(observed):
        raise ValueError("Atom identities differ; supply an explicit complete identity mapping.")
    lookup = {key: i for i, key in enumerate(observed)}
    return [lookup[key] for key in expected]


def reference_coordinates(path, inputs, mapping):
    """NPZ uses saved identities; mmCIF requires explicit author-chain selections."""
    names = inputs.atom_names
    if Path(path).suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            data = {key: archive[key] for key in archive.files}
    else:
        selection = mapping.get("reference_selection")
        if not selection:
            raise ValueError("A coordinate reference requires reference_selection in --mapping JSON.")
        data = {}
        for role in ("poi", "e3"):
            cloud = read_cloud(path, selection[role])
            data[role], data[role + "_keys"] = cloud.xyz, cloud.keys
        cloud = read_cloud(path, [selection["ligand_chain"]], selection["ligand"], selection["ligand_residue"])
        data["ligand"], data["atom_names"] = cloud.xyz, [key[-1] for key in cloud.keys]
    result = {}
    for role in ("poi", "e3"):
        expected = [[str(v) for v in key] for key in getattr(inputs, role).keys]
        observed = [[str(v) for v in key] for key in data[role + "_keys"]]
        pairs = mapping.get("protein_atom_map", {}).get(role)
        if pairs is not None:
            source = [pair["model"] for pair in pairs]
            order = identity_order(expected, source)
            expected = [pairs[i]["reference"] for i in order]
        result[role] = np.asarray(data[role])[identity_order(expected, observed)]
    ligand_map = mapping.get("ligand_atom_map", {n: n for n in names})
    if set(ligand_map) != set(names):
        raise ValueError("ligand_atom_map must cover every model heavy atom exactly once.")
    order = identity_order([(ligand_map[n],) for n in names], [(str(n),) for n in data["atom_names"]])
    result["ligand"] = np.asarray(data["ligand"])[order]
    for role, xyz in result.items():
        n = len(names) if role == "ligand" else len(getattr(inputs, role).xyz)
        if xyz.shape != (n, 3) or not np.isfinite(xyz).all():
            raise ValueError(f"Invalid or incomplete reference {role} coordinates.")
    return result


def head_groups(mol, names, mapping):
    """Validate user-curated chemical partitions; anchor rings are not whole heads."""
    groups = mapping.get("head_groups")
    if groups is None:
        return None
    if set(groups) != {"poi", "e3", "linker"} or not mapping.get("head_definition_source"):
        raise ValueError("head_groups needs poi/e3/linker and a head_definition_source citation.")
    flattened = [n for atoms in groups.values() for n in atoms]
    if len(flattened) != len(set(flattened)) or set(flattened) != set(names):
        raise ValueError("Head/linker groups must partition every ligand heavy atom without overlap.")
    indices = {role: [names.index(n) for n in atoms] for role, atoms in groups.items()}
    for role, atoms in indices.items():
        if len(atoms) < (3 if role != "linker" else 1):
            raise ValueError("Chemical heads need at least three atoms; the linker must be nonempty.")
        reached, pending = set(), [atoms[0]]
        while pending:
            i = pending.pop()
            if i in reached:
                continue
            reached.add(i)
            pending.extend(a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors()
                           if a.GetIdx() in atoms and a.GetIdx() not in reached)
        if reached != set(atoms):
            raise ValueError(f"Disconnected chemical group: {role}.")
    return indices


def symmetry_maps(mol, limit=4096):
    matches = mol.GetSubstructMatches(mol, uniquify=False, useChirality=True, maxMatches=limit + 1)
    if not matches or len(matches) > limit:
        raise ValueError("Exact symmetry enumeration exceeds its budget or has no valid mapping.")
    return np.asarray(matches, dtype=int)


def stereo_valid(mol, xyz):
    """Compare defined tetrahedral and double-bond stereo to geometry, not SDF tags."""
    original = Chem.Mol(mol)
    Chem.AssignStereochemistry(original, cleanIt=True, force=True)
    posed = Chem.Mol(original)
    posed.RemoveAllConformers()
    conformer = Chem.Conformer(posed.GetNumAtoms())
    for i, point in enumerate(xyz):
        conformer.SetAtomPosition(i, point)
    posed.AddConformer(conformer)
    Chem.RemoveStereochemistry(posed)
    Chem.AssignStereochemistryFrom3D(posed, replaceExistingTags=True)
    Chem.AssignStereochemistry(posed, cleanIt=True, force=True)
    for a, b in zip(original.GetAtoms(), posed.GetAtoms()):
        if a.HasProp("_CIPCode") and (not b.HasProp("_CIPCode") or a.GetProp("_CIPCode") != b.GetProp("_CIPCode")):
            return False
    for a, b in zip(original.GetBonds(), posed.GetBonds()):
        if a.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ) and a.GetStereo() != b.GetStereo():
            return False
    return True


def pair_metrics(a, b, backend):
    contacts, legacy = interface(a, b, backend)
    table = Chem.GetPeriodicTable()
    radii_a = np.asarray([table.GetRvdw(e) for e in a.elements])
    radii_b = np.asarray([table.GetRvdw(e) for e in b.elements])
    typed = 0
    for start, distances in backend.blocks(a.xyz, b.xyz):
        bounds = 0.75 * (radii_a[start:start + len(distances), None] + radii_b[None, :])
        typed += int(np.count_nonzero(distances < bounds))
    return {"contacts_5A": len(contacts), "clashes_2A": legacy, "vdw_clashes": typed}


def write_pair(path, poi, e3):
    if any(len({k[0] for k in cloud.keys}) != 1 for cloud in (poi, e3)):
        raise ValueError("DockQ assessment currently requires one explicitly selected chain per protein role.")
    write_complex(path, {"POI": poi, "E3": e3})
    structure = gemmi.read_structure(str(path))
    for chain, name in zip(structure[0], ("A", "B")):
        chain.name = name
        # Identity matching has already happened. Canonical numbers keep DockQ's
        # numbering alignment unambiguous for insertion codes and author offsets.
        for number, residue in enumerate(chain, 1):
            residue.seqid = gemmi.SeqId(number, " ")
    structure.make_mmcif_document().write_file(str(path))


def success_summary(rows):
    result = {}
    for metric in ("ppi_success", "ligand_success", "chemical_valid", "steric_valid", "joint_success"):
        values = [row.get(metric) for row in rows]
        available = sum(v is not None for v in values)
        result[metric] = {"assessed": available, "passed": sum(v is True for v in values),
                          "candidate_denominator": len(rows),
                          **{f"success_at_{k}": any(v is True for v in values[:k]) if available else None
                             for k in (1, 5, 10)},
                          "oracle_success": any(v is True for v in values) if available else None}
    return result


def assess(run_dir, config_path, reference_path, output, mapping_path=None, evaluation_python=None, timeout=1800):
    run_dir, output = Path(run_dir), Path(output)
    if output.exists():
        raise ValueError("Assessment requires a new output directory.")
    if output.resolve().is_relative_to(run_dir.resolve()):
        raise ValueError("Save assessment outside the original run directory.")
    if not 1 <= timeout <= 86400:
        raise ValueError("Evaluation timeout must be between 1 and 86400 seconds.")
    manifest = json.loads((run_dir / "manifest.json").read_text())
    if manifest.get("status") != "complete":
        raise ValueError("Only complete saved ensembles can be assessed; retain failed runs in the benchmark denominator.")
    inputs = load_inputs(config_path)
    if inputs.provenance["hashes"] != manifest["provenance"]["hashes"]:
        raise ValueError("Input content hashes differ from the recorded prediction inputs.")
    hashes = {}
    for name in ("candidates.csv", "ensemble.npz", "ensemble.sdf", "best.cif"):
        hashes[name] = digest(run_dir / name)
        if hashes[name] != manifest["output_hashes"].get(name):
            raise ValueError(f"Saved output hash mismatch: {name}.")
    for role in ("poi", "e3"):
        if [list(k) for k in getattr(inputs, role).keys] != manifest["atom_keys"][role]:
            raise ValueError(f"Input {role} identities differ from the ensemble.")
    if inputs.atom_names != manifest["atom_keys"]["ligand"]:
        raise ValueError("Ligand atom order differs from the ensemble.")
    with (run_dir / "candidates.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows or len({row["candidate"] for row in rows}) != len(rows):
        raise ValueError("Missing or duplicate candidate identifiers.")
    if [int(row["rank"]) for row in rows] != list(range(1, len(rows) + 1)):
        raise ValueError("Saved candidates must have consecutive ranks in ensemble order.")
    with np.load(run_dir / "ensemble.npz", allow_pickle=False) as archive:
        arrays = {role: archive[role] for role in ("poi", "e3", "ligand")}
    for role, xyz in arrays.items():
        n = len(inputs.atom_names) if role == "ligand" else len(getattr(inputs, role).xyz)
        shape = (n, 3) if role == "poi" else (len(rows), n, 3)
        if xyz.shape != shape or not np.isfinite(xyz).all():
            raise ValueError(f"Invalid saved {role} coordinates or candidate count.")
    if not np.array_equal(arrays["poi"], inputs.poi.xyz):
        raise ValueError("Saved POI differs from its hashed input coordinates.")
    molecules = list(Chem.SDMolSupplier(str(run_dir / "ensemble.sdf"), removeHs=True))
    if len(molecules) != len(rows):
        raise ValueError("SDF and coordinate ensemble counts differ.")
    def ordered_graph(molecule):
        return ([ (a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsotope()) for a in molecule.GetAtoms()],
                sorted((min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                        str(b.GetBondType())) for b in molecule.GetBonds()))
    graph = ordered_graph(inputs.mol)
    for row, mol, xyz in zip(rows, molecules, arrays["ligand"]):
        if (mol is None or not mol.HasProp("_Name") or mol.GetProp("_Name") != row["candidate"]
                or ordered_graph(mol) != graph
                or mol.GetNumAtoms() != len(inputs.atom_names)
                or not np.allclose(mol.GetConformer().GetPositions(), xyz, atol=0.000051, rtol=0)):
            raise ValueError("SDF graph, identity or coordinates differ from the saved ensemble.")
    mapping = json.loads(Path(mapping_path).read_text()) if mapping_path else {}
    reference = reference_coordinates(reference_path, inputs, mapping)
    groups = head_groups(inputs.mol, inputs.atom_names, mapping)
    maps = symmetry_maps(inputs.mol)
    native = {role: Cloud(reference[role], getattr(inputs, role).keys, getattr(inputs, role).elements)
              for role in ("poi", "e3")}
    native_ligand = ligand_cloud(inputs, reference["ligand"])
    backend = Distances("cpu")
    native_contacts, _ = interface(native["poi"], native["e3"], backend)
    ca = np.array([k[-1] == "CA" for k in inputs.e3.keys])
    if not ca.any():
        raise ValueError("The selected E3 has no alpha carbon atoms.")
    output.mkdir(parents=True)
    result = {"schema_version": 1, "status": "running", "protocol": "topics-assessment-v1",
              "reference_used_for_ranking": False, "candidate_count": len(rows),
              "input_hashes": inputs.provenance["hashes"], "prediction_hashes": hashes,
              "prediction_manifest_sha256": digest(run_dir / "manifest.json"),
              "reference_sha256": digest(reference_path),
              "mapping_sha256": digest(mapping_path) if mapping_path else None,
              "mapping": mapping, "versions": {p: importlib.metadata.version(p) for p in ("numpy", "rdkit", "gemmi", "scipy")},
              "implementation_sha256": {p.name: digest(p) for p in Path(__file__).parent.glob("*.py")},
              "generation_statistics": manifest["conformers"], "prediction_parameters": manifest["parameters"],
              "native_ppi_contacts": len(native_contacts), "native_chemistry_stereo_valid": stereo_valid(inputs.mol, reference["ligand"]),
              "symmetry_mappings": len(maps), "head_metric_applicable": groups is not None,
              "criteria": {"ppi_dockq_minimum": 0.23, "ligand_symmetry_rmsd_maximum_A": 2.0,
                           "vdw_distance_fraction": 0.75, "head_local_rmsd_maximum_A": 2.0,
                           "no_native_ppi_fallback": "E3 CA RMSD <= 5 A after fitting POI; exploratory, separately stratified"},
              "limitations": ["Bound-component controls do not measure blind prediction accuracy.",
                              "mol_fast is a PoseBusters ligand-geometry subset, not full PB-valid.",
                              "Sterics exclude hydrogens, cofactors, solvent and covalent cross-component bonds.",
                              "Joint success requires curated complete heads; absent definitions produce null, not a pass.",
                              "All atoms must map explicitly; missing residues require curated matched structures.",
                              "Success@k uses saved candidates, not yet distinct structural clusters."]}
    target = output / "assessment.json"
    target.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    start = time.perf_counter()
    evaluated = []
    try:
        with threadpool_limits(limits=1), tempfile.TemporaryDirectory(prefix="topics-assessment-") as folder:
            folder = Path(folder)
            job = {"native": str(folder / "native.cif"), "poses": [], "sdf": str(run_dir.resolve() / "ensemble.sdf"),
                   "native_sdf": str(folder / "native.sdf"), "ppi_applicable": bool(native_contacts)}
            worker_python = evaluation_python or os.environ.get("TOPICS_EVALUATION_PYTHON")
            if worker_python:
                write_pair(folder / "native.cif", native["poi"], native["e3"])
                mol = Chem.Mol(inputs.mol)
                mol.RemoveAllConformers()
                conformer = Chem.Conformer(mol.GetNumAtoms())
                for i, xyz in enumerate(reference["ligand"]):
                    conformer.SetAtomPosition(i, xyz)
                mol.AddConformer(conformer)
                with Chem.SDWriter(str(folder / "native.sdf")) as writer:
                    writer.write(mol)
            result["native_pairs"] = {name: pair_metrics(a, b, backend) for name, a, b in (
                ("ppi", native["poi"], native["e3"]), ("poi_ligand", native["poi"], native_ligand),
                ("e3_ligand", native["e3"], native_ligand))}
            rotation, translation, _ = fit(inputs.poi.xyz, native["poi"].xyz)
            for index, row in enumerate(rows):
                proteins = {"poi": inputs.poi, "e3": Cloud(arrays["e3"][index], inputs.e3.keys, inputs.e3.elements)}
                ligand = ligand_cloud(inputs, arrays["ligand"][index])
                moved_e3 = proteins["e3"].moved(rotation, translation)
                moved_ligand = ligand.xyz @ rotation + translation
                contacts, _ = interface(native["poi"], moved_e3, backend)
                entry = {"candidate": row["candidate"], "rank": int(row["rank"]), "seed": int(row["seed"]),
                         "e3_ca_rmsd_A": rmsd(moved_e3.xyz[ca], native["e3"].xyz[ca]),
                         "ligand_rmsd_A": rmsd(moved_ligand, reference["ligand"]),
                         "ligand_symmetry_rmsd_A": min(rmsd(moved_ligand, reference["ligand"][m]) for m in maps),
                         "native_contact_recall": len(contacts & native_contacts) / len(native_contacts) if native_contacts else None,
                         "contact_precision": len(contacts & native_contacts) / len(contacts) if contacts else 0.0,
                         "stereo_valid": stereo_valid(inputs.mol, ligand.xyz), "errors": []}
                for name, a, b in (("ppi", proteins["poi"], proteins["e3"]),
                                   ("poi_ligand", proteins["poi"], ligand), ("e3_ligand", proteins["e3"], ligand)):
                    entry.update({f"{name}_{key}": value for key, value in pair_metrics(a, b, backend).items()})
                for role in ("poi", "e3"):
                    r, t, _ = fit(proteins[role].xyz, native[role].xyz)
                    atoms = getattr(inputs, role + "_indices")
                    entry[role + "_anchor_local_rmsd_A"] = rmsd((ligand.xyz @ r + t)[atoms], reference["ligand"][atoms])
                    entry[role + "_head_local_rmsd_A"] = (rmsd((ligand.xyz @ r + t)[groups[role]], reference["ligand"][groups[role]])
                                                            if groups else None)
                if worker_python:
                    filename = folder / f"candidate_{index}.cif"
                    write_pair(filename, proteins["poi"], proteins["e3"])
                    job["poses"].append({"candidate": row["candidate"], "path": str(filename)})
                evaluated.append(entry)
            result["external_tools"] = {"status": "not_requested"}
            if worker_python:
                job_file, response = folder / "job.json", folder / "result.json"
                job_file.write_text(json.dumps(job))
                env = dict(os.environ)
                env.update({name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")})
                try:
                    completed = subprocess.run([str(worker_python), str(Path(__file__).with_name("assessment_worker.py")),
                                                str(job_file), str(response)], capture_output=True, text=True,
                                               timeout=timeout, env=env, check=False)
                    if completed.returncode != 0 or not response.exists():
                        raise RuntimeError("Evaluation worker failed; verify its isolated dependencies.")
                    external = json.loads(response.read_text())
                    if [e["candidate"] for e in external["candidates"]] != [e["candidate"] for e in evaluated]:
                        raise ValueError("Evaluation worker returned inconsistent candidate identities.")
                    result["external_tools"] = {k: v for k, v in external.items() if k != "candidates"}
                    for entry, metrics in zip(evaluated, external["candidates"]):
                        entry["errors"].extend(metrics.pop("errors"))
                        entry.update(metrics)
                except (OSError, subprocess.TimeoutExpired, RuntimeError, ValueError, KeyError):
                    result["external_tools"] = {"status": "failed", "reason": "Worker failed, timed out or returned invalid data."}
                    for entry in evaluated:
                        entry["errors"].append("external_evaluation_failed")
            for entry in evaluated:
                dockq = entry.get("dockq")
                entry["ppi_success"] = dockq >= 0.23 if dockq is not None else None
                entry["ligand_success"] = entry["ligand_symmetry_rmsd_A"] <= 2.0
                entry["chemical_valid"] = (entry["stereo_valid"] and entry["pb_mol_fast_valid"]
                                           if entry.get("pb_mol_fast_valid") is not None else None)
                entry["steric_valid"] = all(entry[name + "_vdw_clashes"] == 0 for name in ("ppi", "poi_ligand", "e3_ligand"))
                geometry = entry["ppi_success"] if native_contacts else entry["e3_ca_rmsd_A"] <= 5.0
                entry["no_ppi_fallback_success"] = None if native_contacts else geometry
                entry["joint_success"] = (bool(geometry and entry["ligand_success"] and entry["chemical_valid"] and entry["steric_valid"]
                                                and all(entry[r + "_head_local_rmsd_A"] <= 2 for r in ("poi", "e3")))
                                          if geometry is not None and entry["chemical_valid"] is not None and groups else None)
            result.update(status="complete" if all(not e["errors"] for e in evaluated) else "partial",
                          candidates=evaluated, summary=success_summary(evaluated),
                          compute_seconds=time.perf_counter() - start)
            target.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
            with (output / "assessment.csv").open("w", newline="") as stream:
                flattened = [{k: json.dumps(v) if isinstance(v, (list, dict)) else v for k, v in e.items()} for e in evaluated]
                writer = csv.DictWriter(stream, fieldnames=sorted(set().union(*(e.keys() for e in flattened))))
                writer.writeheader(); writer.writerows(flattened)
            from .report import render_assessment
            render_assessment(output, result)
    except Exception:
        result["status"] = "failed"
        result["failure"] = "Assessment interrupted; incomplete output must not count as a successful case."
        target.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        raise
    return result
