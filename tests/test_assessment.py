"""Independent assessment invariants and subprocess failure boundaries."""

import json
from itertools import combinations
import os
from pathlib import Path
import sys

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from topics.assessment import assess, head_groups, identity_order, stereo_valid, success_summary, symmetry_maps
from topics.benchmark import prepare
from topics.engine import digest, ligand_cloud, load_inputs, run
from topics.geometry import fit, rmsd
from topics.inventory import check_inventory
from topics.structures import Cloud, write_complex


def test_complete_identity_mapping_rejects_duplicates_and_missing_atoms():
    expected = [("A", "1", "", "ALA", "CA"), ("A", "1", "A", "ALA", "CA")]
    assert identity_order(expected, expected[::-1]) == [1, 0]
    for invalid in (expected[:1], expected[:1] * 2):
        with pytest.raises(ValueError):
            identity_order(expected, invalid)


def test_symmetry_and_geometry_stereochemistry():
    symmetric = Chem.MolFromSmiles("ClCCCl")
    maps = symmetry_maps(symmetric)
    xyz = np.array([[0., 0., 0.], [1., 0., 0.], [2., 1., 0.], [3., 1., 1.]])
    permuted = xyz[::-1]
    assert rmsd(xyz, permuted) > 0
    assert min(rmsd(xyz, permuted[m]) for m in maps) == 0
    molecule = Chem.AddHs(Chem.MolFromSmiles("F[C@](Cl)(Br)I"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=9) == 0
    coordinates = molecule.GetConformer().GetPositions()
    assert stereo_valid(molecule, coordinates)
    assert not stereo_valid(molecule, coordinates * [-1, 1, 1])
    alkene = Chem.AddHs(Chem.MolFromSmiles("F/C=C/Cl"))
    assert AllChem.EmbedMolecule(alkene, randomSeed=9) == 0
    assert stereo_valid(alkene, alkene.GetConformer().GetPositions())
    opposite = Chem.AddHs(Chem.MolFromSmiles("F/C=C\\Cl"))
    assert not stereo_valid(opposite, alkene.GetConformer().GetPositions())
    with pytest.raises(ValueError, match="budget"):
        symmetry_maps(symmetric, limit=1)


def test_chemical_partitions_and_missing_metric_denominators():
    molecule = Chem.MolFromSmiles("CCCCCCCC")
    names = [str(i) for i in range(8)]
    mapping = {"head_groups": {"poi": names[:3], "linker": names[3:5], "e3": names[5:]},
               "head_definition_source": "Synthetic connectivity control"}
    assert head_groups(molecule, names, mapping)["linker"] == [3, 4]
    mapping["head_groups"]["poi"] = names[:2] + [names[4]]
    mapping["head_groups"]["linker"] = names[2:4]
    with pytest.raises(ValueError, match="Disconnected"):
        head_groups(molecule, names, mapping)
    result = success_summary([{"ppi_success": None}, {"ppi_success": True}])["ppi_success"]
    assert result["candidate_denominator"] == 2 and result["assessed"] == 1
    assert result["success_at_1"] is False and result["oracle_success"] is True


def test_inventory_detects_protected_group_leakage(tmp_path):
    data = json.loads((Path(__file__).resolve().parents[1] / "benchmarks/inventory-v1.json").read_text())
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(data))
    assert check_inventory(path)["status"] == "valid"
    eligible = [c for c in data["cases"] if c["status"] == "eligible"]
    group = next(c["target_family"] for c in eligible if c["split"] == "test")
    one = next(c for c in eligible if c["target_family"] == group)
    one["split"] = "validation"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="crosses"):
        check_inventory(path)


@pytest.mark.skipif(not os.getenv("TOPICS_TEST_DATA"), reason="Requires local PDB/CCD files")
def test_native_assessment_mapping_invariance_and_failure_reporting(tmp_path):
    specification = next((Path(__file__).resolve().parents[1] / "examples").glob("*.json"))
    config = prepare(specification, tmp_path / "inputs", os.environ["TOPICS_TEST_DATA"])
    run_path = run(config, tmp_path / "run", count=1, seeds=[7], pause=0, progress=lambda _: None)
    inputs = load_inputs(config)
    reference_path = config.parent / "reference.npz"
    with np.load(reference_path) as file:
        reference = {k: file[k] for k in file.files}
    # A native pose is a mathematical positive control, never a scored prediction.
    rotation, translation, _ = fit(reference["poi"], inputs.poi.xyz)
    ligand = reference["ligand"] @ rotation + translation
    np.savez_compressed(run_path / "ensemble.npz", poi=inputs.poi.xyz,
                        e3=(reference["e3"] @ rotation + translation)[None], ligand=ligand[None])
    mol = next(iter(Chem.SDMolSupplier(str(run_path / "ensemble.sdf"))))
    for i, xyz in enumerate(ligand):
        mol.GetConformer().SetAtomPosition(i, xyz)
    with Chem.SDWriter(str(run_path / "ensemble.sdf")) as writer:
        writer.write(mol)
    manifest = json.loads((run_path / "manifest.json").read_text())
    for name in manifest["output_hashes"]:
        manifest["output_hashes"][name] = digest(run_path / name)
    (run_path / "manifest.json").write_text(json.dumps(manifest))
    before = {p.name: digest(p) for p in run_path.iterdir() if p.is_file()}
    worker = os.getenv("TOPICS_TEST_EVALUATION_PYTHON")
    first = assess(run_path, config, reference_path, tmp_path / "first", evaluation_python=worker)
    row = first["candidates"][0]
    assert first["status"] == "complete"
    assert row["e3_ca_rmsd_A"] < 1e-5 and row["ligand_symmetry_rmsd_A"] < 1e-5
    assert row["joint_success"] is None and first["head_metric_applicable"] is False
    if worker:
        assert row["dockq"] > 0.999 and row["pb_mol_fast_valid"] is not None
    native_file = tmp_path / "native.cif"
    write_complex(native_file, {**{role.upper(): Cloud(reference[role], getattr(inputs, role).keys, getattr(inputs, role).elements)
                                 for role in ("poi", "e3")}, "PROTAC": ligand_cloud(inputs, reference["ligand"])})
    mapping = {"reference_selection": {"poi": [inputs.poi.keys[0][0]], "e3": [inputs.e3.keys[0][0]],
                "ligand_chain": "PROTAC_L", "ligand": inputs.mol.GetProp("ccd_id"), "ligand_residue": "1"}}
    for cuts in combinations(range(inputs.mol.GetNumBonds()), 2):
        parts = Chem.GetMolFrags(Chem.FragmentOnBonds(inputs.mol, list(cuts), addDummies=False))
        if len(parts) != 3:
            continue
        poi = next((p for p in parts if set(inputs.poi_indices) <= set(p)), None)
        e3 = next((p for p in parts if set(inputs.e3_indices) <= set(p)), None)
        if poi is not None and e3 is not None and poi != e3:
            linker = next(p for p in parts if p not in (poi, e3))
            mapping["head_groups"] = {role: [inputs.atom_names[i] for i in atoms]
                                      for role, atoms in (("poi", poi), ("e3", e3), ("linker", linker))}
            mapping["head_definition_source"] = "Connectivity-only test partition; not a chemical benchmark annotation."
            break
    assert "head_groups" in mapping
    mapping_file = tmp_path / "mapping.json"
    mapping_file.write_text(json.dumps(mapping))
    explicit = assess(run_path, config, native_file, tmp_path / "explicit", mapping_file, worker)
    assert explicit["head_metric_applicable"] is True
    assert all(explicit["candidates"][0][role + "_head_local_rmsd_A"] < 1e-5 for role in ("poi", "e3"))
    if worker:
        assert explicit["candidates"][0]["joint_success"] is True
    for role in ("poi", "e3", "ligand"):
        reference[role] = (reference[role] @ np.diag([-1., -1., 1.]) + [20, -40, 60])[::-1]
        key = "atom_names" if role == "ligand" else role + "_keys"
        reference[key] = reference[key][::-1]
    transformed = tmp_path / "transformed.npz"
    np.savez_compressed(transformed, **reference)
    second = assess(run_path, config, transformed, tmp_path / "second")
    for field in ("e3_ca_rmsd_A", "ligand_symmetry_rmsd_A", "poi_anchor_local_rmsd_A", "e3_anchor_local_rmsd_A"):
        assert second["candidates"][0][field] == pytest.approx(row[field], abs=1e-5)
    # Empty native interfaces are an explicit separate stratum, not a division by zero.
    reference["e3"] += [1000, 0, 0]
    np.savez_compressed(transformed, **reference)
    empty = assess(run_path, config, transformed, tmp_path / "empty", evaluation_python=worker)
    assert empty["native_ppi_contacts"] == 0
    assert empty["candidates"][0]["native_contact_recall"] is None
    assert empty["candidates"][0]["ppi_success"] is None
    assert empty["candidates"][0]["no_ppi_fallback_success"] is False
    # A missing optional environment must retain failed candidates and exit as partial.
    failed = assess(run_path, config, reference_path, tmp_path / "failed", evaluation_python=str(tmp_path / "absent-python"))
    assert failed["status"] == "partial" and failed["candidate_count"] == 1
    assert failed["candidates"][0]["errors"] == ["external_evaluation_failed"]
    assert before == {p.name: digest(p) for p in run_path.iterdir() if p.is_file()}
    report = (tmp_path / "first/report.html").read_text()
    assert "plotly-graph-div" in report and "<script src=" not in report
    with pytest.raises(ValueError, match="new output"):
        assess(run_path, config, reference_path, tmp_path / "first")
    (run_path / "ensemble.sdf").write_text("corrupted")
    with pytest.raises(ValueError, match="hash mismatch"):
        assess(run_path, config, reference_path, tmp_path / "corrupt")
