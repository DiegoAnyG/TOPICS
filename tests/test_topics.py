"""Checks scientific invariants and the application boundaries, not preferred scores."""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from topics.benchmark import evaluate, prepare
from topics.chemistry import generate
from topics.engine import Inputs, assemble, load_inputs, rank_candidates, run
from topics.geometry import Distances, fit, interface, rmsd
from topics.structures import Cloud, read_cloud, write_complex


def cloud(xyz, chain="A"):
    return Cloud(np.array(xyz, dtype=float), [(chain, i + 1, "", "ALA", "CA") for i in range(len(xyz))], ["C"] * len(xyz))


def rotation():
    return np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])


def test_rigid_fit_and_reflection():
    a = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
    b = a @ rotation() + [5., -3., 8.]
    r, t, error = fit(a, b)
    assert error < 1e-10
    assert np.linalg.det(r) == pytest.approx(1.)
    mirrored = a * [-1, 1, 1]
    assert fit(a, mirrored)[2] > 0.1


def test_fit_rejects_bad_geometry():
    with pytest.raises(ValueError, match="collinear"):
        fit(np.zeros((3, 3)), np.zeros((3, 3)))
    with pytest.raises(ValueError, match="finite"):
        fit(np.full((3, 3), np.nan), np.eye(3))


def test_independent_binary_frames_do_not_change_assembly():
    ligand = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [8., 0., 0.], [8., 1., 0.], [8., 0., 1.]])
    data = Inputs(cloud(ligand[:3]), cloud(ligand[3:]), None, [], [0, 1, 2], [3, 4, 5], ligand[:3], ligand[3:], "test", {})
    original, pose, _ = assemble(data, ligand)
    r, t = rotation(), np.array([80., -40., 20.])
    data.e3 = data.e3.moved(r, t)
    data.e3_head = data.e3_head @ r + t
    transformed, new_pose, _ = assemble(data, ligand)
    assert rmsd(original.xyz, transformed.xyz) < 1e-10
    np.testing.assert_allclose(pose, new_pose, atol=1e-10)


def test_interfaces_preserve_chain_identity_and_detect_negative_control():
    a, b = cloud([[0, 0, 0], [2, 0, 0]], "A"), cloud([[0, 0, 3]], "B")
    contacts, clashes = interface(a, b, Distances("cpu"))
    assert len(contacts) == 2 and clashes == 0
    assert all(x[0] == "A" and y[0] == "B" for x, y in contacts)
    assert interface(a, b.moved(np.eye(3), [100, 0, 0]), Distances("cpu")) == (set(), 0)


def test_mmcif_roundtrip(tmp_path):
    original = cloud([[0, 0, 0], [2, 1, 0]])
    path = tmp_path / "structure.cif"
    write_complex(path, {"POI": original})
    read = read_cloud(path)
    np.testing.assert_allclose(read.xyz, original.xyz)
    assert read.keys[0][0] == "POI_A"
    with pytest.raises(ValueError, match="Unknown author chain"):
        read_cloud(path, ["missing"])


def test_conformer_repeatability_and_stereochemistry():
    molecule = Chem.MolFromSmiles("C[C@H](O)CCNC(=O)C")
    a, stats = generate(molecule, 3, 17)
    b, other = generate(molecule, 3, 17)
    assert stats == other and len(a) > 0
    for first, second in zip(a, b):
        np.testing.assert_allclose(first[1], second[1], atol=1e-10)
        assert first[2] == pytest.approx(second[2], abs=1e-10)
    assert molecule.GetNumConformers() == 0


def test_rank_ignores_reference_error():
    rows = [{"candidate": "a", "clash_pairs": 2, "head_rmsd_A": 0.1, "ppi_contacts": 5, "energy_kcal_mol": 1, "native_rmsd": 0},
            {"candidate": "b", "clash_pairs": 0, "head_rmsd_A": 0.2, "ppi_contacts": 3, "energy_kcal_mol": 2, "native_rmsd": 99}]
    rank_candidates(rows)
    assert rows[0]["candidate"] == "b"


def test_cli_menu_and_error_exit():
    menu = subprocess.run([sys.executable, "-m", "topics"], input="4\n0\n", text=True, capture_output=True)
    assert menu.returncode == 0 and "Available CPU cores" in menu.stdout
    bad = subprocess.run([sys.executable, "-m", "topics", "assemble", "missing.json", "--out", "unused"], text=True, capture_output=True)
    assert bad.returncode == 1 and "Traceback" not in bad.stderr


def test_cuda_distance_parity():
    try:
        backend = Distances("cuda")
    except ValueError:
        pytest.skip("CUDA runtime/device unavailable")
    rng = np.random.default_rng(7)
    a, b = rng.normal(size=(150, 3)), rng.normal(size=(45, 3))
    cpu = np.concatenate([d for _, d in Distances("cpu").blocks(a, b)])
    gpu = np.concatenate([d for _, d in backend.blocks(a, b)])
    np.testing.assert_allclose(cpu, gpu, rtol=1e-12, atol=1e-12)


@pytest.mark.skipif(not os.getenv("TOPICS_TEST_DATA"), reason="Set TOPICS_TEST_DATA to local PDB/CCD files")
def test_crystal_geometry_controls(tmp_path):
    examples = Path(__file__).resolve().parents[1] / "examples"
    for spec in sorted(examples.glob("*.json")):
        folder = tmp_path / spec.stem
        config = prepare(spec, folder, os.environ["TOPICS_TEST_DATA"])
        positive = json.loads((folder / "positive_control.json").read_text())
        assert positive["e3_rmsd_A"] < 0.001
        assert positive["ligand_rmsd_A"] < 0.001
        loaded = load_inputs(config)
        assert not set(loaded.poi_indices) & set(loaded.e3_indices)
        assert "reference" not in json.loads(config.read_text())


@pytest.mark.skipif(not os.getenv("TOPICS_TEST_DATA"), reason="Requires local PDB/CCD control data")
def test_end_to_end_report_and_reference_validation(tmp_path):
    spec = next((Path(__file__).resolve().parents[1] / "examples").glob("*.json"))
    config = prepare(spec, tmp_path / "inputs", os.environ["TOPICS_TEST_DATA"])
    output = run(config, tmp_path / "run", count=2, seeds=[7], pause=0, progress=lambda _: None)
    summary = evaluate(output, tmp_path / "inputs/reference.npz")
    assert summary["reference_used_for_ranking"] is False
    assert json.loads((output / "manifest.json").read_text())["status"] == "complete"
    report = (output / "report.html").read_text()
    assert "plotly-graph-div" in report and "<script src=" not in report
    for name in ("assessment.svg", "assessment.pdf", "assessment.png", "assembly.svg", "best.cif", "ensemble.sdf"):
        assert (output / name).stat().st_size > 100
    with pytest.raises(ValueError, match="not empty"):
        run(config, output, count=2)
    manifest = json.loads((output / "manifest.json").read_text())
    manifest["status"] = "running"
    (output / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="completed"):
        evaluate(output, tmp_path / "inputs/reference.npz")
