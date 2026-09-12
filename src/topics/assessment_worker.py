"""Optional isolated DockQ/PoseBusters worker; never imported by prediction code."""

import importlib.metadata
import json
from pathlib import Path
import sys


def main(job_file, output):
    from DockQ import DockQ as dockq_module
    from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
    from posebusters import PoseBusters
    from rdkit import Chem
    from threadpoolctl import threadpool_limits

    job = json.loads(Path(job_file).read_text())
    native = load_PDB(job["native"])
    molecules = list(Chem.SDMolSupplier(job["sdf"], removeHs=True))
    if len(molecules) != len(job["poses"]):
        raise ValueError("SDF and candidate counts differ.")
    buster = PoseBusters(config="mol_fast", max_workers=0)

    def clear_pose_caches():
        # Upstream caches retain model structures; release them between candidates.
        for function in vars(dockq_module).values():
            if getattr(function, "__module__", None) == dockq_module.__name__ and hasattr(function, "cache_clear"):
                function.cache_clear()

    def geometry(molecule):
        frame = buster.bust(molecule, full_report=False)
        if len(frame) != 1 or frame.empty or frame.isna().any().any():
            raise ValueError("Missing PoseBusters checks.")
        checks = {str(k): bool(v) for k, v in frame.iloc[0].items()}
        return {"pb_mol_fast_valid": all(checks.values()), "pb_checks": checks}

    def score(path):
        if not job["ppi_applicable"]:
            return {"dockq": None, "dockq_status": "no_native_ppi_contacts"}
        interfaces, _ = run_on_all_native_interfaces(load_PDB(path), native,
                                                    chain_map={"A": "A", "B": "B"}, no_align=True, low_memory=False)
        if len(interfaces) != 1:
            raise ValueError("Expected one explicitly mapped native interface.")
        metric = next(iter(interfaces.values()))
        return {"dockq": float(metric["DockQ"]), "dockq_irmsd_A": float(metric["iRMSD"]),
                "dockq_lrmsd_A": float(metric["LRMSD"]), "dockq_fnat": float(metric["fnat"]), "dockq_status": "assessed"}

    result = {"status": "complete", "versions": {p: importlib.metadata.version(p) for p in
               ("DockQ", "posebusters", "rdkit", "numpy", "biopython")}, "threads": 1,
              "mapping": {"native_A": "model_A (POI)", "native_B": "model_B (E3)"},
              "posebusters_config": "mol_fast", "candidates": []}
    with threadpool_limits(limits=1):
        result["native_control"] = {**score(job["native"]), **geometry(next(iter(Chem.SDMolSupplier(job["native_sdf"]))))}
        if job["ppi_applicable"] and result["native_control"]["dockq"] < 0.999:
            raise ValueError("Native DockQ identity control failed.")
        clear_pose_caches()
        for pose, molecule in zip(job["poses"], molecules):
            entry = {"candidate": pose["candidate"], "errors": []}
            try:
                if molecule is None or molecule.GetProp("_Name") != pose["candidate"]:
                    raise ValueError("Invalid SDF molecule or candidate identity.")
                entry.update(geometry(molecule))
            except Exception as exc:
                entry["pb_mol_fast_valid"] = None
                entry["errors"].append("posebusters_failed:" + type(exc).__name__)
            try:
                entry.update(score(pose["path"]))
            except Exception as exc:
                entry["dockq"] = None
                entry["errors"].append("dockq_failed:" + type(exc).__name__)
            finally:
                clear_pose_caches()
            result["candidates"].append(entry)
    Path(output).write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
