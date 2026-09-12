"""Validate a frozen benchmark manifest without exposing references to prediction."""

from collections import Counter
import json
from pathlib import Path


def check_inventory(path):
    data = json.loads(Path(path).read_text())
    if data.get("schema_version") != 1 or not data.get("protocol", {}).get("frozen_on"):
        raise ValueError("Expected a version-1 inventory with a frozen protocol.")
    cases = data["cases"]
    if not cases or len({c["pdb_id"] for c in cases}) != len(cases):
        raise ValueError("Inventory must contain unique PDB systems, not repeated crystal copies.")
    groups = {}
    for case in cases:
        if case["status"] not in {"eligible", "excluded"}:
            raise ValueError("Each inventory case must be explicitly eligible or excluded.")
        if case["status"] == "excluded" and not case.get("reasons"):
            raise ValueError("Every excluded case requires an explicit reason.")
        if case["status"] == "eligible":
            if case.get("reasons") or not case.get("selection") or not case.get("ccd_sha256"):
                raise ValueError("Eligible cases require complete selections, CCD provenance and no exclusion reason.")
            if case["split"] not in {"development", "validation", "test"}:
                raise ValueError("Invalid dataset split.")
            if case.get("prediction_evaluated_by_TOPICS") and case["split"] != "development":
                raise ValueError("Previously inspected TOPICS predictions must stay in development.")
        if case.get("split"):
            for field in ("target_family", "whole_ligand_scaffold_sha256"):
                if case.get(field):
                    key = (field, case[field])
                    if key in groups and groups[key] != case["split"]:
                        raise ValueError(f"Protected {field} group crosses dataset splits.")
                    groups[key] = case["split"]
    counts = {"inventory": len(cases), "eligible": sum(c["status"] == "eligible" for c in cases),
              "excluded": sum(c["status"] == "excluded" for c in cases),
              "eligible_by_split": dict(Counter(c["split"] for c in cases if c["status"] == "eligible"))}
    if counts != data["counts"]:
        raise ValueError("Recorded inventory counts disagree with case records.")
    return {"status": "valid", "protocol": data["protocol"]["id"], **counts,
            "readiness": data["protocol"]["readiness"]}
