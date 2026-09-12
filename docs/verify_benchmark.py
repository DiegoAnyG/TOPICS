"""Recheck frozen inventory evidence using local PDB mmCIF archives and CCD files.

Usage: python docs/verify_benchmark.py --manifest benchmarks/inventory-v1.json
       --data-dir <PDB-and-CCD-downloads> --out <new-verification.json>
No predictions, ranking, parameter tuning or network access occur here.
"""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import tempfile

import gemmi
from rdkit.Chem.Scaffolds import MurckoScaffold
from threadpoolctl import threadpool_limits

from topics.chemistry import read_ccd
from topics.engine import digest
from topics.geometry import Distances, interface
from topics.inventory import check_inventory
from topics.structures import read_cloud


def verify(manifest_path, data_dir):
    summary = check_inventory(manifest_path)
    manifest = json.loads(Path(manifest_path).read_text())
    records = []
    with threadpool_limits(limits=1), tempfile.TemporaryDirectory(prefix="topics-inventory-") as folder:
        for case in manifest["cases"]:
            code = case["pdb_id"]
            source = data_dir / f"{code}.cif"
            if source.exists():
                raw = source.read_bytes()
            else:
                raw = gzip.decompress((data_dir / f"pdb_0000{code.lower()}.cif.gz").read_bytes())
            if hashlib.sha256(raw).hexdigest() != case["coordinate_sha256"]:
                raise ValueError(f"Coordinate snapshot changed for {code}; do not overwrite the frozen inventory.")
            checked = {"pdb_id": code, "coordinate_hash_verified": True, "status": case["status"]}
            if case["status"] == "eligible":
                ccd = data_dir / (case["ligand"] + ".cif")
                if digest(ccd) != case["ccd_sha256"]:
                    raise ValueError(f"CCD snapshot changed for {code}.")
                molecule, names = read_ccd(ccd)
                scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule, includeChirality=False)
                if hashlib.sha256(scaffold.encode()).hexdigest() != case["whole_ligand_scaffold_sha256"]:
                    raise ValueError(f"Scaffold identity changed for {code}.")
                path = Path(folder) / "reference.cif"
                path.write_bytes(raw)
                selection = case["selection"]
                ligand = read_cloud(path, [selection["ligand_chain"]], case["ligand"], selection["ligand_residue"])
                if {key[-1] for key in ligand.keys} != set(names):
                    raise ValueError(f"Heavy-atom completeness failed for {code}.")
                poi, e3 = (read_cloud(path, [selection[role + "_chain"]]) for role in ("poi", "e3"))
                contacts, _ = interface(poi, e3, Distances("cpu"))
                if len(contacts) != case["native_ppi_contacts_5A"]:
                    raise ValueError(f"Native interface selection changed for {code}.")
                block = gemmi.cif.read_string(raw.decode()).sole_block()
                operators = block.get_mmcif_category("_pdbx_struct_assembly_gen.")
                atom = block.get_mmcif_category("_atom_site.")
                labels = {label for chain, label, comp, number in zip(atom["auth_asym_id"], atom["label_asym_id"],
                          atom["label_comp_id"], atom["auth_seq_id"])
                          if (chain in (selection["poi_chain"], selection["e3_chain"])
                              and gemmi.find_tabulated_residue(comp).is_amino_acid())
                          or (chain == selection["ligand_chain"] and comp == case["ligand"] and number == selection["ligand_residue"])}
                matches = [identifier for identifier, operation, chains in zip(operators["assembly_id"], operators["oper_expression"], operators["asym_id_list"])
                           if operation in ("1", "(1)") and labels <= set(chains.split(","))]
                if matches != case["assembly_ids_with_identity_operation"] or not matches:
                    raise ValueError(f"Biological assembly selection changed for {code}.")
                checked.update(ccd_hash_verified=True, complete_heavy_atoms=len(names), native_ppi_contacts=len(contacts),
                               assembly_verified=True, scaffold_verified=True)
            records.append(checked)
    return {"manifest_sha256": digest(manifest_path), "verifier_sha256": digest(__file__),
            "summary": summary, "cases": records}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("Use a new verification output file.")
    result = verify(args.manifest, args.data_dir)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result["summary"], indent=2))
