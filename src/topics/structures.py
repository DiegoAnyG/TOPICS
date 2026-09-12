"""Explicit heavy-atom selections; author chain and residue identifiers are retained."""

from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np


@dataclass
class Cloud:
    xyz: np.ndarray
    keys: list[tuple[str, int, str, str, str]]
    elements: list[str]

    def moved(self, rotation, translation):
        return Cloud(self.xyz @ rotation + translation, self.keys, self.elements)

    def subset(self, mask):
        indices = np.flatnonzero(mask)
        return Cloud(self.xyz[indices], [self.keys[i] for i in indices], [self.elements[i] for i in indices])


def read_cloud(path, chains=None, ligand=None, residue=None):
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise ValueError("Exactly one coordinate model is required; select a model before running.")
    selected = {}
    available = {chain.name for chain in structure[0]}
    if chains is not None and not set(chains) <= available:
        raise ValueError("Unknown author chain(s): " + ", ".join(sorted(set(chains) - available)))
    for chain in structure[0]:
        if chains is not None and chain.name not in chains:
            continue
        for res in chain:
            if ligand is not None:
                if res.name != ligand or (residue is not None and str(res.seqid) != str(residue)):
                    continue
            elif not gemmi.find_tabulated_residue(res.name).is_amino_acid():
                continue
            for atom in res:
                if atom.element.is_hydrogen or atom.occ <= 0:
                    continue
                key = (chain.name, res.seqid.num, res.seqid.icode.strip(), res.name, atom.name)
                priority = (-atom.occ, atom.altloc)
                if key not in selected or priority < selected[key][0]:
                    selected[key] = (priority, [atom.pos.x, atom.pos.y, atom.pos.z], atom.element.name)
    if not selected:
        raise ValueError(f"No atoms match the selection in {Path(path).name}.")
    keys = sorted(selected)
    xyz = np.array([selected[k][1] for k in keys], dtype=float)
    if not np.isfinite(xyz).all():
        raise ValueError("Coordinates must be finite.")
    if ligand and len({k[:4] for k in keys}) != 1:
        raise ValueError("Ligand selection is ambiguous; specify one chain and residue.")
    return Cloud(xyz, keys, [selected[k][2] for k in keys])


def write_complex(path, clouds):
    """Write mmCIF with role-prefixed chains to avoid identifier collisions."""
    structure = gemmi.Structure()
    structure.name = "TOPICS"
    model = gemmi.Model("1")
    for role, cloud in clouds.items():
        chains, residues = {}, {}
        for xyz, key, element in zip(cloud.xyz, cloud.keys, cloud.elements):
            chain_id, num, insertion, resname, atomname = key
            if not chain_id.startswith(role + "_"):
                chain_id = f"{role}_{chain_id}"
            if chain_id not in chains:
                chains[chain_id] = model.add_chain(gemmi.Chain(chain_id))
            rkey = (chain_id, num, insertion, resname)
            if rkey not in residues:
                res = gemmi.Residue()
                res.name = resname
                res.seqid = gemmi.SeqId(num, insertion or " ")
                res.het_flag = "H" if role == "PROTAC" else "A"
                residues[rkey] = chains[chain_id].add_residue(res)
            atom = gemmi.Atom()
            atom.name, atom.element = atomname, gemmi.Element(element)
            atom.pos = gemmi.Position(*xyz)
            atom.occ = 1.0
            residues[rkey].add_atom(atom)
    structure.add_model(model)
    structure.make_mmcif_document().write_file(str(path))
