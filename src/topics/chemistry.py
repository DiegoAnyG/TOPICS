"""CCD chemistry and unconstrained ETKDG ensembles; no native ternary coordinates."""

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


def read_ccd(path):
    block = gemmi.cif.read_file(str(path)).sole_block()
    atoms = block.get_mmcif_category("_chem_comp_atom.")
    bonds = block.get_mmcif_category("_chem_comp_bond.")
    molecule, indices = Chem.RWMol(), {}
    for name, element, charge in zip(atoms["atom_id"], atoms["type_symbol"], atoms["charge"]):
        atom = Chem.Atom(element.title())
        atom.SetFormalCharge(int(charge or 0))
        atom.SetProp("atom_name", name)
        indices[name] = molecule.AddAtom(atom)
    types = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for a, b, order in zip(bonds["atom_id_1"], bonds["atom_id_2"], bonds["value_order"]):
        if order not in types:
            raise ValueError(f"Unsupported CCD bond order: {order}")
        molecule.AddBond(indices[a], indices[b], types[order])
    mol = molecule.GetMol()
    Chem.SanitizeMol(mol)
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        xyz = [float(atoms[f"pdbx_model_Cartn_{axis}_ideal"][i]) for axis in "xyz"]
        conformer.SetAtomPosition(i, xyz)
    mol.AddConformer(conformer)
    Chem.AssignStereochemistryFrom3D(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for i, expected in enumerate(atoms["pdbx_stereo_config"]):
        if expected in {"R", "S"}:
            atom = mol.GetAtomWithIdx(i)
            if not atom.HasProp("_CIPCode") or atom.GetProp("_CIPCode") != expected:
                raise ValueError(f"CCD stereochemistry disagrees with ideal geometry at {atoms['atom_id'][i]}.")
    mol = Chem.RemoveHs(mol)
    mol.SetProp("ccd_id", block.find_value("_chem_comp.id"))
    mol.RemoveAllConformers()
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("The PROTAC must be one connected molecule.")
    return mol, [a.GetProp("atom_name") for a in mol.GetAtoms()]


def descriptors(mol):
    return {"molecular_weight": Descriptors.MolWt(mol), "crippen_logp": Descriptors.MolLogP(mol),
            "tpsa_A2": rdMolDescriptors.CalcTPSA(mol), "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True)}


def generate(mol, count, seed, threads=1):
    molecule = Chem.AddHs(mol)
    molecule.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed, params.numThreads = seed, threads
    params.enforceChirality, params.useRandomCoords = True, True
    params.maxIterations = 500
    ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=count, params=params))
    if not ids:
        raise ValueError("No conformers generated. Check chemistry or increase the sampling budget.")
    if not AllChem.MMFFHasAllMoleculeParams(molecule):
        raise ValueError("MMFF parameters are missing for this molecule; no energy fallback is applied.")
    results = AllChem.MMFFOptimizeMoleculeConfs(molecule, numThreads=threads, maxIters=1500, mmffVariant="MMFF94s")
    expected = {a.GetIdx(): a.GetProp("_CIPCode") for a in mol.GetAtoms() if a.HasProp("_CIPCode")}
    ensemble = []
    for cid, (status, energy) in zip(ids, results):
        check = Chem.Mol(molecule)
        Chem.AssignStereochemistryFrom3D(check, confId=cid, replaceExistingTags=True)
        Chem.AssignStereochemistry(check, cleanIt=True, force=True)
        if any(not check.GetAtomWithIdx(i).HasProp("_CIPCode") or check.GetAtomWithIdx(i).GetProp("_CIPCode") != cip for i, cip in expected.items()):
            continue
        xyz = np.array(molecule.GetConformer(cid).GetPositions())[:mol.GetNumAtoms()]
        if status == 0 and np.isfinite(xyz).all() and np.isfinite(energy):
            ensemble.append((cid, xyz, float(energy)))
    if not ensemble:
        raise ValueError("No converged, stereochemistry-preserving conformers; increase the budget or inspect chemistry.")
    return ensemble, {"requested": count, "embedded": len(ids), "accepted": len(ensemble)}
