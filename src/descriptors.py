from rdkit.Chem import Descriptors, rdMolDescriptors
from dataclasses import dataclass, asdict
from rdkit import Chem

@dataclass
class MolecularDescriptors:
    mw: float
    logp: float
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    aromatic_rings: int
    fraction_csp3: float
    heavy_atom_count: int

def compute_descriptors(mol: Chem.Mol) -> MolecularDescriptors:
    """Compute full descriptor set. Raises ValueError if mol is None."""
    if mol is None:
        raise ValueError("Cannot compute descriptors for None molecule")
    return MolecularDescriptors(
        mw=Descriptors.MolWt(mol),
        logp=Descriptors.MolLogP(mol),
        tpsa=Descriptors.TPSA(mol),
        hbd=rdMolDescriptors.CalcNumHBD(mol),
        hba=rdMolDescriptors.CalcNumHBA(mol),
        rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(mol),
        aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(mol),
        fraction_csp3=rdMolDescriptors.CalcFractionCSP3(mol),
        heavy_atom_count=mol.GetNumHeavyAtoms()
    )

def descriptor_delta(a: MolecularDescriptors, 
                     b: MolecularDescriptors) -> dict:
    """
    Returns signed difference (b - a) for each numeric field.
    Positive = b is larger. Used to explain structural changes.
    """
    a_dict = asdict(a)
    b_dict = asdict(b)
    return {
        k: round(b_dict[k] - a_dict[k], 4) 
        for k in a_dict
    }
