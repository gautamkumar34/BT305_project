"""
Molecular embedding module.
Handles 3D conformation generation using ETKDG v3 and MMFF94 optimization.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMoleculeConfs
import logging
from typing import Tuple, List, Optional

logger = logging.getLogger(__name__)

class MolecularEmbedder:
    """
    Generates high-quality 3D conformations for molecules.
    """
    
    def __init__(self, n_confs: int = 20, random_seed: int = 42):
        self.n_confs = n_confs
        self.random_seed = random_seed

    def embed_molecule(self, smiles: str) -> Tuple[Chem.Mol, int, List[float]]:
        """
        Generate multiple conformers using ETKDG v3, optimize with MMFF94,
        and return the molecule with the lowest energy conformer identified.
        
        Returns:
            mol: RDKit Mol with all conformers embedded
            best_conf_id: Index of the lowest-energy conformer
            energies: List of all conformer energies in kcal/mol
            
        Raises:
            ValueError: If no conformers could be generated.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Failed to parse SMILES: {smiles}")
            return None, -1, []

        # Add hydrogens for realistic 3D geometry
        mol = Chem.AddHs(mol)
        
        # ETKDG v3 parameters
        params = AllChem.ETKDGv3()
        params.randomSeed = self.random_seed
        
        # Generate multiple conformers
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=self.n_confs, params=params)
        
        if not conf_ids:
            raise ValueError(f"Failed to generate any conformers for SMILES: {smiles}")
        
        # Optimize geometry using MMFF94
        # MMFFOptimizeMoleculeConfs returns a list of (success, energy)
        try:
            results = MMFFOptimizeMoleculeConfs(mol)
        except Exception as e:
            logger.warning(f"MMFF optimization failed: {e}. Falling back to unoptimized conformers.")
            # If MMFF fails, we can't get energies easily, so we return the first one
            return mol, 0, [0.0] * len(conf_ids)

        energies = []
        for success, energy in results:
            # If optimization failed for a specific conformer, assign a high energy
            energies.append(energy if success else 9999.9)
            
        best_conf_id = energies.index(min(energies))
        
        logger.info(f"Generated {len(conf_ids)} conformers for {smiles}. "
                    f"Energy range: {min(energies):.2f} to {max(energies):.2f} kcal/mol. "
                    f"Best conf ID: {best_conf_id}")
        
        return mol, best_conf_id, energies

if __name__ == "__main__":
    # Quick test
    embedder = MolecularEmbedder()
    try:
        mol, best_id, energies = embedder.embed_molecule("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
        print(f"Successfully embedded. Best conformer ID: {best_id}")
    except Exception as e:
        print(f"Error: {e}")
