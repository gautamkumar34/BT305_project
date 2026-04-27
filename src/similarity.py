"""
Molecular similarity computation module.
Computes both 2D (fingerprint-based) and 3D (shape-based) similarity.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdShapeHelpers
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
import numpy as np
from typing import Tuple, Optional, Dict, Any, List
import logging

logger = logging.getLogger(__name__)

class MolecularSimilarity:
    """
    Computes molecular similarity using multiple methods.
    
    Features:
    - 2D Tanimoto similarity using Morgan fingerprints
    - 3D shape similarity (Shape Tanimoto)
    - 3D alignment and RMSD calculation using O3A
    """
    
    def __init__(self, 
                 fingerprint_radius: int = 2,
                 fingerprint_nbits: int = 2048):
        """
        Initialize similarity calculator.
        """
        self.fingerprint_radius = fingerprint_radius
        self.fingerprint_nbits = fingerprint_nbits
        self.morgan_gen = GetMorganGenerator(radius=fingerprint_radius)
        
    def compute_2d_similarity(self, 
                             smiles_a: str, 
                             smiles_b: str,
                             mol_a: Optional[Chem.Mol] = None,
                             mol_b: Optional[Chem.Mol] = None) -> float:
        """
        Compute 2D Tanimoto similarity using Morgan fingerprints.
        """
        try:
            if mol_a is None:
                mol_a = Chem.MolFromSmiles(smiles_a)
            if mol_b is None:
                mol_b = Chem.MolFromSmiles(smiles_b)
            
            if mol_a is None or mol_b is None:
                logger.error("Failed to parse SMILES for 2D similarity")
                return 0.0
            
            fp_a = self.morgan_gen.GetFingerprint(mol_a)
            fp_b = self.morgan_gen.GetFingerprint(mol_b)
            
            similarity = AllChem.DataStructs.TanimotoSimilarity(fp_a, fp_b)
            return similarity
            
        except Exception as e:
            logger.error(f"Error computing 2D similarity: {e}")
            return 0.0
    
    def align_and_score_3d(self, mol_ref, mol_probe, conf_id_ref=0, conf_id_probe=0):
        """
        Align mol_probe onto mol_ref using O3A.
        Returns a dict with:
          - o3a_score: RDKit O3A alignment score (higher is better)
          - rmsd: RMSD after alignment (lower is better)
          - shape_tanimoto: shape similarity after alignment (higher is better)
        """
        try:
            # Step 1: compute MMFF properties (REQUIRED before GetO3A)
            mmff_ref   = AllChem.MMFFGetMoleculeProperties(mol_ref,   mmffVariant="MMFF94s")
            mmff_probe = AllChem.MMFFGetMoleculeProperties(mol_probe, mmffVariant="MMFF94s")

            if mmff_ref is None or mmff_probe is None:
                raise ValueError("MMFF94s properties returned None")

            # Step 2: O3A alignment
            o3a = rdMolAlign.GetO3A(
                mol_probe, mol_ref,
                mmff_probe, mmff_ref,
                prbCid=conf_id_probe,
                refCid=conf_id_ref
            )
            o3a_score = o3a.Align()   # .Align() does the alignment in-place, returns score

            # RMSD after alignment
            try:
                rmsd = rdMolAlign.CalcRMS(
                    mol_ref, mol_probe,
                    confId1=conf_id_ref,
                    confId2=conf_id_probe,
                )
            except Exception:
                rmsd = 99.0
            
            # Step 3: shape AFTER alignment (molecules now superimposed)
            shape_dist = rdShapeHelpers.ShapeTanimotoDist(
                mol_ref, mol_probe,
                confId1=conf_id_ref,
                confId2=conf_id_probe
            )
            shape_tanimoto = 1.0 - shape_dist

            logger.info(
                f"O3A success — Score: {o3a_score:.3f}, RMSD: {rmsd:.3f}, Shape: {shape_tanimoto:.3f}"
            )
            return {
                "o3a_score": round(float(o3a_score), 4),
                "rmsd": round(float(rmsd), 4),
                "shape_tanimoto": round(float(shape_tanimoto), 4),
            }

        except Exception as e:
            # Do not suppress exception, print full traceback as requested
            import traceback
            logger.error(f"O3A failed:\n{traceback.format_exc()}")
            # Fallback: shape only, no O3A score
            try:
                shape_dist = rdShapeHelpers.ShapeTanimotoDist(mol_ref, mol_probe)
                return {
                    "o3a_score": 99.0,
                    "rmsd": 99.0,
                    "shape_tanimoto": round(1.0 - shape_dist, 4),
                }
            except:
                return {"o3a_score": 99.0, "rmsd": 99.0, "shape_tanimoto": 0.0}

    def compute_3d_similarity(self,
                             mol_a: Chem.Mol,
                             mol_b: Chem.Mol,
                             conf_id_a: int = 0,
                             conf_id_b: int = 0) -> Dict[str, float]:
        """
        Compute 3D shape similarity and RMSD using the aligned scoring method.
        """
        try:
            if mol_a.GetNumConformers() == 0 or mol_b.GetNumConformers() == 0:
                logger.error("Molecules must have 3D coordinates for shape similarity")
                return {"shape_tanimoto": 0.0, "o3a_score": 99.0, "rmsd": 99.0, "overlap_score": 0.0}
            
            # Use the new align_and_score_3d method
            align = self.align_and_score_3d(mol_a, mol_b, conf_id_a, conf_id_b)
            
            # Protrude distance for overlap score
            try:
                protrude_dist = rdShapeHelpers.ShapeProtrudeDist(mol_a, mol_b)
                overlap_score = 1.0 / (1.0 + protrude_dist)
            except Exception as e:
                logger.warning(f"ProtrudeDist failed: {e}")
                overlap_score = 0.0
            
            return {
                "shape_tanimoto": align["shape_tanimoto"],
                "o3a_score": align["o3a_score"],
                "rmsd": align["rmsd"],
                "overlap_score": overlap_score
            }
            
        except Exception as e:
            logger.error(f"Error computing 3D similarity: {e}")
            return {"shape_tanimoto": 0.0, "o3a_score": 99.0, "rmsd": 99.0, "overlap_score": 0.0}

    def compute_combined_similarity(self,
                                    smiles_a: str,
                                    smiles_b: str,
                                    mol_a_3d: Optional[Chem.Mol] = None,
                                    mol_b_3d: Optional[Chem.Mol] = None,
                                    conf_id_a: int = 0,
                                    conf_id_b: int = 0,
                                    alpha: float = 0.5,
                                    beta: float = 0.5) -> Dict[str, Any]:
        """
        Compute combined 2D and 3D similarity and return alignment SDFs.
        """
        similarity_2d = self.compute_2d_similarity(smiles_a, smiles_b)
        
        similarity_3d = {"shape_tanimoto": 0.0, "o3a_score": 99.0, "rmsd": 99.0, "overlap_score": 0.0}
        alignment_sdfs = {"mol_a_sdf": "", "mol_b_sdf": ""}
        
        if mol_a_3d is not None and mol_b_3d is not None:
            similarity_3d = self.compute_3d_similarity(mol_a_3d, mol_b_3d, conf_id_a, conf_id_b)
            
            # Generate SDF strings for the aligned molecules
            try:
                # mol_b_3d was aligned to mol_a_3d in-place by align_and_score_3d
                alignment_sdfs["mol_a_sdf"] = Chem.MolToMolBlock(mol_a_3d)
                alignment_sdfs["mol_b_sdf"] = Chem.MolToMolBlock(mol_b_3d)
            except Exception as e:
                logger.error(f"Error generating alignment SDFs: {e}")

        total_weight = alpha + beta
        alpha_norm = alpha / total_weight if total_weight > 0 else 0.5
        beta_norm = beta / total_weight if total_weight > 0 else 0.5
        
        combined_score = (alpha_norm * similarity_2d + 
                         beta_norm * similarity_3d["shape_tanimoto"])
        
        return {
            "molecule_a": smiles_a,
            "molecule_b": smiles_b,
            "tanimoto_2d": similarity_2d,
            "shape_3d": similarity_3d["shape_tanimoto"],
            "o3a_score": similarity_3d["o3a_score"],
            "rmsd": similarity_3d.get("rmsd", 99.0),
            "overlap_score": similarity_3d["overlap_score"],
            "final_score": combined_score,
            "weights": {"alpha": alpha_norm, "beta": beta_norm},
            "alignment": alignment_sdfs
        }
