"""
Molecular scoring module.
Combines Chemprop MPNN toxicity predictions with rule-based risk flagging
and chemically specific explanations.
"""

from typing import Dict, Any, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors
import logging

logger = logging.getLogger(__name__)

# Lazy-loaded Chemprop predictor (loaded once on first use)
_predictor = None

def _get_tox_predictor():
    """Get or initialize the Chemprop toxicity predictor (singleton)."""
    global _predictor
    if _predictor is None:
        try:
            from src.tox_model import get_predictor
            _predictor = get_predictor()
        except Exception as e:
            logger.warning(f"Chemprop predictor unavailable: {e}")
            _predictor = None
    return _predictor

class SimilarityScorer:
    """
    Combines 2D and 3D similarity scores and applies rule-based risk heuristics.
    """
    
    def __init__(self, alpha: float = 0.5, beta: float = 0.5):
        self.alpha = alpha
        self.beta = beta
        self._normalize_weights()

    def _normalize_weights(self):
        total = self.alpha + self.beta
        if total > 0:
            self.alpha /= total
            self.beta /= total
        else:
            self.alpha = 0.5
            self.beta = 0.5

    def compute_descriptors(self, smiles: str) -> Dict[str, float]:
        """Compute key molecular descriptors for risk flagging."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        
        return {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "tpsa": Descriptors.TPSA(mol),
            "rotb": Descriptors.NumRotatableBonds(mol),
            "fsp3": Descriptors.FractionCSP3(mol)
        }

    def rule_based_risk_flag(self, descriptors: Any) -> Tuple[str, List[str], float]:
        """
        Determine risk flag and estimated probability based on physicochemical properties.
        Returns (risk_level, rules, probability).
        """
        if hasattr(descriptors, "logp"):  # Dataclass
            logp = descriptors.logp
            tpsa = descriptors.tpsa
            mw = descriptors.mw
            hbd = descriptors.hbd
        else:  # Dict
            logp = descriptors.get("logp", 0.0)
            tpsa = descriptors.get("tpsa", 0.0)
            mw = descriptors.get("mw", 0.0)
            hbd = descriptors.get("hbd", 0)

        # LOW_RISK: high polarity overrides lipophilicity
        if tpsa >= 75 and hbd >= 2:
            return "LOW_RISK", ["High TPSA (≥75 Å²) + HBD≥2 reduces membrane permeability and hERG risk"], 0.12

        # HIGH_RISK: lipophilic + low polarity = classic hERG profile  
        if logp > 4.5 and tpsa < 75:
            rules = [f"High logP ({logp:.2f}) with low TPSA ({tpsa:.1f} Å²)"]
            if mw > 450:
                rules.append(f"Large molecule (MW={mw:.0f}) amplifies hERG concern")
            return "HIGH_RISK", rules, 0.88

        # MODERATE_RISK: everything else
        return "MODERATE_RISK", ["Intermediate lipophilicity/polarity profile"], 0.55

    def compute_final_score(self, similarity_results: Dict[str, Any], 
                            mol_a_smiles: str, mol_b_smiles: str) -> Dict[str, Any]:
        """
        Compute final combined score, risk flags for both molecules, and detailed explanation.
        """
        score_2d = similarity_results.get("tanimoto_2d", 0.0)
        score_3d = similarity_results.get("shape_3d", 0.0)
        
        final_score = (self.alpha * score_2d) + (self.beta * score_3d)
        
        # Compute descriptors for both molecules
        desc_a = self.compute_descriptors(mol_a_smiles)
        desc_b = self.compute_descriptors(mol_b_smiles)
        
        # Risk flags for both molecules (rule-based tier labels)
        risk_level_a, rules_a, _ = self.rule_based_risk_flag(desc_a)
        risk_level_b, rules_b, _ = self.rule_based_risk_flag(desc_b)

        # ML toxicity predictions via Chemprop MPNN
        predictor = _get_tox_predictor()
        if predictor is not None:
            prob_a = predictor.predict_toxicity_prob(mol_a_smiles)
            prob_b = predictor.predict_toxicity_prob(mol_b_smiles)
            tox_disclaimer = (
                "Chemprop MPNN prediction trained on ChEMBL/Tox21 IC50 data. "
                "Not for clinical use."
            )
        else:
            # Fallback to rule-based probabilities
            _, _, prob_a = self.rule_based_risk_flag(desc_a)
            _, _, prob_b = self.rule_based_risk_flag(desc_b)
            tox_disclaimer = "Prediction based on physicochemical heuristics (model not loaded)."

        # Generate explanation
        explanation = self.generate_explanation(
            "Molecule A", "Molecule B", desc_a, desc_b, 
            score_2d, score_3d, 
            risk_level_a, risk_level_b, 
            {"tpsa": desc_b.get("tpsa", 0) - desc_a.get("tpsa", 0), 
             "logp": desc_b.get("logp", 0) - desc_a.get("logp", 0)}
        )

        from src.descriptors import MolecularDescriptors, descriptor_delta
        
        # Build strict descriptors for delta calculation
        struct_a = MolecularDescriptors(
            mw=desc_a['mw'], logp=desc_a['logp'], tpsa=desc_a['tpsa'],
            hbd=desc_a['hbd'], hba=desc_a['hba'], 
            rotatable_bonds=desc_a['rotb'], fraction_csp3=desc_a['fsp3'],
            aromatic_rings=0, heavy_atom_count=0 # Placeholders for delta
        )
        struct_b = MolecularDescriptors(
            mw=desc_b['mw'], logp=desc_b['logp'], tpsa=desc_b['tpsa'],
            hbd=desc_b['hbd'], hba=desc_b['hba'], 
            rotatable_bonds=desc_b['rotb'], fraction_csp3=desc_b['fsp3'],
            aromatic_rings=0, heavy_atom_count=0
        )
        delta = descriptor_delta(struct_a, struct_b)
            
        return {
            **similarity_results,
            "final_score": round(final_score, 3),
            "risk_flag_a": risk_level_a,
            "risk_rules_a": rules_a,
            "risk_flag_b": risk_level_b,
            "risk_rules_b": rules_b,
            "explanation": explanation,
            "descriptors_a": {
                "MW": desc_a['mw'],
                "logP": desc_a['logp'],
                "HBD": desc_a['hbd'],
                "HBA": desc_a['hba'],
                "TPSA": desc_a['tpsa']
            },
            "descriptors_b": {
                "MW": desc_b['mw'],
                "logP": desc_b['logp'],
                "HBD": desc_b['hbd'],
                "HBA": desc_b['hba'],
                "TPSA": desc_b['tpsa']
            },
            "descriptor_delta": delta,
            "ml_toxicity": {
                "prob_a": round(float(prob_a), 4),
                "prob_b": round(float(prob_b), 4),
                "correct_ranking": prob_a > prob_b,
                "disclaimer": tox_disclaimer
            }
        }

    def generate_explanation(self, name_a: str, name_b: str, desc_a: Any, desc_b: Any,
                             tanimoto_2d: float, shape_3d: float, 
                             risk_a: str, risk_b: str, delta: Dict[str, float]) -> str:
        """
        Generate a scaffold-aware explanation based on physicochemical changes.
        """
        if tanimoto_2d >= 0.75:
            scaffold_str = "share a highly similar core scaffold"
        elif tanimoto_2d >= 0.50:
            scaffold_str = "share partial structural similarity"
        else:
            scaffold_str = "have low structural similarity"
        
        abs_delta_tpsa = abs(delta.get('tpsa', 0))
        abs_delta_logp = abs(delta.get('logp', 0))
        
        tpsa_a = desc_a.get("tpsa", 0) if isinstance(desc_a, dict) else getattr(desc_a, "tpsa", 0)
        tpsa_b = desc_b.get("tpsa", 0) if isinstance(desc_b, dict) else getattr(desc_b, "tpsa", 0)
        logp_a = desc_a.get("logp", 0) if isinstance(desc_a, dict) else getattr(desc_a, "logp", 0)
        logp_b = desc_b.get("logp", 0) if isinstance(desc_b, dict) else getattr(desc_b, "logp", 0)

        if abs_delta_tpsa >= 20:
            structural_change = (
                f"a significant TPSA change of {delta.get('tpsa', 0):+.1f} Å² "
                f"(from {tpsa_a:.1f} to {tpsa_b:.1f} Å²)"
            )
        elif abs_delta_logp >= 1.0:
            structural_change = (
                f"a logP change of {delta.get('logp', 0):+.2f} "
                f"(from {logp_a:.2f} to {logp_b:.2f})"
            )
        else:
            structural_change = (
                f"modest physicochemical changes "
                f"(ΔTPSA={delta.get('tpsa', 0):+.1f} Å², ΔlogP={delta.get('logp', 0):+.2f})"
            )
        
        explanation_text = (
            f"{name_a} ({risk_a}, TPSA={tpsa_a:.1f} Å², logP={logp_a:.2f}) and "
            f"{name_b} ({risk_b}, TPSA={tpsa_b:.1f} Å², logP={logp_b:.2f}) "
            f"{scaffold_str} (2D Tanimoto={tanimoto_2d:.3f}). "
            f"The key structural difference is {structural_change}, "
            f"which {'reduces' if delta.get('tpsa', 0) > 0 else 'increases'} "
            f"predicted membrane permeability."
        )
        return explanation_text
