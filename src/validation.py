"""
Validation module for the Drug Similarity and Cardiotoxicity Pipeline.
Contains test cases and validation logic to ensure scientific correctness.
"""

import logging
from typing import List, Dict, Any, Optional

from src.constants import VERIFIED_SMILES
from src.embedding import MolecularEmbedder
from src.similarity import MolecularSimilarity
from src.scoring import SimilarityScorer

logger = logging.getLogger(__name__)

class ToxicityValidator:
    """
    Validates the pipeline using known clinical cases of cardiotoxicity.
    """
    
    def __init__(self, embedder: MolecularEmbedder, sim_calc: MolecularSimilarity, scorer: SimilarityScorer):
        self.embedder = embedder
        self.sim_calc = sim_calc
        self.scorer = scorer

    def get_cases(self) -> List[Dict[str, Any]]:
        """
        Canonical validation cases based on VERIFIED_SMILES.
        """
        return [
            {
                "id": "terfenadine_fexofenadine",
                "label": "Terfenadine vs Fexofenadine",
                "a_key": "terfenadine",
                "b_key": "fexofenadine",
                "expected": "High scaffold similarity; Fexofenadine should be lower risk",
            },
            {
                "id": "cisapride_domperidone",
                "label": "Cisapride vs Domperidone",
                "a_key": "cisapride",
                "b_key": "domperidone",
                "expected": "Low-to-moderate similarity; both can present risk in practice",
            },
            {
                "id": "aspirin_ibuprofen",
                "label": "Aspirin vs Ibuprofen",
                "a_key": "aspirin",
                "b_key": "ibuprofen",
                "expected": "Low similarity; both generally not hERG-classic blockers",
            },
        ]

    def resolve_case(self, case_id: str) -> Dict[str, Any]:
        """
        Resolve a case by id (stable API contract).
        """
        for c in self.get_cases():
            if c["id"] == case_id:
                return c
        raise ValueError(f"Unknown case_id: {case_id}")

    def run_case(self, case: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run a single validation case.
        """
        a = VERIFIED_SMILES[case["a_key"]]
        b = VERIFIED_SMILES[case["b_key"]]

        smiles_a = a["smiles"]
        smiles_b = b["smiles"]

        logger.info(f"Running validation case: {case['label']}")

        # 1. Embed
        res_a = self.embedder.embed_molecule(smiles_a)
        res_b = self.embedder.embed_molecule(smiles_b)

        mol_a, conf_a, _ = res_a
        mol_b, conf_b, _ = res_b

        # 2. Similarity
        sim = self.sim_calc.compute_combined_similarity(
            smiles_a, smiles_b, mol_a, mol_b,
            conf_id_a=conf_a, conf_id_b=conf_b
        )

        # 3. Scoring & Risk
        final = self.scorer.compute_final_score(sim, smiles_a, smiles_b)

        return {
            "case_id": case["id"],
            "case": case["label"],
            "expected": case.get("expected"),
            "molecule_a": {"key": case["a_key"], **a},
            "molecule_b": {"key": case["b_key"], **b},
            "result": final,
        }

    def run_all_tests(self) -> List[Dict[str, Any]]:
        """
        Runs a suite of validation tests.
        """
        results = []
        for case in self.get_cases():
            try:
                out = self.run_case(case)
                final = out["result"]
                out["status"] = "PASS" if final.get("final_score", 1.0) < 0.8 else "REVIEW"
                results.append(out)
            except Exception as e:
                logger.error(f"Validation case {case.get('id')} failed: {e}")
                results.append({"case_id": case.get("id"), "case": case.get("label"), "status": "FAIL", "error": str(e)})
                
        return results

if __name__ == "__main__":
    # Setup for standalone test
    import logging as py_logging
    py_logging.basicConfig(level=py_logging.INFO)
    
    embedder = MolecularEmbedder()
    sim_calc = MolecularSimilarity()
    scorer = SimilarityScorer()
    
    validator = ToxicityValidator(embedder, sim_calc, scorer)
    print(validator.run_all_tests())