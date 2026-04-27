"""
Main entry point for the Drug Similarity and Cardiotoxicity Pipeline.
Provides CLI interface for molecule comparison and validation.
"""

import argparse
import json
import logging
import sys
from src.embedding import MolecularEmbedder
from src.similarity import MolecularSimilarity
from src.scoring import SimilarityScorer
from src.validation import ToxicityValidator

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("DrugSimilarityPipeline")

def run_comparison(smiles_a: str, smiles_b: str, alpha: float, beta: float):
    """
    Run a full similarity analysis between two molecules.
    """
    logger.info(f"Comparing {smiles_a} and {smiles_b}...")
    
    # Initialize components
    embedder = MolecularEmbedder()
    sim_calc = MolecularSimilarity()
    scorer = SimilarityScorer(alpha=alpha, beta=beta)
    
    # 1. Generate 3D Embeddings
    mol_a_3d = embedder.embed_molecule(smiles_a)
    mol_b_3d = embedder.embed_molecule(smiles_b)
    
    if mol_a_3d is None or mol_b_3d is None:
        return {"error": "Failed to generate 3D embeddings for one or both molecules"}
    
    # 2. Compute Combined Similarity
    sim_results = sim_calc.compute_combined_similarity(
        smiles_a, smiles_b, mol_a_3d, mol_b_3d
    )
    
    # 3. Apply Scoring
    final_results = scorer.compute_final_score(
        sim_results,
        mol_a_smiles=smiles_a,
        mol_b_smiles=smiles_b,
    )
    
    return final_results

def run_validation():
    """
    Run the built-in validation cases.
    """
    logger.info("Running validation pipeline...")
    embedder = MolecularEmbedder()
    sim_calc = MolecularSimilarity()
    scorer = SimilarityScorer()
    validator = ToxicityValidator(embedder, sim_calc, scorer)
    
    return validator.run_all_tests()

def main():
    parser = argparse.ArgumentParser(description="Ligand-based Drug Similarity & Cardiotoxicity Pipeline")
    
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Compare command
    compare_parser = subparsers.add_parser("compare", help="Compare two molecules")
    compare_parser.add_argument("--a", required=True, help="SMILES of molecule A")
    compare_parser.add_argument("--b", required=True, help="SMILES of molecule B")
    compare_parser.add_argument("--alpha", type=float, default=0.5, help="Weight for 2D similarity")
    compare_parser.add_argument("--beta", type=float, default=0.5, help="Weight for 3D similarity")
    
    # Validate command
    subparsers.add_parser("validate", help="Run validation test cases")

    # Serve command
    subparsers.add_parser("serve", help="Start the FastAPI server")
    
    args = parser.parse_args()
    
    if args.command == "compare":
        result = run_comparison(args.a, args.b, args.alpha, args.beta)
        print(json.dumps(result, indent=2))
    elif args.command == "validate":
        results = run_validation()
        print(json.dumps(results, indent=2))
    elif args.command == "serve":
        import uvicorn
        from src.api import app
        logger.info("Starting API server on http://0.0.0.0:8000")
        uvicorn.run(app, host="0.0.0.0", port=8000)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
