"""
FastAPI wrapper for the Drug Similarity and Cardiotoxicity Pipeline.
Provides an API for molecule comparison and validation.
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import Optional, Dict, Any, List
import uvicorn

from src.embedding import MolecularEmbedder
from src.similarity import MolecularSimilarity
from src.scoring import SimilarityScorer
from src.validation import ToxicityValidator

app = FastAPI(
    title="Drug Similarity & Cardiotoxicity API",
    description="Production-grade pipeline for ligand-based drug similarity analysis and cardiotoxicity validation",
    version="1.0.0"
)

# Configure CORS to allow requests from the frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, replace with specific frontend URL
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize components globally for reuse
embedder = MolecularEmbedder()
sim_calc = MolecularSimilarity()
scorer = SimilarityScorer()

class ComparisonRequest(BaseModel):
    smiles_a: str
    smiles_b: str
    alpha: Optional[float] = 0.5
    beta: Optional[float] = 0.5

class ComparisonResponse(BaseModel):
    molecule_a: str
    molecule_b: str
    tanimoto_2d: float
    shape_3d: float
    o3a_score: float
    rmsd: Optional[float] = None
    overlap_score: float
    final_score: float
    weights: Dict[str, float]
    explanation: str
    # Risk Assessment
    risk_flag_a: Optional[str] = None
    risk_flag_b: Optional[str] = None
    risk_rules_a: List[str] = Field(default_factory=list)
    risk_rules_b: List[str] = Field(default_factory=list)
    # 3D Alignment Data
    alignment: Optional[Dict[str, str]] = None
    # Physicochemical Descriptors
    descriptors_a: Optional[Dict[str, Any]] = None
    descriptors_b: Optional[Dict[str, Any]] = None
    descriptor_delta: Optional[Dict[str, Any]] = None
    ml_toxicity: Optional[Dict[str, Any]] = None

@app.get("/")
async def root():
    return {"message": "Drug Similarity API is running. Use /compare or /validate endpoints."}

@app.get("/health")
async def health_check():
    return {"status": "online", "version": "1.0.0"}

@app.post("/analyze", response_model=ComparisonResponse)
@app.post("/compare", response_model=ComparisonResponse)
async def compare_molecules(request: ComparisonRequest):
    """
    Compare two molecules using 2D and 3D similarity.
    """
    try:
        # 1. Generate 3D Embeddings
        res_a = embedder.embed_molecule(request.smiles_a)
        res_b = embedder.embed_molecule(request.smiles_b)
        
        if res_a is None or res_b is None:
            raise HTTPException(status_code=400, detail="Failed to generate 3D embeddings for one or both molecules")
        
        mol_a_3d, conf_id_a, _ = res_a
        mol_b_3d, conf_id_b, _ = res_b
        
        # 2. Compute Combined Similarity
        sim_results = sim_calc.compute_combined_similarity(
            request.smiles_a, request.smiles_b, mol_a_3d, mol_b_3d,
            conf_id_a=conf_id_a, conf_id_b=conf_id_b,
            alpha=request.alpha, beta=request.beta
        )
        
        # 3. Apply Scoring and Risk Flagging
        temp_scorer = SimilarityScorer(alpha=request.alpha, beta=request.beta)
        final_results = temp_scorer.compute_final_score(
            sim_results, 
            mol_a_smiles=request.smiles_a, 
            mol_b_smiles=request.smiles_b
        )
        
        return final_results
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/validate")
async def validate_pipeline():
    """
    Run the built-in validation cases (e.g., Terfenadine vs Fexofenadine).
    """
    try:
        validator = ToxicityValidator(embedder, sim_calc, scorer)
        results = validator.run_all_tests()
        return results
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/validate/{case_id}")
async def validate_case(case_id: str):
    """
    Run a single validation case by case_id.
    """
    try:
        validator = ToxicityValidator(embedder, sim_calc, scorer)
        case = validator.resolve_case(case_id)
        return validator.run_case(case)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)