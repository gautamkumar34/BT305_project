"""
FastAPI wrapper for the Drug Similarity and Cardiotoxicity Pipeline.
Provides an API for molecule comparison, validation, and metabolism prediction.
"""

import os
# macOS Apple Silicon OpenBLAS / Python 3.13 Thread/Fork Deadlock Workarounds
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import Optional, Dict, Any, List
import uvicorn

# Internal Engine Imports
from src.embedding import MolecularEmbedder
from src.similarity import MolecularSimilarity
from src.scoring import SimilarityScorer
from src.validation import ToxicityValidator
from src.metabolism import MetabolismEngine

app = FastAPI(
    title="BioTransformer: Drug Similarity & Cardiotoxicity API",
    description="Production-grade pipeline for ligand-based drug similarity analysis, cardiotoxicity validation, and metabolic trajectory simulation.",
    version="1.1.0"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize core pipeline components globally so they stay warm in RAM
embedder = MolecularEmbedder()
sim_calc = MolecularSimilarity()
scorer = SimilarityScorer()
metabolism_engine = MetabolismEngine()

# --- Request & Response Models ---

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
    risk_flag_a: Optional[str] = None
    risk_flag_b: Optional[str] = None
    risk_rules_a: List[str] = Field(default_factory=list)
    risk_rules_b: List[str] = Field(default_factory=list)
    alignment: Optional[Dict[str, str]] = None
    descriptors_a: Optional[Dict[str, Any]] = None
    descriptors_b: Optional[Dict[str, Any]] = None
    descriptor_delta: Optional[Dict[str, Any]] = None
    ml_toxicity: Optional[Dict[str, Any]] = None

class MetabolismRequest(BaseModel):
    smiles: str
    max_depth: Optional[int] = 2
    mode: Optional[str] = "liver"

# --- API Endpoints ---

@app.get("/")
async def root():
    return {"message": "BioTransformer API is running securely."}

@app.get("/health")
async def health_check():
    return {"status": "online", "version": "1.1.0"}

@app.post("/analyze", response_model=ComparisonResponse)
@app.post("/compare", response_model=ComparisonResponse)
async def compare_molecules(request: ComparisonRequest):
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
        
        # 3. Apply Scoring and MPNN Toxicity checks
        temp_scorer = SimilarityScorer(alpha=request.alpha, beta=request.beta)
        return temp_scorer.compute_final_score(sim_results, request.smiles_a, request.smiles_b)
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/validate")
async def validate_pipeline():
    try:
        validator = ToxicityValidator(embedder, sim_calc, scorer)
        return validator.run_all_tests()
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/validate/{case_id}")
async def validate_case(case_id: str):
    try:
        validator = ToxicityValidator(embedder, sim_calc, scorer)
        case = validator.resolve_case(case_id)
        return validator.run_case(case)
    except Exception as e:
        raise HTTPException(status_code=404, detail=str(e))

@app.post("/metabolize")
async def metabolize_molecule(request: MetabolismRequest):
    try:
        # Pass the mode to the engine!
        trajectory_tree = metabolism_engine.predict_metabolites(
            request.smiles, 
            max_depth=request.max_depth, 
            mode=request.mode
        )
        return trajectory_tree
        
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=7860)