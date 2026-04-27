"""
Toxicity Prediction Model
-------------------------
IMPORTANT DISCLAIMER: This model uses SYNTHETIC training data 
constructed from known physicochemical heuristics for hERG blockade.
It is NOT trained on real experimental IC50 data.

For production use, replace X_train/y_train with data from:
ChEMBL assay CHEMBL1614027 (hERG IC50 measurements)

Reference for hERG pharmacophore rules used to generate 
synthetic data: Aronov AM (2008) Drug Discov Today 13:149-155
"""

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
import joblib
import logging
from src.descriptors import MolecularDescriptors

FEATURE_NAMES = [
    "tanimoto_2d", "shape_3d", "mw", "logp", "tpsa",
    "hbd", "hba", "rotatable_bonds", "aromatic_rings", "fraction_csp3"
]

def build_feature_vector(tanimoto_2d: float, shape_3d: float,
                         desc: MolecularDescriptors) -> np.ndarray:
    return np.array([
        tanimoto_2d, shape_3d,
        desc.mw, desc.logp, desc.tpsa,
        desc.hbd, desc.hba, desc.rotatable_bonds,
        desc.aromatic_rings, desc.fraction_csp3
    ], dtype=float)

def generate_synthetic_training_data() -> tuple[np.ndarray, np.ndarray]:
    """
    Synthetic data based on hERG pharmacophore rules.
    
    CARDIOTOXIC proxy (label=1): 
      High logP (4-7), low TPSA (20-60), MW 350-550, HBD 0-2
      Characteristic of flat, lipophilic hERG channel blockers
    
    SAFE proxy (label=0):
      High TPSA (70-120), logP 1-4, HBD 2-4
      Characteristic of polar, hydrophilic molecules
    
    50 samples each for balance. seed=42 for reproducibility.
    """
    rng = np.random.default_rng(42)
    
    # Cardiotoxic-like: [tanimoto_2d, shape_3d, mw, logp, tpsa,
    #                    hbd, hba, rot_bonds, arom_rings, fsp3]
    toxic = np.column_stack([
        rng.uniform(0.6, 0.95, 50),   # tanimoto_2d
        rng.uniform(0.4, 0.8,  50),   # shape_3d
        rng.uniform(350, 550,  50),   # mw
        rng.uniform(3.5, 7.0,  50),   # logp  ← HIGH
        rng.uniform(15,  80,   50),   # tpsa  ← LOW
        rng.integers(0, 2,     50),   # hbd
        rng.integers(2, 5,     50),   # hba
        rng.integers(3, 8,     50),   # rot_bonds
        rng.integers(2, 4,     50),   # arom_rings
        rng.uniform(0.1, 0.4,  50),   # fsp3
    ])
    
    safe = np.column_stack([
        rng.uniform(0.3, 0.75, 50),
        rng.uniform(0.2, 0.6,  50),
        rng.uniform(250, 500,  50),
        rng.uniform(0.5, 5.0,  50),   # logp  ← LOW
        rng.uniform(50,  130,  50),   # tpsa  ← HIGH
        rng.integers(2, 5,     50),   # hbd   ← HIGH
        rng.integers(3, 7,     50),
        rng.integers(2, 7,     50),
        rng.integers(1, 3,     50),
        rng.uniform(0.3, 0.7,  50),
    ])
    
    X = np.vstack([toxic, safe])
    y = np.array([1]*50 + [0]*50)
    return X, y

def get_training_data(chembl_csv: str = "data/herg_chembl.csv"):
    """
    Load real ChEMBL data if available, else fall back to synthetic.
    """
    import os
    import logging
    if os.path.exists(chembl_csv):
        import pandas as pd
        from scripts.fetch_herg_data import compute_features
        df = pd.read_csv(chembl_csv)
        X, y = compute_features(df)
        logging.info(f"Using REAL ChEMBL data: {len(y)} samples")
        data_source = "ChEMBL CHEMBL240 hERG IC50"
    else:
        X, y = generate_synthetic_training_data()
        logging.warning(
            "Using SYNTHETIC data. Run: "
            "python scripts/fetch_herg_data.py"
        )
        data_source = "SYNTHETIC (rule-based proxy)"

    return X, y, data_source

def train_toxicity_model() -> dict:
    """Train RF classifier. Returns model dict with metadata."""
    X, y, data_source = get_training_data()
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    clf = RandomForestClassifier(
        n_estimators=100, random_state=42, 
        max_depth=5, min_samples_leaf=3
    )
    
    cv_scores = cross_val_score(clf, X_scaled, y, cv=5, scoring="roc_auc")
    clf.fit(X_scaled, y)
    
    importance = dict(zip(FEATURE_NAMES, 
                         clf.feature_importances_.tolist()))
    
    logging.info(f"Model trained using {data_source}. CV ROC-AUC: {cv_scores.mean():.3f} "
                 f"± {cv_scores.std():.3f}")
    
    return {
        "model": clf,
        "scaler": scaler,
        "cv_roc_auc": round(cv_scores.mean(), 3),
        "feature_importances": importance,
        "disclaimer": (
            f"Trained on {data_source}. "
            "Not for clinical use. For production, use validated assay data."
        )
    }

def predict_toxicity(model_dict: dict, 
                     feature_vector: np.ndarray) -> dict:
    """
    Predict cardiotoxicity probability for a single molecule.
    Returns probability, label, and disclaimer.
    """
    clf    = model_dict["model"]
    scaler = model_dict["scaler"]
    
    x_scaled = scaler.transform(feature_vector.reshape(1, -1))
    prob = clf.predict_proba(x_scaled)[0][1]  # P(toxic)
    
    return {
        "probability": round(float(prob), 4),
        "label": "HIGH_RISK" if prob >= 0.5 else "LOW_RISK",
        "disclaimer": model_dict["disclaimer"]
    }

def save_model(model_dict: dict, path: str = "model.pkl"):
    """Save model and scaler. Exclude non-serializable keys."""
    joblib.dump({
        "model": model_dict["model"],
        "scaler": model_dict["scaler"],
        "cv_roc_auc": model_dict["cv_roc_auc"],
        "feature_importances": model_dict["feature_importances"],
        "disclaimer": model_dict["disclaimer"]
    }, path)
    logging.info(f"Model saved to {path}")

def load_model(path: str = "model.pkl") -> dict:
    return joblib.load(path)
