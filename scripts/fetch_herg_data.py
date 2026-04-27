"""
Fetch real hERG IC50 data from ChEMBL.
Replaces synthetic training data in src/model.py.

Target: hERG (KCNH2) potassium channel
ChEMBL Target ID: CHEMBL240
Classification threshold: IC50 <= 1000 nM = cardiotoxic (label=1)

Run once to download, saves to data/herg_chembl.csv
Requires: pip install chembl_webresource_client pandas
"""

import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import logging
import os

logging.basicConfig(level=logging.INFO)

def fetch_herg_data(max_records: int = 2000) -> pd.DataFrame:
    """
    Fetch hERG IC50 data from ChEMBL.
    Returns DataFrame with SMILES, IC50, and binary label.
    """
    activity = new_client.activity

    logging.info("Fetching hERG IC50 data from ChEMBL...")
    # CHEMBL240 = hERG (human ether-a-go-go related gene) channel
    results = activity.filter(
        target_chembl_id="CHEMBL240",
        standard_type="IC50",
        standard_units="nM",
        standard_relation="="       # exact measurements only
    ).only([
        "molecule_chembl_id",
        "canonical_smiles",
        "standard_value",
        "standard_units"
    ])

    records = []
    for i, r in enumerate(results):
        if i >= max_records:
            break
        smiles = r.get("canonical_smiles")
        value  = r.get("standard_value")
        if smiles and value:
            try:
                records.append({
                    "smiles": smiles,
                    "ic50_nM": float(value),
                    "chembl_id": r.get("molecule_chembl_id")
                })
            except (ValueError, TypeError):
                continue

    df = pd.DataFrame(records)
    logging.info(f"Fetched {len(df)} records")

    # Label: IC50 <= 1000 nM = cardiotoxic blocker (label=1)
    # This threshold follows literature convention
    # Reference: Aronov (2008), Sanguinetti & Tristani-Firouzi (2006)
    df["label"] = (df["ic50_nM"] <= 1000).astype(int)
    logging.info(
        f"Class balance — Toxic (<=1µM): {df['label'].sum()} | "
        f"Safe (>1µM): {(df['label']==0).sum()}"
    )
    return df


def compute_features(df: pd.DataFrame) -> tuple:
    """
    Compute RDKit descriptors for each SMILES.
    Returns X (feature matrix) and y (labels).
    Drops rows where mol parsing fails.
    """
    from src.descriptors import compute_descriptors
    import dataclasses

    X_rows, y_rows = [], []
    failed = 0

    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None:
            failed += 1
            continue
        try:
            desc = compute_descriptors(mol)
            # Use descriptor features only (no tanimoto — 
            # those are pairwise, not per-molecule)
            features = [
                desc.mw, desc.logp, desc.tpsa,
                desc.hbd, desc.hba, desc.rotatable_bonds,
                desc.aromatic_rings, desc.fraction_csp3,
                desc.heavy_atom_count
            ]
            X_rows.append(features)
            y_rows.append(row["label"])
        except Exception:
            failed += 1

    logging.info(f"Feature extraction: {len(X_rows)} ok, {failed} failed")
    return np.array(X_rows, dtype=float), np.array(y_rows)


def save_dataset(df: pd.DataFrame, path: str = "data/herg_chembl.csv"):
    os.makedirs("data", exist_ok=True)
    df.to_csv(path, index=False)
    logging.info(f"Saved to {path}")


if __name__ == "__main__":
    df = fetch_herg_data(max_records=2000)
    save_dataset(df)
    X, y = compute_features(df)
    logging.info(f"Final dataset: {X.shape[0]} molecules, "
                 f"{y.mean()*100:.1f}% cardiotoxic")
    print(f"\nDataset ready: {X.shape[0]} samples")
    print(f"Positive (toxic): {y.sum()} | Negative (safe): {(y==0).sum()}")
