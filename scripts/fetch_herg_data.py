"""
Fetch and assemble hERG cardiotoxicity training dataset.

Data sources:
  1. ChEMBL REST API — three hERG-related targets (direct HTTP, limit=1000/page)
  2. Tox21 — NIH screening dataset (via requests to bypass SSL issues)
  3. Curated ground truth — pharmacologically verified molecules

Labeling:
  IC50 <= 1000 nM  → label=1 (toxic / hERG blocker)
  IC50 >= 10000 nM → label=0 (safe / non-blocker)
  1000 < IC50 < 10000 → EXCLUDED (ambiguous zone)

Output: data/herg_training_final.csv  (smiles, label, source)
"""

import logging
import os
import gzip
import io
from concurrent.futures import ThreadPoolExecutor, as_completed

# macOS Apple Silicon OpenBLAS / Python 3.13 Thread/Fork Deadlock Workarounds
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"

import pandas as pd
import numpy as np
import requests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────
# SMILES utilities
# ─────────────────────────────────────────────────────────

def canonicalize_smiles(smi: str):
    """Return canonical SMILES or None if invalid."""
    from rdkit import Chem
    if not smi or not isinstance(smi, str):
        return None
    mol = Chem.MolFromSmiles(smi.strip())
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def deduplicate_df(df: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate by canonical SMILES, keeping first occurrence."""
    before = len(df)
    df = df.drop_duplicates(subset=["smiles"], keep="first").reset_index(drop=True)
    logger.info(f"Deduplication: {before} → {len(df)} rows")
    return df


# ─────────────────────────────────────────────────────────
# SOURCE 1 — ChEMBL REST API (direct HTTP, fast pagination)
# ─────────────────────────────────────────────────────────

CHEMBL_TARGETS = {
    "CHEMBL240":  "hERG primary assay",
    "CHEMBL4282": "hERG patch clamp",
    "CHEMBL3307": "IKr potassium channel",
}

CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"


def fetch_chembl_target(target_id: str, max_records: int = 8000) -> pd.DataFrame:
    """Fetch IC50 data via ChEMBL REST API with large page sizes."""
    logger.info(f"Fetching {target_id} ({CHEMBL_TARGETS.get(target_id, '')})...")

    records = []
    offset = 0
    page_size = 1000  # Much faster than default 20

    while len(records) < max_records:
        params = {
            "target_chembl_id": target_id,
            "standard_type": "IC50",
            "standard_units": "nM",
            "limit": page_size,
            "offset": offset,
            "format": "json",
        }

        try:
            resp = requests.get(CHEMBL_API_BASE, params=params, timeout=60)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            logger.warning(f"  API error at offset {offset}: {e}")
            break

        activities = data.get("activities", [])
        if not activities:
            break

        for r in activities:
            if len(records) >= max_records:
                break

            smi = r.get("canonical_smiles")
            val = r.get("standard_value")
            rel = r.get("standard_relation", "=")

            if not smi or val is None:
                continue
            if rel not in ("=", "<", ">", "<=", ">="):
                continue

            try:
                ic50 = float(val)
            except (ValueError, TypeError):
                continue

            # Ambiguous zone exclusion
            if 1000 < ic50 < 10000:
                continue

            can_smi = canonicalize_smiles(smi)
            if can_smi is None:
                continue

            if ic50 <= 1000:
                label = 1  # toxic
            else:
                label = 0  # safe (>= 10000)

            records.append({
                "smiles": can_smi,
                "label": label,
                "source": f"ChEMBL_{target_id}",
            })

        logger.info(f"  {target_id}: page offset={offset}, total records so far: {len(records)}")

        # Check if there are more pages
        next_url = data.get("page_meta", {}).get("next")
        if not next_url:
            break
        offset += page_size

    df = pd.DataFrame(records)
    if len(df) > 0:
        toxic_count = df['label'].sum()
        safe_count = (df['label'] == 0).sum()
        logger.info(f"  {target_id}: {len(df)} total (toxic={toxic_count}, safe={safe_count})")
    else:
        logger.warning(f"  {target_id}: 0 records fetched")
    return df


def fetch_chembl_multi_assay(max_per_target: int = 8000) -> pd.DataFrame:
    """Fetch from all three ChEMBL targets in parallel."""
    dfs = []

    # Parallel fetch for speed
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {
            executor.submit(fetch_chembl_target, tid, max_per_target): tid
            for tid in CHEMBL_TARGETS
        }
        for future in as_completed(futures):
            tid = futures[future]
            try:
                df = future.result()
                dfs.append(df)
            except Exception as e:
                logger.error(f"Failed to fetch {tid}: {e}")

    if not dfs:
        logger.error("No ChEMBL data fetched!")
        return pd.DataFrame(columns=["smiles", "label", "source"])

    combined = pd.concat(dfs, ignore_index=True)
    combined = deduplicate_df(combined)
    logger.info(
        f"ChEMBL combined: {len(combined)} unique molecules "
        f"(toxic={combined['label'].sum()}, safe={(combined['label']==0).sum()})"
    )
    return combined


# ─────────────────────────────────────────────────────────
# SOURCE 2 — Tox21 hERG Dataset
# ─────────────────────────────────────────────────────────

def fetch_tox21_herg() -> pd.DataFrame:
    """Download Tox21 dataset using requests (handles SSL properly)."""
    urls = [
        "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz",
        "https://raw.githubusercontent.com/deepchem/deepchem/master/datasets/tox21.csv.gz",
    ]

    tox21_df = None
    for url in urls:
        try:
            logger.info(f"Downloading Tox21 from {url}...")
            resp = requests.get(url, timeout=60, verify=True)
            resp.raise_for_status()
            data = resp.content

            if url.endswith(".gz"):
                data = gzip.decompress(data)
            tox21_df = pd.read_csv(io.BytesIO(data))
            logger.info(f"Tox21 downloaded: {len(tox21_df)} rows, columns: {list(tox21_df.columns)}")
            break
        except requests.exceptions.SSLError:
            # Retry without SSL verification as fallback
            try:
                logger.warning(f"SSL error, retrying without verification...")
                resp = requests.get(url, timeout=60, verify=False)
                resp.raise_for_status()
                data = resp.content
                if url.endswith(".gz"):
                    data = gzip.decompress(data)
                tox21_df = pd.read_csv(io.BytesIO(data))
                logger.info(f"Tox21 downloaded (no SSL verify): {len(tox21_df)} rows")
                break
            except Exception as e2:
                logger.warning(f"Tox21 download failed from {url}: {e2}")
                continue
        except Exception as e:
            logger.warning(f"Tox21 download failed from {url}: {e}")
            continue

    if tox21_df is None:
        logger.error("All Tox21 download URLs failed")
        return pd.DataFrame(columns=["smiles", "label", "source"])

    # Find SMILES and hERG columns
    smiles_col = next((c for c in tox21_df.columns if "smiles" in c.lower()), None)
    herg_col = next((c for c in tox21_df.columns if "herg" in c.lower()), None)

    if smiles_col is None:
        logger.error(f"No SMILES column found. Columns: {list(tox21_df.columns)}")
        return pd.DataFrame(columns=["smiles", "label", "source"])

    if herg_col is None:
        logger.warning(f"No hERG column found. Columns: {list(tox21_df.columns)}")
        return pd.DataFrame(columns=["smiles", "label", "source"])

    logger.info(f"Using Tox21 columns: smiles={smiles_col}, herg={herg_col}")

    # Filter valid rows
    valid = tox21_df[[smiles_col, herg_col]].dropna(subset=[herg_col])
    valid = valid[valid[herg_col].isin([0.0, 1.0])]

    records = []
    for _, row in valid.iterrows():
        can_smi = canonicalize_smiles(row[smiles_col])
        if can_smi is None:
            continue
        records.append({
            "smiles": can_smi,
            "label": int(row[herg_col]),
            "source": "Tox21",
        })

    df = pd.DataFrame(records)
    df = deduplicate_df(df)
    logger.info(
        f"Tox21 hERG: {len(df)} molecules "
        f"(toxic={df['label'].sum()}, safe={(df['label']==0).sum()})"
    )
    return df


# ─────────────────────────────────────────────────────────
# SOURCE 3 — Curated Ground Truth Correction Set
# ─────────────────────────────────────────────────────────

CURATED_CORRECTIONS = [
    # KNOWN SAFE (label=0)
    {"smiles": "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc4cccnc24)CC1",
     "label": 0, "name": "Loratadine"},
    {"smiles": "CC(C)(C(=O)O)c1ccc(C(O)CCCN2CCC(CC2)C(O)(c2ccccc2)c2ccccc2)cc1",
     "label": 0, "name": "Fexofenadine"},
    {"smiles": "OC(=O)CN1CCN(CC1)c1ccc(Cl)cc1C(c1ccccc1)n1ccnc1",
     "label": 0, "name": "Cetirizine"},
    {"smiles": "CC(=O)Oc1ccccc1C(=O)O",
     "label": 0, "name": "Aspirin"},
    {"smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
     "label": 0, "name": "Ibuprofen"},
    {"smiles": "CC(C)NCC(O)COc1ccc(CCOC)cc1",
     "label": 0, "name": "Metoprolol"},
    {"smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C=C[C@@H]3[C@H]1C5",
     "label": 0, "name": "Morphine"},
    {"smiles": "Clc1ccc2c(c1)CC1=CCCCN1c2=O",
     "label": 0, "name": "Desloratadine"},
    # KNOWN TOXIC (label=1)
    {"smiles": "CC(C)(C)c1ccc(C(O)CCCN2CCC(CC2)C(O)(c2ccccc2)c2ccccc2)cc1",
     "label": 1, "name": "Terfenadine"},
    {"smiles": "COc1ccc(CCN(C)C2CCN(Cc3nc4ccccc4[nH]3)CC2)cc1",
     "label": 1, "name": "Astemizole"},
    {"smiles": "OC(CCN1CCC(CC1)Nc1ncc(Cl)cc1OC)c1ccc(F)cc1",
     "label": 1, "name": "Cisapride"},
    {"smiles": "OC(CCCN1CCC(CC1)c2ccc(Cl)cc2)c1ccc(F)cc1",
     "label": 1, "name": "Haloperidol"},
    {"smiles": "CC(CS(N)(=O)=O)NCC(O)c1ccc(N)cc1",
     "label": 1, "name": "Sotalol"},
    {"smiles": "CN(CCOc1ccc(NS(=O)(=O)C)cc1)CCOc1ccc(NS(=O)(=O)C)cc1",
     "label": 1, "name": "Dofetilide"},
    {"smiles": "COc1ccc2nccc(C(O)C3CC4CCN3CC4C=C)c2c1",
     "label": 1, "name": "Quinidine"},
]


def build_correction_set(weight: int = 15) -> pd.DataFrame:
    """Build the curated correction set with validation."""
    valid_records = []
    for entry in CURATED_CORRECTIONS:
        can_smi = canonicalize_smiles(entry["smiles"])
        if can_smi is None:
            logger.warning(f"Invalid SMILES for {entry['name']}, skipping")
            continue
        valid_records.append({
            "smiles": can_smi,
            "label": entry["label"],
            "source": "Curated",
        })
        logger.info(f"  ✓ {entry['name']}: valid (label={entry['label']})")

    correction_df = pd.DataFrame(valid_records)
    logger.info(f"Curated set: {len(correction_df)} valid molecules")

    # Weight by repeating
    weighted = pd.concat([correction_df] * weight, ignore_index=True)
    logger.info(f"Weighted correction set: {len(weighted)} rows (×{weight})")
    return weighted


# ─────────────────────────────────────────────────────────
# FINAL ASSEMBLY
# ─────────────────────────────────────────────────────────

def assemble_final_dataset(
    chembl_df: pd.DataFrame,
    tox21_df: pd.DataFrame,
    correction_weighted_df: pd.DataFrame,
    output_path: str = "data/herg_training_final.csv"
) -> pd.DataFrame:
    """Combine all sources, deduplicate, add corrections, shuffle."""

    base = pd.concat([chembl_df, tox21_df], ignore_index=True)
    base = deduplicate_df(base)

    final = pd.concat([base, correction_weighted_df], ignore_index=True)
    final = final.sample(frac=1, random_state=42).reset_index(drop=True)

    total = len(final)
    toxic = int(final["label"].sum())
    safe = int((final["label"] == 0).sum())
    unique = final["smiles"].nunique()

    logger.info("=" * 60)
    logger.info("FINAL DATASET SUMMARY")
    logger.info("=" * 60)
    logger.info(f"  Total samples:     {total}")
    logger.info(f"  Toxic (label=1):   {toxic} ({100*toxic/total:.1f}%)")
    logger.info(f"  Safe  (label=0):   {safe} ({100*safe/total:.1f}%)")
    logger.info(f"  Unique molecules:  {unique}")
    logger.info(f"  Sources:           ChEMBL, Tox21, Curated")
    logger.info("=" * 60)

    if total < 15000:
        logger.warning(
            f"Dataset below 15,000 target ({total}). "
            f"ChEMBL: {len(chembl_df)}, Tox21: {len(tox21_df)}"
        )

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    final.to_csv(output_path, index=False)
    logger.info(f"Saved to {output_path}")
    return final


# ─────────────────────────────────────────────────────────
# Legacy compatibility — keep compute_features for model.py
# ─────────────────────────────────────────────────────────

def compute_features(df: pd.DataFrame) -> tuple:
    """Compute RDKit descriptors. Kept for backward compatibility."""
    from rdkit import Chem
    from src.descriptors import compute_descriptors
    import numpy as np

    X_rows, y_rows = [], []
    failed = 0
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None:
            failed += 1
            continue
        try:
            desc = compute_descriptors(mol)
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
    logger.info(f"Feature extraction: {len(X_rows)} ok, {failed} failed")
    return np.array(X_rows, dtype=float), np.array(y_rows)


# ─────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────

if __name__ == "__main__":
    logger.info("Phase 1: Fetching ChEMBL multi-assay data...")
    chembl_df = fetch_chembl_multi_assay(max_per_target=8000)

    logger.info("Phase 2: Fetching Tox21 data...")
    tox21_df = fetch_tox21_herg()

    logger.info("Phase 3: Building curated correction set...")
    correction_df = build_correction_set(weight=15)

    logger.info("Phase 4: Assembling final dataset...")
    final_df = assemble_final_dataset(
        chembl_df, tox21_df, correction_df,
        output_path="data/herg_training_final.csv"
    )

    logger.info("Phase 5: Training Chemprop model...")
    try:
        from src.tox_model import train_and_save_model
        train_and_save_model(
            csv_path="data/herg_training_final.csv",
            model_dir="models/chemprop_herg",
            epochs=30
        )
        logger.info("Training complete. Restart server to load new model.")
    except Exception as e:
        logger.error(f"Training failed: {e}")
        import traceback
        traceback.print_exc()
        logger.info("Dataset saved. You can train manually later.")
