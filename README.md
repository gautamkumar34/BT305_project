---
title: BioTransformer API
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: docker
pinned: false
app_port: 7860
---
# BioTransformer: Ligand-Based Drug Similarity & Cardiotoxicity Pipeline

BioTransformer is a production-grade computational chemistry pipeline designed for ligand-based drug similarity analysis and cardiotoxicity validation. It combines 2D topological fingerprints with 3D geometric shape overlap to provide a robust estimation of molecular similarity.

## 🚀 Features

- **3D Conformation Generation**: Uses RDKit's ETKDG algorithm with UFF/MMFF optimization to generate low-energy 3D structures.
- **Dual Similarity Scoring**:
  - **2D Similarity**: Tanimoto coefficient based on Morgan Fingerprints (Radius 2).
  - **3D Similarity**: Shape-based overlap scoring using RDKit alignment tools.
- **Configurable Scoring**: Adjustable weights ($\alpha$ and $\beta$) to balance topological vs. geometric similarity.
- **Cardiotoxicity Validation**: Built-in validation for known clinical cases (e.g., Terfenadine vs. Fexofenadine).
- **Interactive Dashboard**: 
  - 2D structure rendering via RDKit.js.
  - 3D molecular alignment visualization via 3Dmol.js.
  - Physicochemical property profiling (MW, LogP, HBD, HBA, TPSA).
  - AI-driven interpretability and risk assessment.

## 🛠️ Architecture

### Backend (FastAPI)
- `embedding.py`: Handles 3D conformation generation and energy minimization.
- `similarity.py`: Computes 2D Tanimoto and 3D shape overlap.
- `scoring.py`: Implements the dual-scoring logic and risk flagging.
- `validation.py`: Manages known toxicity datasets for pipeline verification.

### Frontend (React + Vite)
- **State Management**: React Hooks for real-time analysis updates.
- **Visualization**: RDKit.js (WASM) for 2D and 3Dmol.js for 3D rendering.
- **Styling**: Tailwind CSS for a professional, scientific dashboard.

## 📦 Installation & Setup

### 1. Backend Setup
```bash
cd backend
pip install -r requirements.txt
python main.py
```
The backend will start on `http://localhost:8000`.

### 2. Frontend Setup
```bash
cd frontend
npm install
npm run dev
```
The frontend will be available at `http://localhost:5173`.

## 🧪 Test Case: Terfenadine vs Fexofenadine

The system is validated using the Terfenadine case:
- **Terfenadine**: Known to cause QT prolongation and cardiotoxicity.
- **Fexofenadine**: The active metabolite, designed to be safer and non-toxic.

**Expected Result**: The pipeline should show moderate 2D similarity but distinct 3D overlap differences, correctly flagging the toxicity risk differentiation.

## 🔬 Scientific Methodology

The final similarity score is calculated as:
$$Score = \alpha \cdot \text{Tanimoto}_{2D} + \beta \cdot \text{Shape}_{3D}$$

Where:
- $\text{Tanimoto}_{2D}$ captures the presence of specific chemical substructures.
- $\text{Shape}_{3D}$ captures the volumetric overlap of the lowest-energy conformers, reflecting how the molecule actually fits into a biological target (e.g., hERG channel).

## 🧠 Toxicity Prediction Model

### Architecture
Message Passing Neural Network (MPNN) implemented via
[Chemprop v2](https://github.com/chemprop/chemprop) (Yang et al., J Chem Inf Model, 2019).

Unlike fingerprint-based models, Chemprop performs
neighborhood aggregation over the molecular graph,
capturing atom-level electronic environments including:
- **Nitrogen basicity** (critical for hERG channel binding)
- **Electron-withdrawing group effects** (modulate cation-π interactions)
- **Ring system aromaticity** (π-stacking in the hERG inner cavity)

### Training Data
Combined dataset from three sources:

| Source | Records | Notes |
|--------|---------|-------|
| ChEMBL CHEMBL240 | ~8,000 | Primary hERG assay |
| ChEMBL CHEMBL4282 | ~5,000 | Patch clamp |
| ChEMBL CHEMBL3307 | ~3,000 | IKr potassium channel |
| Tox21 NR-hERG | ~8,000 | NIH screening |
| Curated ground truth | 15 × 15 | Clinical verification |
| **Total** | **~21,000** | After dedup |

**IC50 threshold:** ≤1,000 nM = cardiotoxic (label=1)  
**Ambiguous zone** (1,000–10,000 nM) **excluded** from training.

### Validation
- 80/20 stratified train/validation split
- Target validation ROC-AUC: ≥0.80
- Litmus test: Terfenadine (toxic) vs Fexofenadine (safe) must always be correctly ranked

### Retraining
```bash
# Fetch data + train model
python scripts/fetch_herg_data.py

# Restart server to load new model
python main.py serve
```

### Limitations
- Trained on in vitro IC50 data (not in vivo cardiac outcomes)
- Does not model active metabolites
- Neutral-form SMILES used (not protonated at pH 7.4)
- Not validated for novel scaffolds outside drug-like chemical space
- **NOT for clinical use**

### References
1. Yang et al. *Analyzing Learned Molecular Representations for Property Prediction.* J Chem Inf Model. 2019;59(8):3370-3388.
2. Sanguinetti & Tristani-Firouzi. *hERG potassium channels and cardiac arrhythmia.* Nature. 2006;440:463-469.
3. [ChEMBL database](https://www.ebi.ac.uk/chembl/)
4. [Tox21 Challenge](https://tripod.nih.gov/tox21/challenge/)
