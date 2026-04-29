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

## 🚀 Key Features

### 1. Dual-Engine Generative Simulator
- **Liver Simulator (Hepatic Clearance):** Maps standard metabolic degradation via destructive Phase I (oxidation, dealkylation) and Phase II (glucuronidation) reaction rules.
- **MedChem Generator:** Explores local chemical space using constructive transformations (e.g., fluorination, methylation) to evolve parent drugs into safer analogs.
- **Dynamic Lipinski Filtering:** Silently drops generated analogs that exceed a maximum of one Rule of 5 violation (MW > 500, LogP > 5, HBD > 5, HBA > 10) to ensure oral bioavailability.

### 2. Interactive DAG Visualization (Metabolism Tree)
- **Real-Time Risk Mapping:** Every node in the generated trajectory is evaluated asynchronously by the MPNN, assigning a visual risk border (Green for <50% hERG risk, Red for ≥50% risk).
- **Shape-Coded Pathways:** UI distinguishes between Initial Molecules (Diamond), Phase I/Minor edits (Circle), and Phase II/Major edits (Hexagon).
- **Interactive Exploration:** Click any generated metabolite to set it as the new active drug and cascade a fresh analysis.

### 3. Dual Similarity Scoring & 3D Conformations
- **2D Similarity:** Tanimoto coefficient based on Morgan Fingerprints (Radius 2). 
- **3D Similarity:** Shape-based volumetric overlap using the Open3DAlign (O3A) algorithm.
- **3D Conformations:** RDKit's ETKDG v3 algorithm with MMFF94s force field optimization. 

### 4. Interactive Dashboard & Analytics
- **Live Visualizations:** 2D structure rendering via RDKit.js and 3D molecular alignment via 3Dmol.js. 
- **Physicochemical Profiling:** Radar and table views of 9 key descriptors (MW, Crippen logP, TPSA, HBD, HBA, Rotatable Bonds, Aromatic Rings, Fraction CSP3). 
- **Vercel Analytics:** Integrated tracking for production usage monitoring.

## 🛠️ Architecture

### Backend (FastAPI / Hugging Face Spaces)
- `api.py`: FastAPI entry point with CORS enabled for cross-origin frontend requests.
- `embedding.py`: Handles 3D conformation generation and energy minimization. 
- `similarity.py`: Computes 2D Tanimoto and 3D shape overlap. 
- `scoring.py`: Implements dual-scoring logic, DAG generation rules, and risk flagging. 
- `tox_model.py`: Handles the Chemprop Message Passing Neural Network (MPNN) inference.

### Frontend (React + Vite + Vercel)
- **State Management:** React Hooks mapped to dynamic backend endpoints.
- **Network Graphs:** `reactflow` paired with `dagre` for automated, hierarchical DAG layout generation.
- **Styling:** Tailwind CSS with Lucide-React iconography.
- **Deployment:** Hosted on Vercel with automatic CI/CD and `@vercel/analytics` integration.

## 📦 Installation & Setup

### 1. Backend Setup
```bash
cd backend
pip install -r requirements.txt
python main.py
```
The backend will start on `http://localhost:7860`.

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

$$
\text{Score} = \alpha \cdot \text{Tanimoto}_{2D} + \beta \cdot \text{Shape}_{3D}
$$

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

## 🔮 Roadmap & Future Work

While the current pipeline demonstrates robust accuracy for standard validation sets, we are actively developing the next iteration to address known edge cases and structural blindspots in the MPNN:

* **Physiological State Modeling (pH 7.4):** Integrating automated SMILES protonation (e.g., via Dimorphite-DL) prior to graph embedding. This will allow the model to recognize formal positive charges on basic nitrogens, which are critical for capturing the cation-π interactions that drive hERG channel blockade.
* **Continuous Risk Scoring (pIC50):** Transitioning the MPNN architecture from binary IC50 classification to continuous pIC50 regression. This will resolve the current "ambiguous zone" (1,000–10,000 nM) blindspot, providing more accurate, gradient-based risk assessments for borderline compounds like Verapamil.
* **3D-Aware Message Passing:** Upgrading the 2D topological graph convolution to incorporate 3D spatial coordinate embeddings. By feeding the ETKDG-generated conformations directly into the neural network, the model will better differentiate between geometrically distinct targets (e.g., L-type calcium channels vs. hERG) that share similar 2D lipophilic pharmacophores.
