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

