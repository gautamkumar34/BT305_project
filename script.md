# 📜 Project: BioTransformer

## 👥 Team Members

- **Gautam Kumar(230106026)**: Lead Developer & Cheminformatics Architect 
- **Dipanshu Subodh Ojha(230106024)**: Frontend Developer
- **Akshat Gupta(230106006)**: UI designer,Collaborative research and validation contributor

## 🎯 Objective

The primary objective of BioTransformer is to provide an end-to-end computational pipeline for **Cardiotoxicity Inference** and **Lead Optimization**. By integrating deep learning with classical 3D geometric analysis, the project aims to:

- **Predict hERG Liability**: Accurately classify drug-like molecules as toxic or safe based on their interaction with the human Ether-à-go-go-Related Gene (hERG) potassium channel.
- **Simulate Metabolic Fate**: Map the "Liver" clearance pathways to see how a drug degrades.
- **Propose Structural Rescue**: Generate "MedChem" analogs that retain therapeutic similarity while reducing cardiotoxic risk.

## 🛠️ Novel or Existing Tool?

**BioTransformer is a Hybrid Novel Pipeline.**

## 🏛️ Comparison with Existing Tools

| Feature | BioTransformer | Common Existing Tools (e.g., SwissADME, DataWarrior) |
|:---|:---|:---|
| **hERG Prediction** | MPNN-based (Chemprop v2) with deep electronic context. | Often use simpler Random Forest or SVM models based on 2D descriptors. |
| **3D Alignment** | Real-time 3Dmol.js rendering with O3A volumetric overlap. | Usually requires desktop software (PyMOL, Schrödinger) rather than a web-based UI. |
| **Metabolic DAG** | Interactive, generative Directed Acyclic Graph (DAG) with real-time risk scoring. | Metabolism predictors (like BioTransformer 3.0—the namesake) exist but are often static lists, not interactive "trees". |
| **MedChem Engine** | Constructive generative rules with Lipinski's Rule of 5 filters. | Requires separate "de novo" design software; rarely linked directly to a live hERG simulator. |

---

## 🛠️ Closest Existing Counterparts

### SwissADME / pkCSM
These are excellent for static physicochemical property calculation and binary toxicity flags. **However**, they lack:
- Interactive 3D alignment
- The generative medchem trajectory that allows you to click through a family tree of metabolites

### BioTransformer 3.0 (Commercial/Academic Tool)
There is a well-known academic tool also named "BioTransformer" which predicts metabolic products. **Our tool is distinct because**:
- It focuses specifically on the **cardiotoxicity (hERG) risk** of those products
- Provides a **visual, interactive network** for exploring metabolites
- Whereas the existing academic version is primarily a database of metabolic rules

### DataWarrior
This is a powerful open-source desktop app that can generate "evolutionary" chemical libraries. **Limitations compared to our tool**:
- Not a web-based tool
- Does not provide integrated 3D-aware similarity scoring (α and β weights)
- Lacks the real-time risk assessment pipeline

---

## 💎 The "BioTransformer" Novelty Gap

The novelty of our project lies in **Unified Interactive Orchestration**.

While a scientist could use three different tools to:
1. Check 3D alignment
2. Calculate LogP
3. Predict hERG risk

**Our tool does all three simultaneously** and then automatically generates safer versions of the molecule based on those results.

### The "Proactive Rescue" Capability

By providing the **"Liver"** and **"MedChem"** engines in a single DAG, we have moved the workflow from:

> ❌ **"Predicting a result"** → ✅ **"Navigating a solution"**

This **"proactive rescue"** capability is rarely found in free, web-based open-source tools.

---

## 📊 Summary: Why BioTransformer Stands Apart

| Aspect | Traditional Workflow | BioTransformer |
|:---|:---|:---|
| **Tool Chain** | 3-4 separate tools | Single unified platform |
| **3D Analysis** | Desktop software required | In-browser with 3Dmol.js |
| **Metabolite Exploration** | Static lists | Interactive DAG with click-through |
| **Lead Optimization** | Manual interpretation | Automatic "safer analog" generation |
| **Cost** | Often commercial | Free & open-source |
| **Platform** | Desktop-bound | Web-based + API accessible |

---

*"From prediction to solution — in one pipeline."*

## 💎 Unique Selling Proposition (USP)

What makes BioTransformer distinct from standard toxicity checkers is its **proactive optimization approach**:

- **Dynamic MPNN Evaluation**: Unlike static checkers, it evaluates an entire "family tree" of metabolites and analogs in seconds.
- **3D-Awareness**: It uses ETKDG v3 and Open3DAlign (O3A) to ensure that similarity isn't just "on paper" (2D) but reflects actual spatial overlap (3D) in the hERG pore.
- **Built-in Pharmacokinetic Filters**: It utilizes a strict Lipinski's Rule of 5 filter during generative mode to ensure that every suggested "safe" analog is actually bioavailable.
- **Visual Risk Mapping**: The use of red/green color-coded nodes in the trajectory allows for immediate visual identification of safe evolutionary paths for toxic lead compounds.

## 💻 Program Tools & Stack

The project is built using a modern, scalable stack designed for high-performance scientific computing:

### Backend (Python)

| Tool | Purpose |
|:---|:---|
| **FastAPI** | The asynchronous core for the API. |
| **RDKit** | Used for 2D/3D structure generation, fingerprinting, and reaction SMARTS. |
| **Chemprop (PyTorch)** | Powering the Message Passing Neural Network (MPNN). |
| **Dagre** | Used for computing the hierarchical layout of the metabolic trees. |

### Frontend (JavaScript/React)

| Tool | Purpose |
|:---|:---|
| **Vite + React** | The framework for a fast, responsive user interface. |
| **ReactFlow** | The library used to render the interactive metabolic networks. |
| **3Dmol.js** | For high-performance 3D molecular alignment rendering in the browser. |
| **RDKit.js (WASM)** | For client-side 2D structure rendering. |

### DevOps & Deployment

| Tool | Purpose |
|:---|:---|
| **Docker** | For standardizing the complex scientific environment. |
| **Hugging Face Spaces** | Hosting the backend API and deep learning models. |
| **Vercel** | Hosting the React frontend with integrated analytics. |

## 📂 Core Module Index

| Module | Description |
|:---|:---|
| `metabolism.py` | The generative engine for DAG trajectories. |
| `tox_model.py` | MPNN inference logic. |
| `descriptors.py` | Molecular property calculator. |
| `similarity.py` | 2D/3D similarity computation. |
| `embedding.py` | 3D conformer generation and minimization. |

## 📊 Key Metrics

| Metric | Value |
|:---|:---|
| **Training Data** | 21,000+ hERG assay records |
| **Model Architecture** | MPNN (Message Passing Neural Network) |
| **Validation ROC-AUC** | ≥0.80 |
| **Inference Latency** | <500ms (p95) |
| **Similarity Scoring** | Tanimoto (2D) + Shape Overlap (3D) |

## 🧪 Validation Case Study

| Compound | Known Outcome | Model Prediction |
|:---|:---|:---|
| **Terfenadine** | Cardiotoxic (withdrawn) | High risk (≥50%) |
| **Fexofenadine** | Safe (active metabolite) | Low risk (<50%) |

*The pipeline successfully differentiates the toxic parent from its safe metabolite using combined 2D/3D similarity + MPNN toxicity scoring.*

## 🔗 Live Demo

- **Frontend**: [bt-305-project.vercel.app](https://bt-305-project.vercel.app)
- **Backend API**: Hosted on Hugging Face Spaces
- **Repository**: [github.com/gautamkumar34/BT305_project](https://github.com/gautamkumar34/BT305_project)

---

*Last Updated: April 2026*