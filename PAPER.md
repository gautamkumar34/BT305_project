# Ligand-Based 3D Molecular Similarity for Cardiotoxicity Inference: A Cheminformatics Pipeline

## 1. Introduction
In early-stage drug discovery, identifying potential cardiotoxicity—specifically the blockade of the hERG (human Ether-à-go-go-Related Gene) potassium channel—is critical to avoid late-stage clinical failure. hERG blockade can lead to QT interval prolongation, increasing the risk of fatal ventricular arrhythmias. While high-throughput screening and *in vitro* patch-clamp assays are gold standards, there is a significant need for fast, *in silico* screening tools that can predict liability based on molecular structure. This work presents a production-grade pipeline that combines 2D structural fingerprints and 3D geometric similarity to infer cardiotoxicity risk.

## 2. Methods

### 2.1 Conformer Generation
To capture the spatial requirements of the hERG pore, 3D conformations were generated using the ETKDG v3 algorithm. For each input SMILES, 20 conformers were generated and optimized using the MMFF94s force field with a fixed seed (42) for reproducibility. The conformer with the lowest energy was selected for subsequent similarity analysis.

### 2.2 Similarity Metrics
A dual-scoring approach was implemented:
- **2D Similarity**: Computed using Morgan fingerprints (radius=2) and the Tanimoto coefficient, capturing topological connectivity.
- **3D Similarity**: Implemented via the Open3DAlign (O3A) algorithm, which performs a shape-based alignment and computes a Shape Tanimoto score based on volumetric overlap.

### 2.3 Descriptor Engineering
Nine physicochemical descriptors were extracted, including Molecular Weight (MW), Crippen logP, Topological Polar Surface Area (TPSA), Hydrogen Bond Donors (HBD), Hydrogen Bond Acceptors (HBA), Rotatable Bonds, Aromatic Rings, and Fraction CSP3. It is noted that Crippen logP represents the neutral form and may overestimate lipophilicity for ionizable groups at physiological pH.

### 2.4 ML Model (Toxicity Prediction)
A Message Passing Neural Network (MPNN) implemented via Chemprop v2 was trained on a comprehensive dataset aggregated from ChEMBL and Tox21 (~21,000 total records). Molecules with an IC50 ≤ 1,000 nM were classified as cardiotoxic, while those with IC50 > 10,000 nM were deemed safe (the ambiguous intermediate zone was excluded). The graph-based neural network performs neighborhood aggregation to capture critical atom-level electronic environments, such as nitrogen basicity and electron-withdrawing group effects, replacing legacy rule-based heuristics.

### 2.5 Generative Trajectory Simulation
To transition from passive evaluation to proactive structural optimization, a dual-engine deterministic Directed Acyclic Graph (DAG) simulator was engineered:
- **Hepatic Clearance Simulation (Liver Engine):** Applies destructive Phase I (e.g., oxidation, dealkylation) and Phase II (e.g., glucuronidation) SMARTS-based reaction rules to map standard metabolic degradation.
- **Generative Multi-Parameter Optimization (MedChem Engine):** Applies constructive transformations (e.g., fluorination, methylation) to evolve the parent drug. 
- **Dynamic Property Filtering:** The MedChem engine utilizes a strict Lipinski's Rule of 5 filter, silently dropping generated analogs that exceed a maximum of one violation (MW > 500, LogP > 5, HBD > 5, HBA > 10) to ensure oral bioavailability.
- **Real-time Toxicity Inference:** Every generated node in the trajectory is evaluated asynchronously by the MPNN, assigning a visual risk border (Green for <50% hERG risk, Red for ≥50% risk) to instantly highlight safe evolutionary pathways.

## 3. Results

### 3.1 Static Pair Validation
The pipeline was validated against four distinct molecular pairs. The results are summarized in Table 1.

**Table 1: Validation Results**
| Case | 2D Tanimoto | Shape 3D | Risk A | Risk B | ML Correct |
|------|-------------|----------|--------|--------|------------|
| Terfenadine/Fexofenadine | 0.804 | 0.323 | HIGH | LOW | ✓ |
| Cisapride/Domperidone    | 0.207 | 0.295 | MODERATE | HIGH | ✓ |
| Aspirin/Ibuprofen        | 0.195 | 0.435 | LOW | LOW | ✓ |
| Astemizole/Loratadine    | 0.120 | 0.000 | HIGH | LOW | ✓ |

For the Terfenadine/Fexofenadine pair, the system correctly identified a high 2D similarity (0.804) but a significantly lower 3D overlap (0.323). The transition from Terfenadine to Fexofenadine involved a significant TPSA increase (+37.3 Å²) and a logP decrease (-0.94), which the pipeline correctly associated with a reduction in predicted membrane permeability and hERG risk. 

Furthermore, the MPNN model demonstrated the ability to learn complex electronic effects by correctly identifying Astemizole as highly toxic (>98% probability) while predicting its structural relative Loratadine, which features a deactivated carbamate nitrogen, as safe (<6% probability).

### 3.2 Generative Pathway Validation
When initialized with toxic parent compounds (e.g., Astemizole), the Generative Trajectory Simulator successfully mapped functional chemical modifications. The integration of dynamic MPNN scoring across the DAG allowed for the instantaneous visual identification of "safe" analogs (green-bordered nodes) amidst predominantly toxic modifications, proving the utility of combining generative rules with real-time deep learning evaluation.

## 4. Discussion
The results demonstrate that 2D similarity alone is insufficient for toxicity inference; the Terfenadine case shows that a small structural change (addition of a carboxylic acid) drastically alters the safety profile despite high topological similarity. The MPNN successfully captured non-linear electronic features that basic lipophilicity descriptors miss.

Limitations include the reliance on neutral-form SMILES and in vitro IC50 training data rather than in vivo clinical outcomes. Unlike docking approaches, this pipeline is purely ligand-based and does not utilize the hERG protein structure. However, its computational efficiency makes it ideal for initial library screening.

## 5. Conclusion
We implemented a robust cheminformatics pipeline that successfully differentiates between cardiotoxic and safe analogs. The integration of 3D shape analysis and a highly-trained graph neural network (MPNN) provides a powerful, fast, and scalable method for hERG liability prediction.

## References
1. Landrum G. RDKit: Open-source cheminformatics. https://www.rdkit.org
2. Riniker S, Landrum GA. J Chem Inf Model. 2015;55(12):2562-74.
3. Tosco P, Stiefl N, Landrum G. J Cheminform. 2014;6:26.
4. Aronov AM. Drug Discov Today. 2008;13(5-6):149-55.
5. Sanguinetti MC, Tristani-Firouzi M. Nature. 2006;440:463-9.
