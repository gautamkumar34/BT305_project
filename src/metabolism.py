# src/metabolism.py
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski
from typing import Dict, Any, Set
from src.scoring import _get_tox_predictor

class MetabolismEngine:
    def __init__(self):
        self.liver_rules = [
            {"name": "Aromatic Hydroxylation", "rxn_smarts": "[c:1]-[H]>>[c:1]-[O]-[H]", "type": "Phase I"},
            {"name": "Aliphatic Hydroxylation", "rxn_smarts": "[C:1]-[H]>>[C:1]-[O]-[H]", "type": "Phase I"},
            {"name": "O-Dealkylation", "rxn_smarts": "[O:1]-[C:2]>>[O:1]-[H]", "type": "Phase I"},
            {"name": "Glucuronidation", "rxn_smarts": "[O:1]-[H]>>[O:1]-C1OC(C(=O)O)C(O)C(O)C1O", "type": "Phase II"}
        ]
        
        self.medchem_rules = [
            {"name": "Aromatic Fluorination", "rxn_smarts": "[c:1]-[H]>>[c:1]-[F]", "type": "Halogenation"},
            {"name": "Aromatic Chlorination", "rxn_smarts": "[c:1]-[H]>>[c:1]-[Cl]", "type": "Halogenation"},
            {"name": "N-Methylation", "rxn_smarts": "[N:1]-[H]>>[N:1]-[C]", "type": "Alkylation"},
            {"name": "O-Methylation", "rxn_smarts": "[O:1]-[H]>>[O:1]-[C]", "type": "Alkylation"},
        ]

        for rule in self.liver_rules + self.medchem_rules:
            rule["rxn"] = AllChem.ReactionFromSmarts(rule["rxn_smarts"])

    # --- PATH 2: Lipinski Filter Logic ---
    def passes_lipinski(self, mol) -> bool:
        mw = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        
        # Count violations
        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        return violations <= 1 # Allow max 1 violation

    # --- PATH 1: Toxicity Integration ---
    def get_toxicity_score(self, smiles: str) -> float:
        predictor = _get_tox_predictor()
        
        if predictor is not None:
            try:
                # Ask your Chemprop model for the real probability
                return float(predictor.predict_toxicity_prob(smiles))
            except Exception as e:
                # If RDKit generates a weird SMILES the model can't parse, fail gracefully
                print(f"Toxicity prediction failed for {smiles}: {e}")
                return 0.50 
        
        # Fallback if the model completely fails to load
        return 0.50
        
    def predict_metabolites(self, smiles: str, max_depth: int = 2, mode: str = "liver") -> Dict[str, Any]:
        parent_mol = Chem.MolFromSmiles(smiles)
        if not parent_mol:
            raise ValueError("Invalid SMILES provided")
        
        parent_mol = Chem.AddHs(parent_mol)
        root_name = "Parent Drug" if mode == "liver" else "Lead Compound"

        tree = {
            "smiles": smiles,
            "name": root_name,
            "transformation": "Initial Molecule",
            "toxicity": self.get_toxicity_score(smiles), # Fetch score for parent
            "children": []
        }

        parent_mol_clean = Chem.RemoveHs(parent_mol)
        for atom in parent_mol_clean.GetAtoms():
            atom.SetAtomMapNum(0)

        parent_smiles = Chem.MolToSmiles(parent_mol_clean, canonical=True)
        visited_smiles = {parent_smiles}
        
        active_rules = self.medchem_rules if mode == "medchem" else self.liver_rules
        self._generate_tree_recursive(parent_mol, tree, 0, max_depth, visited_smiles, active_rules, mode)
        return tree

    def _generate_tree_recursive(self, mol, current_node, depth, max_depth, visited_smiles: Set[str], active_rules, mode: str):
        if depth >= max_depth:
            return

        MAX_PRODUCTS_PER_RULE = 1

        for rule in active_rules:
            rxn = rule["rxn"]
            products = rxn.RunReactants((mol,))
            rule_success_count = 0
            
            for product_tuple in products:
                if rule_success_count >= MAX_PRODUCTS_PER_RULE:
                    break
                    
                prod_mol = product_tuple[0]
                try:
                    Chem.SanitizeMol(prod_mol)
                    prod_mol_no_h = Chem.RemoveHs(prod_mol)

                    # PATH 2: Generative Filter (Only apply to MedChem mode)
                    if mode == "medchem" and not self.passes_lipinski(prod_mol_no_h):
                        continue # Silently drop molecule if it violates Lipinski

                    for atom in prod_mol_no_h.GetAtoms():
                        atom.SetAtomMapNum(0)

                    new_smiles = Chem.MolToSmiles(prod_mol_no_h, canonical=True)
                    
                    if new_smiles and new_smiles not in visited_smiles:
                        visited_smiles.add(new_smiles)
                        rule_success_count += 1
                        
                        child_node = {
                            "smiles": new_smiles,
                            "name": rule["name"],
                            "transformation": rule["type"],
                            "toxicity": self.get_toxicity_score(new_smiles), # Fetch score for child
                            "children": []
                        }
                        current_node["children"].append(child_node)
                        
                        if depth + 1 < max_depth:
                            self._generate_tree_recursive(Chem.AddHs(prod_mol_no_h), child_node, depth + 1, max_depth, visited_smiles, active_rules, mode)
                            
                except Exception:
                    continue