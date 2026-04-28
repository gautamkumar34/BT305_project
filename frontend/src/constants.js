export const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:8000'

export const ENDPOINTS = {
  HEALTH: '/health',
  SIMILARITY: '/compare',
  VALIDATE: '/validate',
}

export const DEFAULT_CONFIG = {
  alpha: 0.5,
  beta: 0.5,
}

export const RISK_LEVELS = {
  LOW_RISK: { label: 'Low', color: 'green' },
  MODERATE_RISK: { label: 'Moderate', color: 'yellow' },
  HIGH_RISK: { label: 'High', color: 'red' },
}

export const TEST_CASES = [
  {
    id: 'terfenadine_fexofenadine',
    label: 'Terfenadine vs Fexofenadine (hERG classic)',
    smiles_a:
      'CC(C)(C)c1ccc(C(O)CCCN2CCC(CC2)C(O)(c2ccccc2)c2ccccc2)cc1',
    smiles_b:
      'CC(C)(C(=O)O)c1ccc(C(O)CCCN2CCC(CC2)C(O)(c2ccccc2)c2ccccc2)cc1',
  },
  {
    id: 'cisapride_domperidone',
    label: 'Cisapride vs Domperidone',
    smiles_a: 'OC(CCN1CCC(CC1)Nc1ncc(Cl)cc1OC)c1ccc(F)cc1',
    smiles_b: 'O=C(CCCN1CCC(CC1)n1c(=O)[nH]c2ccccc21)c1ccc(Cl)cc1',
  },
  {
    id: 'aspirin_ibuprofen',
    label: 'Aspirin vs Ibuprofen',
    smiles_a: 'CC(=O)Oc1ccccc1C(=O)O',
    smiles_b: 'CC(C)Cc1ccc(C(C)C(=O)O)cc1',
  },
  {
    id: 'astemizole_loratadine',
    label: 'hERG Litmus Test (Astemizole vs. Loratadine)',
    smiles_a: 'COc1ccc(CCN(C)C2CCN(Cc3nc4ccccc4[nH]3)CC2)cc1',
    smiles_b: 'CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc4cccnc24)CC1',
  },
]
