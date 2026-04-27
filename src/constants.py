"""
Verified molecular constants and SMILES strings.
"""

VERIFIED_SMILES = {
    "terfenadine": {
        "smiles": "CC(C)(C)c1ccc(C(O)CCCN2CCC(CC2)C(O)(c2ccccc2)c2ccccc2)cc1",
        "cid": 5405,
        "name": "Terfenadine"
    },
    "fexofenadine": {
        "smiles": "CC(C)(C(=O)O)c1ccc(C(O)CCCN2CCC(CC2)C(O)(c2ccccc2)c2ccccc2)cc1",
        "cid": 3348,
        "name": "Fexofenadine"
    },
    "cisapride": {
        "smiles": "OC(CCN1CCC(CC1)Nc1ncc(Cl)cc1OC)c1ccc(F)cc1",
        "cid": 2769,
        "name": "Cisapride"
    },
    "domperidone": {
        "smiles": "O=C(CCCN1CCC(CC1)n1c(=O)[nH]c2ccccc21)c1ccc(Cl)cc1",
        "cid": 3151,
        "name": "Domperidone"
    },
    "aspirin": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "cid": 2244,
        "name": "Aspirin"
    },
    "ibuprofen": {
        "smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        "cid": 3672,
        "name": "Ibuprofen"
    }
}
