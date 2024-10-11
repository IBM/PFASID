#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

from typing import Any
from .echa.classify import *
from .epa.classify import *
from .chemical.classify import *
from .ontochem.classify import *
from .OECD.classify import *
import multiprocessing

CLASSIFICATION_METHODS = {
    "ECHA_2021_struct": ECHA_2021_struct,
    # 'EPA_2021': EPA_2021,
    "contains_CF": F_only,
    "Ontochem_2022": ontochem_2022,
    "OECD_2021": OECD_2021,
    "EPA_2023": EPA_2023,
}

CLASSIFICATION_DATA = {
    # 'EPA_2021':{
    #     "description":"",
    #     "reference": r"https://pubs.acs.org/doi/10.1021/cen-10040-polcon1"},
    "contains_CF": {
        "description": "At least one C-F bond in the molecule.",
        "reference": "Dan Sanders, PhD",
    },
    "EPA_2023": {
        "description": "Under the finalized rule, “PFAS” were defined as any chemical that includes at least one of the following three structures: (1) R-(CF2)-CF(R’)R’’, where both the CF2 and CF moieties are saturated carbons; (2) R-CF2OCF2-R’, where R and R’ can either be F, O, or saturated carbons; and (3) CF3C(CF3)R’R’’, where R’ and R’’ can either be F or saturated carbons.",
        "reference": r"https://www.epa.gov/system/files/documents/2023-09/prepublicationcopy_7902-02_fr-doc_aa_esignatureverified_2023-09-28.pdf",
    },
    "OntoChem_2022": {
        "description": "CF3C(CF3)RR’ where R and R' are any atom other than H. Or -(CF2)-CF- bonded to any atoms.",
        "reference": r"https://pubs.rsc.org/en/content/articlelanding/2022/dd/d2dd00019a",
    },
    "OECD_2021": {
        "description": "Contains −CF3 or −CF2− without a H/Cl/Br/I atom attached to it.",
        "reference": r"https://pubs.acs.org/doi/10.1021/acs.est.1c06896",
    },
    "ECHA_2021_struct": {
        "description": "er- and polyfluoroalkyl substances (PFASs) defined as: Any substance that contains at least one fully fluorinated methyl (CF3-) or methylene (-CF2-) carbon atom (without any H/Cl/Br/I attached to it). A substance that only contains the following structural elements is excluded from the scope of the proposed restriction: CF3-X or X-CF2-X’, where X = -OR or -NRR’ and X’ = methyl (-CH3), methylene (-CH2-), an aromatic group, a carbonyl group (-C(O)-), -OR’’, -SR’’ or –NR’’R’’’, and where R/R’/R’’/R’’’ is a hydrogen (-H), methyl (-CH3), methylene (-CH2-), an aromatic group or a carbonyl group (-C(O)-).",
        "reference": r"https://echa.europa.eu/de/registry-of-restriction-intentions/-/dislist/details/0b0236e18663449b",
    },
}


def classify(
    smiles_list: list[str],
    classification_methods: list[str] = list(CLASSIFICATION_METHODS.keys()),
    multiprocesss: bool = True,
) -> dict[str, Any]:
    result = {}
    # ensure correct types
    assert isinstance(
        smiles_list, list
    ), f"All smiles must be listed in as strings in the format list[str]. You have given type {type(smiles_list)}"
    assert isinstance(
        classification_methods, list
    ), f"All methods must be listed in as strings in the format list[str]. You have given type {type(classification_methods)}"
    # check if smiles can be processed by rdkit
    good_smiles = []
    bad_smiles = []
    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            good_smiles.append(smi)
        except:
            bad_smiles.append(smi)
    # run one method at a time, on all smiles with multi processing
    for method in classification_methods:
        classification_function = CLASSIFICATION_METHODS.get(method, None)
        if multiprocesss == True:
            assert (
                classification_function is not None
            ), f"Classification method {method} not yet supported!"
            with multiprocessing.Pool() as pool:
                # apply function to all smiles
                multi_res = pool.map(classification_function, good_smiles)
                # format results
                result[method] = {
                    good_smiles[i]: multi_res[i] for i in range(len(good_smiles))
                }
        else:
            result[method] = {smi: classification_function(smi) for smi in good_smiles}
    # record bad smiles
    result["bad_smiles"] = bad_smiles
    return result


def get_classifiers():
    return CLASSIFICATION_DATA
