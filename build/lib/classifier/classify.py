#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

from typing import Any
from echa.classify import *
from epa.classify import *
import multiprocessing

CLASSIFICATION_METHODS = {
    "ECHA_2021_struct": ECHA_2021_struct,
    "EPA_2024_0": EPA_2024_0,
}

CLASSIFICATION_DATA = {
    "EPA_2024_0": {
        "description": "",
        "reference": r"https://pubs.acs.org/doi/10.1021/cen-10040-polcon1",
    },
    "EPA_2023_list": {
        "description": "",
        "reference": r"https://www.epa.gov/ccl/ccl-5-chemical-contaminants",
    },
    "ECHA_2021_struct": {
        "description": "",
        "reference": r"https://echa.europa.eu/de/registry-of-restriction-intentions/-/dislist/details/0b0236e18663449b",
    },
}


def classify(
    smiles_list: list[str],
    classification_methods: list[str] = list(CLASSIFICATION_METHODS.keys()),
) -> dict[str, Any]:
    result = {}

    for method in classification_methods:
        classification_function = CLASSIFICATION_METHODS.get(method, None)
        assert (
            classification_function is not None
        ), f"Classification method {method} not yet supported!"
        with multiprocessing.Pool() as pool:
            multi_res = pool.map(classification_function, smiles_list)
            result[method] = {
                smiles_list[i]: multi_res[i] for i in range(len(smiles_list))
            }

    return result


def get_classifiers():
    return CLASSIFICATION_DATA


# if __name__ == '__main__':
#     test_smiles = ['CC(C(=O)O)c1ccc(C(F)(F)N(C)C)cc1','FC1(F)N=C2N(c3c4ccccc4c(N4CCCN5CCCN=C54)c4ccccc34)CCCN2C(F)(F)C1(F)F','O=S1(=O)OC(F)(F)C(F)=C(F)C(F)(F)O1','CC(C(=O)O)c1ccc(C(F)(F)N(C)C)cc1','CNC(F)(F)F','CC(=O)OC(F)(F)c1ccccc1','CC(C(=O)O)c1ccc(C(F)(F)N(C)C)cc1','CC(=O)OC(F)(F)c1ccccc1','CCC(F)(F)F','CC(=O)C(F)(F)c1ccccc1','CC(=O)OC(F)(F)c1ccccc1','OC(F)(F)c1ccccc1','CC(F)(F)OC(F)(F)CCC(C(F)(F)(F))(C(F)(F)(F))CCC(C)(F)C(F)(F)CC','CNC(F)(F)F','COC(F)(F)F','OC(F)(F)F']
#     print(get_classifiers())
#     print(classify(test_smiles,['EPA_2024_0','ECHA_2021_struct']))
