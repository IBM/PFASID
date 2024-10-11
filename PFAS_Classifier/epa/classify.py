#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

from typing import Any
import numpy as np
from rdkit import Chem  # type: ignore
from .utils import *


# retired/incorrect
def EPA_2021(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    EPA_SUBS = ["C(F)(F)C(F)", "C(F)(F)OC(F)(F)", "C(C(F)(F)(F))(C(F)(F)(F))"]
    EPA_PFAS = [Chem.MolFromSmiles(p) for p in EPA_SUBS]
    # substructure matches for CF2
    subms = np.array([mol.HasSubstructMatch(p) for p in EPA_PFAS])

    subms_ats_tup = []
    for p in np.array(EPA_PFAS)[subms]:
        subms_ats_tup.extend(list(mol.GetSubstructMatches(p)))

    subms_ats_list = []
    for pfas_struct in subms_ats_tup:
        subms_ats_list.extend(list(pfas_struct))

    subms_ats = list(set(subms_ats_list))

    # Do any substructures have Hydrogen
    # (which would make the substructure NOT PFAS)
    subms_Hs = [
        mol.GetAtomWithIdx(int(i)).GetNumImplicitHs()
        for i in subms_ats
        if mol.GetAtomWithIdx(int(i)).GetAtomicNum() == 6
    ]

    # if there are some Hs, the substructure is not PFAS
    if (~np.any(subms_Hs)) & (np.any(subms)):
        return True
    else:
        return False


def EPA_2023(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)

    EPA_1_sub = Chem.MolFromSmiles("C(F)(F)C(F)")
    # check that both Cs have are saturated
    EPA_1_pass = np.all(
        [mol.HasSubstructMatch(EPA_1_sub), saturatedCCheck(mol, EPA_1_sub)]
    )

    EPA_2_sub = Chem.MolFromSmiles("C(F)(F)OC(F)(F)")
    # return atom idx of the R groups attached to the Cs in the substructure
    EPA_2_neighbors = getRNeighbors(mol, EPA_2_sub)
    # Are the R groups saturated carbon
    EPA_2_R_saturated = saturatedRCheck(mol, EPA_2_neighbors)
    # Are the R groups F or O
    EPA_2_R_elements = elementCheck(mol, EPA_2_neighbors, ["F", "O"])
    # R groups pass at least 1/2 tests
    EPA_2_R_group_test = np.all(
        [np.any([EPA_2_R_elements[i], EPA_2_R_saturated[i]]) for i in EPA_2_neighbors]
    )
    EPA_2_pass = np.all([mol.HasSubstructMatch(EPA_2_sub), EPA_2_R_group_test])

    EPA_3_sub = Chem.MolFromSmiles("C(C(F)(F)(F))(C(F)(F)(F))")
    # return atom idx of the R groups
    EPA_3_neighbors = getRNeighbors(mol, EPA_3_sub)
    # Are the R groups saturated carbon
    EPA_3_R_saturated = saturatedRCheck(mol, EPA_3_neighbors)
    # Are the R groups F or O
    EPA_3_R_elements = elementCheck(mol, EPA_3_neighbors, ["F"])
    # R groups pass at least 1/2 tests
    EPA_3_R_group_test = np.all(
        [np.any([EPA_3_R_elements[i], EPA_3_R_saturated[i]]) for i in EPA_3_neighbors]
    )
    EPA_3_pass = np.all([mol.HasSubstructMatch(EPA_3_sub), EPA_3_R_group_test])

    # if any of the definitions are met...
    if np.any([EPA_1_pass, EPA_2_pass, EPA_3_pass]):
        return True
    else:
        return False
