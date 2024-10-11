#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

from typing import Any
import numpy as np
from rdkit import Chem  # type: ignore


def ontochem_2022(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    # substructure match for 1st subgroup
    small_sub = Chem.MolFromSmiles("C(C(F)(F)(F))(C(F)(F)(F))")
    small_subms = mol.HasSubstructMatch(small_sub)
    small_subms_ats_tup = list(mol.GetSubstructMatches(small_sub))
    # get stoms for the substructure
    small_subms_ats_list = []
    for pfas_struct in small_subms_ats_tup:
        small_subms_ats_list.extend(list(pfas_struct))
    small_subms_ats = list(set(small_subms_ats_list))

    # Do any substructure matches have Hydrogen (which would make the substructure NOT PFAS)
    small_subms_Hs = [
        mol.GetAtomWithIdx(int(i)).GetNumImplicitHs()
        for i in small_subms_ats
        if mol.GetAtomWithIdx(int(i)).GetAtomicNum() == 6
    ]

    # check 2nd subgroup match, any matches would be a PFAS
    big_sub = Chem.MolFromSmiles("C(F)(F)CF")
    big_sub_match = mol.HasSubstructMatch(big_sub)

    # if either substructure matches, it is PFAS
    if ((~np.any(small_subms_Hs)) & (small_subms)) | (big_sub_match):
        return True
    else:
        return False
