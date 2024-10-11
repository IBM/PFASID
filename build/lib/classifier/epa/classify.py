#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

from typing import Any
import numpy as np
from rdkit import Chem  # type: ignore


def EPA_2024_0(smiles: str) -> bool:
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
        # pfas_epa_patterns = list(np.array(EPA_SUBS, dtype=object)[subms])
        return True
    else:
        return False
