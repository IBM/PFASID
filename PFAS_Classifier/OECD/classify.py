#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

from typing import Any
import numpy as np
from rdkit import Chem  # type: ignore


def OECD_2021(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    OECD_SUBS = ["C(F)(F)(F)", "C(F)(F)"]
    OECD_PFAS = [Chem.MolFromSmiles(p) for p in OECD_SUBS]
    # are there substructure matches for the two sub groups
    subms = np.array([mol.HasSubstructMatch(p) for p in OECD_PFAS])

    # get the atoms for the substructure matches
    subms_ats_tup = []
    for p in np.array(OECD_PFAS)[subms]:
        subms_ats_tup.extend(list(mol.GetSubstructMatches(p)))
    subms_ats_list = []
    for pfas_struct in subms_ats_tup:
        subms_ats_list.extend(list(pfas_struct))
    subms_ats = list(set(subms_ats_list))

    # Do any substructures have H/Cl/Br/I (which would make the substructure NOT PFAS)
    halogens = ["H", "Cl", "Br", "I"]
    subms_neighbors = []
    # get all of the neighbors of the Cs as the atom symbols
    [
        subms_neighbors.extend(
            [
                neighbor.GetSymbol()
                for neighbor in mol.GetAtomWithIdx(int(i)).GetNeighbors()
            ]
        )
        for i in subms_ats
        if mol.GetAtomWithIdx(int(i)).GetAtomicNum() == 6
    ]
    # return if the symbol is in the halogens list
    subms_halogen = [i for i in subms_neighbors if i in halogens]

    # if there are some no halogens, the substructure is PFAS
    if (len(subms_halogen) == 0) & (np.any(subms)):
        return True
    else:
        return False
