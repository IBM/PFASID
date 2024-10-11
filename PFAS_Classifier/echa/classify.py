#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

from typing import Any
import numpy as np
from rdkit import Chem  # type: ignore
from .sidechainchecker import _aromatic_SubstructMatches
from .substructures import ECHA_2021_SUB, ECHA_2021_SUB_ARO


def ECHA_2021_struct(smiles: str) -> bool:
    # convert smiles to mol object
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    # initialize response in case there are no exceptions found
    echa_pfas_found: bool = False

    # find all CF2 groups, CF2 will be a part of each possible exception group
    CF2_ats = list(mol.GetSubstructMatches(Chem.MolFromSmiles("C(F)(F)")))

    # check if aromatic exception subgroups are found in the mol
    subats_aromatic_ragged = []
    for pfas_aro_sub in ECHA_2021_SUB_ARO:
        subats_aromatic_ragged.extend(_aromatic_SubstructMatches(mol, pfas_aro_sub))
    subats_aromatic = [i for i in subats_aromatic_ragged if i is not None]

    # check if non-aromatic exception subgroups are found in the mol
    CF2_subs = np.array(
        [pfas_sub for pfas_sub in ECHA_2021_SUB if mol.HasSubstructMatch(pfas_sub)]
    )

    # if there was an exception by any substructure...
    # check if all cf2 groups are within the exception groups
    if len(CF2_subs) > 0 or len(subats_aromatic) > 0:

        subats_ragged = [list(mol.GetSubstructMatches(cf2_sub)) for cf2_sub in CF2_subs]
        subats = [ats for tup in subats_ragged for ats in tup]
        subats_all = list(set(subats) | set(subats_aromatic))

        echa_pfas_found = bool(
            np.all(
                [
                    np.any([np.all([atom in i for atom in j]) for i in subats_all])
                    for j in CF2_ats
                ]
            )
        )
    # if there were no cf2 groups found outside of exception groups and there are cf2 groups on the molecule
    if (echa_pfas_found == False) and (len(CF2_ats) > 0):
        return True
    else:
        return False
