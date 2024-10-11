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
    mol = Chem.MolFromSmiles(smiles)
    # here we need to see if the fluorinated carbons are part
    # of the "exceptions" in ECHA regulations
    mol = Chem.AddHs(mol)
    echa_pfas_found: bool = False

    # get atom numbers if they are part of a PFAS group on the
    # molecule of interest
    CF2_ats = list(mol.GetSubstructMatches(Chem.MolFromSmiles("C(F)(F)")))

    # FOR EACH SUBSTRUCTURE SEE IF IT HAS AN EXCEPTION
    subats_aromatic_ragged = []
    for pfas_aro_sub in ECHA_2021_SUB_ARO:
        subats_aromatic_ragged.extend(_aromatic_SubstructMatches(mol, pfas_aro_sub))
    subats_aromatic = [i for i in subats_aromatic_ragged if i is not None]

    # find exception substructure matches for the molecule
    CF2_subs = np.array(
        [pfas_sub for pfas_sub in ECHA_2021_SUB if mol.HasSubstructMatch(pfas_sub)]
    )

    # if there was an exception by any substructure...
    if len(CF2_subs) > 0 or len(subats_aromatic) > 0:
        # find the atoms on the original molecules which are in the
        # "exceptions" groups
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

    if (echa_pfas_found == False) and (len(CF2_ats) > 0):
        return True
        # no exceptions, and at least 1 CF2 group
        # return "True" if there are no subgroups found from the exception list
        # CF3_ats = list(mol.GetSubstructMatches(
        #     Chem.MolFromSmiles('C(F)(F)(F)')))
        # if 3*len(CF3_ats) < len(CF2_ats) and len(CF3_ats) > 0:
        #     return (True, ['C(F)(F)(F)', 'C(F)(F)'])
        # elif 3*len(CF3_ats) == len(CF2_ats) and len(CF3_ats) > 0:
        #     return (True, ['C(F)(F)(F)'])
        # else:
        #     return (True, ['C(F)(F)'])
    else:
        # return (False, [])
        return False
