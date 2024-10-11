#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

from rdkit import Chem
import numpy as np


def saturatedCCheck(mol, sub):
    # get all atoms which match the substructure
    subms_ats_tup = list(mol.GetSubstructMatches(sub))
    subms_ats_list = []
    for pfas_struct in subms_ats_tup:
        subms_ats_list.extend(list(pfas_struct))
    subms_ats = list(set(subms_ats_list))

    # get the neighbot count for each C atom in the substructure
    subms_neighbors = [
        mol.GetAtomWithIdx(int(i)).GetTotalNumHs()
        + len(mol.GetAtomWithIdx(int(i)).GetNeighbors())
        for i in subms_ats
        if mol.GetAtomWithIdx(int(i)).GetAtomicNum() == 6
    ]
    # Return true if all of the C atoms are saturated.
    subms_saturated = np.all([True if i == 4 else False for i in subms_neighbors])

    return subms_saturated


def getRNeighbors(mol, sub):
    # get all atoms which match the substructure
    subms_ats_tup = list(mol.GetSubstructMatches(sub))
    subms_ats_list = []
    for pfas_struct in subms_ats_tup:
        subms_ats_list.extend(list(pfas_struct))
    subms_ats = list(set(subms_ats_list))

    # get neighbors of carbon atoms within the subgroup
    subms_neighbors = []
    [
        subms_neighbors.extend(list(mol.GetAtomWithIdx(int(i)).GetNeighbors()))
        for i in subms_ats
        if mol.GetAtomWithIdx(int(i)).GetAtomicNum() == 6
    ]
    # isolate atom indicies which are only the R groups and not in the substructure
    subms_nonsub_neighbors = [
        a.GetIdx() for a in subms_neighbors if a.GetIdx() not in subms_ats
    ]
    return subms_nonsub_neighbors


def elementCheck(mol, atom_idx, symbols):
    # Get atom symbol for each atom
    at_symbols = {at: mol.GetAtomWithIdx(at).GetSymbol() for at in atom_idx}
    # Return if atoms are one of the element types in 'symbols'
    elements = {at: True if at_symbols[at] in symbols else False for at in atom_idx}
    return elements


def saturatedRCheck(mol, atom_idx):
    # count number of neighbors on an atom
    at_neighbors = {
        at: (
            mol.GetAtomWithIdx(at).GetTotalNumHs()
            + len(mol.GetAtomWithIdx(at).GetNeighbors())
            if mol.GetAtomWithIdx(at).GetSymbol() == "C"
            else 0
        )
        for at in atom_idx
    }
    # return if number of neighbors is 4
    elements = {at: True if at_neighbors[at] == 4 else False for at in atom_idx}
    return elements
