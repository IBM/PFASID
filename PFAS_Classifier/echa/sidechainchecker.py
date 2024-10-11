#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

from rdkit import Chem  # type: ignore


class SidechainChecker(object):
    matchers = {"aro": lambda at: at.GetIsAromatic()}

    def __init__(self, query, pName="queryType"):
        # identify the atoms that have the properties we care about
        self._atsToExamine = [
            (x.GetIdx(), x.GetProp(pName)) for x in query.GetAtoms() if x.HasProp(pName)
        ]
        self._pName = pName

    def __call__(self, mol, vect):
        aro_rings = getAromaticRings(mol)
        # loop over the atoms we care about:
        for idx, qtyp in self._atsToExamine:
            # idx is the atom in the substructure
            # get index of the matching atom in the mol
            midx = vect[idx]
            if midx in aro_rings:
                return True
            else:
                return False


def getAromaticRings(mol):
    ri = mol.GetRingInfo()
    rings_aro = []
    for ring in ri.AtomRings():
        if isRingAromatic(mol, ring):
            rings_aro.extend(list(ring))
    return rings_aro


def isRingAromatic(mol, AtomRing):
    for id in AtomRing:
        if not mol.GetAtomWithIdx(id).GetIsAromatic():
            return False
    return True


def _aromatic_SubstructMatches(mol, sub):
    for atom in sub.GetAtoms():
        if atom.GetAtomicNum() == 0:
            # label substructure atoms with * as 'aromatic'
            sub.GetAtomWithIdx(atom.GetIdx()).SetProp("queryType", "aro")

    # Check for aromatic side chains
    params = Chem.SubstructMatchParameters()
    params.removeHs = False
    checker = SidechainChecker(sub)
    params.setExtraFinalCheck(checker)

    # return substructure matches if there are any
    if len(list(mol.GetSubstructMatches(sub, params))) > 0:
        return [i for i in mol.GetSubstructMatches(sub, params)]
    else:
        return [None]
