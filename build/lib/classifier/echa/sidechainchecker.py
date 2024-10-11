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
        seen = [0] * mol.GetNumAtoms()
        for idx in vect:
            seen[idx] = 1
        # loop over the atoms we care about:
        for idx, qtyp in self._atsToExamine:
            aro_rings = getAromaticRings(mol)
            midx = vect[idx]
            if midx in aro_rings:
                return True
            else:
                return False


def getAromaticRings(mol):
    ri = mol.GetRingInfo()
    rings_aro = []
    for ring in ri.BondRings():
        if isRingAromatic(mol, ring):
            rings_aro.extend(list(ring))
    return rings_aro


def isRingAromatic(mol, bondRing):
    for id in bondRing:
        if not mol.GetBondWithIdx(id).GetIsAromatic():
            return False
    return True


def _aromatic_SubstructMatches(mol, sub):
    # this returns substructure matches for exceptions.
    # p = Chem.MolFromSmarts(sub)

    # Set property type search of "aro" for all *
    for atom in sub.GetAtoms():
        if atom.GetAtomicNum() == 0:
            sub.GetAtomWithIdx(atom.GetIdx()).SetProp("queryType", "aro")

    # Check for side chains
    params = Chem.SubstructMatchParameters()
    params.removeHs = False
    checker = SidechainChecker(sub)
    params.setExtraFinalCheck(checker)
    # format return
    if len(list(mol.GetSubstructMatches(sub, params))) > 0:
        return [i for i in mol.GetSubstructMatches(sub, params)]
    else:
        return [None]
