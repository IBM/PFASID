#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

from typing import Any
import numpy as np
from rdkit import Chem  # type: ignore


def F_only(smiles: str) -> bool:
    # convert smiles to mol object
    mol = Chem.MolFromSmiles(smiles)
    # Check C-F bond
    return mol.HasSubstructMatch(Chem.MolFromSmiles("CF"))
