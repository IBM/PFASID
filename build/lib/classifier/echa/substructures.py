#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

import json
import os
from rdkit import Chem  # type: ignore

# defining and loading constant parameters

dir_path = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(dir_path, "ECHA_2021_struct_subs.json"), "r") as f:
    echa_subs = json.load(f)

params_echa_smiles = Chem.SmilesParserParams()
params_echa_smiles.removeHs = False
ECHA_2021_SUB = [
    Chem.MolFromSmiles(sub, params_echa_smiles) for sub in echa_subs["standard"]
]
ECHA_2021_SUB_ARO = [Chem.MolFromSmarts(sub) for sub in echa_subs["aromatic"]]
