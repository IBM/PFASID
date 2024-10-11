#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

from PFAS_Classifier import (
    classify,
    CLASSIFICATION_METHODS,
    CLASSIFICATION_DATA,
    get_classifiers,
)
import pandas as pd
from typing import Any
import multiprocessing

TEST_SMILES = [
    "CC(C(=O)O)c1ccc(C(F)(F)N(C)C)cc1",
    "CC(C)(C(F)(F)(F))(C(F)(F)(F))",
    "CC1=NC=C(COP(O)(O)=O)C(CN)=C1O",
]


def check_classification(smiles, classification, classification_methods):

    assert isinstance(classification, dict)
    assert len(classification_methods) == len(classification.keys()) - 1
    for method in classification_methods:
        for smile in smiles:
            print(f"Classified {smile} as {classification[method][smile]}")
            assert classification.get(method, None) is not None
            assert (
                classification[method][smile] is True
                or classification[method][smile] is False
            )


def check_classification_with_all_methods(smiles, classification):
    check_classification(smiles, classification, list(CLASSIFICATION_METHODS.keys()))


def check_classification_with_one_method(smiles, classification, method):
    check_classification(smiles, classification, [method])


def test_classify_all_implicit() -> None:
    classification = classify(TEST_SMILES)
    check_classification_with_all_methods(TEST_SMILES, classification)


def test_classify_all_explicit() -> None:
    classification = classify(TEST_SMILES, list(CLASSIFICATION_METHODS.keys()))
    check_classification_with_all_methods(TEST_SMILES, classification)


def test_classify_only_first_method() -> None:
    method = list(CLASSIFICATION_METHODS.keys())[0]
    classification = classify(TEST_SMILES, [method])
    check_classification_with_one_method(TEST_SMILES, classification, method)


def test_classify_only_last_method() -> None:
    method = list(CLASSIFICATION_METHODS.keys())[-1]
    classification = classify(TEST_SMILES, [method])
    check_classification_with_one_method(TEST_SMILES, classification, method)
