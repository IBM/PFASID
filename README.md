# PFAS_Classifier

Package for classifying molecules (based on SMILES) under different PFAS regulations.

## Installation

To use this library in you python project you might install it with pip in a fresh environment:

```cli
python3 -m venv pfas_classifier
source pfas_classifier/bin/activate
python3 -m pip install git+ssh://git@github.ibm.com/pfacts-nsf/PFAS_Classifier.git

```

OR

```cli
pipenv install -e git+ssh://git@github.ibm.com/pfacts-nsf/PFAS_Classifier.git@main#egg=PFAS_Classifier
pipenv shell
```

Check if the package was installed with `pip list`.

## Use

```cli
>>> import PFAS_Classifier as classifier
>>> classifier.classify(['CC(C(=O)O)c1ccc(C(F)(F)N(C)C)cc1'],['EPA_2024_0'])
>>> {'EPA_2024_0': {'CC(C(=O)O)c1ccc(C(F)(F)N(C)C)cc1': False}}

>>> classifier.get_classifiers()
>>> {'EPA_2024_0': {'description': '', 'reference': 'https://pubs.acs.org/doi/10.1021/cen-10040-polcon1'}, 'EPA_2023_list': {'description': '', 'reference': 'https://www.epa.gov/ccl/ccl-5-chemical-contaminants'}, 'ECHA_2021_struct': {'description': '', 'reference': 'https://echa.europa.eu/de/registry-of-restriction-intentions/-/dislist/details/0b0236e18663449b'}}
```

## Testing

1. Install test dependencies at [requirements-test.txt](requirements-test.txt)

```cli
pip install -r requirements-test.txt
```

2. Run automated tests with PyTest

```cli
pytest tests/test.py
```
