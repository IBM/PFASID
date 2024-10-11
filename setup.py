#
# Copyright IBM Corp. 2024 - 2024
# SPDX-License-Identifier: MIT
#

import setuptools

MODULE_NAME = "PFAS_Classifier"


with open("requirements.txt") as f:
    requires = f.read().splitlines()


def get_version() -> str:
    with open("./VERSION", "r+") as f:
        return f.readlines()[-1]


def main():
    setuptools.setup(
        name=MODULE_NAME,
        version=get_version(),
        packages=setuptools.find_packages(),
        package_dir={MODULE_NAME: MODULE_NAME},
        package_data={"": ["ECHA_2021_struct_subs.json"]},
        install_requires=requires,
        include_package_data=True,
    )


if __name__ == "__main__":
    print(f"Installing {MODULE_NAME} (version {get_version()})")
    print(setuptools.find_packages())
    main()
