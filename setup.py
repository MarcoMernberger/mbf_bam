#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Use setup.cfg to configure your project.

    This file was generated with PyScaffold 3.1.
    PyScaffold helps you to put up the scaffold of your new Python project.
    Learn more under: https://pyscaffold.org/
"""
import sys

from pkg_resources import require, VersionConflict
from setuptools import setup

from setuptools_rust import Binding, RustExtension

try:
    require("setuptools>=38.3")
except VersionConflict:
    print("Error: version of setuptools is too old (<38.3)!")
    sys.exit(1)


if __name__ == "__main__":
    setup(
        #use_pyscaffold={'git_describe_command': 
                        #r'cat setup.cfg | grep "[metadata]" -A10 | grep "version ?=" | sed "s/version = //"'},
        use_pyscaffold = True,
        rust_extensions=[
            RustExtension("mbf_bam.mbf_bam", binding=Binding.PyO3, debug=False)
        ],

    )
