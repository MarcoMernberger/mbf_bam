#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import pytest
import sys
from pathlib import Path

# from pypipegraph.testing.fixtures import new_pipegraph, pytest_runtest_makereport  # noqa:F401

root = Path(__file__).parent.parent
sys.path.append(str(root / "src"))
from pypipegraph.testing.fixtures import new_pipegraph  # noqa:F401
