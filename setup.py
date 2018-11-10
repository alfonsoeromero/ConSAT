#!/usr/bin/env python
"""ConSAT -- installer script"""

from setuptools import setup, find_packages
from ez_setup import use_setuptools

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"

use_setuptools()

params = {}
params["name"] = "consat"
params["version"] = "1.0"
params["description"] = "The Consensus Signature Architecture Tool"
params["packages"] = find_packages(exclude='tests')
params["scripts"] = ["bin/gfam", "bin/consat", "bin/automated_consat"]

setup(**params)
