#!/usr/bin/env python
"""ConSAT -- installer script"""

__author__  = "Alfonso E. Romero"
__email__   = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"

from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

params = {}
params["name"] = "consat"
params["version"] = "1.0"
params["description"] = "The Consensus Signature Architecture Tool"

params["packages"] = find_packages(exclude='tests')
params["scripts"] = ["bin/gfam", "bin/consat", "bin/automated_consat"]

setup(**params)
