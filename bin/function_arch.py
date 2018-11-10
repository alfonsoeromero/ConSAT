#!/usr/bin/env python
"""Binary for the 'function_arch' app"""

import sys
from gfam.scripts.function_arch import TransferFunctionFromDomainArch

if __name__ == "__main__":
    sys.exit(TransferFunctionFromDomainArch().run())
