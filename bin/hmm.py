#!/usr/bin/env python
"""Binary for the 'hmm' app"""

import sys
from gfam.scripts.hmm import HMM

if __name__ == "__main__":
    sys.exit(HMM().run())
