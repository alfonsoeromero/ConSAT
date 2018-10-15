#!/usr/bin/env python
"""Binary for the 'hmm_scan' app"""

import sys
from gfam.scripts.hmm_scan import HMMScanApp

if __name__ == "__main__":
    sys.exit(HMMScanApp().run())
