#!/usr/bin/env python
"""Binary for the 'overrep' app"""

import sys
from gfam.scripts.overrep import OverrepresentationAnalysisApp

if __name__ == "__main__":
    sys.exit(OverrepresentationAnalysisApp().run())
