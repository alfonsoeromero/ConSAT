#!/usr/bin/env python
"""Binary for the 'coverage' app"""

import sys
from gfam.scripts.coverage import CoverageApp


if __name__ == "__main__":
    sys.exit(CoverageApp().run())
