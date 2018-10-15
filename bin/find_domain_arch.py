#!/usr/bin/env python
"""Binary for the 'find_domain_arch' app"""

import sys
from gfam.scripts.find_domain_arch import FindDomainArchitectureApp

if __name__ == "__main__":
    sys.exit(FindDomainArchitectureApp().run())
