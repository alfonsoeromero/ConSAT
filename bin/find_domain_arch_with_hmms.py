#!/usr/bin/env python
"""Binary for the 'find_domain_arch_with_hmms' app"""

import sys
from gfam.scripts.find_domain_arch_with_hmms\
    import FindDomainArchitectureWithHMMsApp

if __name__ == "__main__":
    sys.exit(FindDomainArchitectureWithHMMsApp().run())
