#!/usr/bin/env python
"""Binary for the 'extract_gene_ids' app"""

import sys
from gfam.scripts.extract_gene_ids import ExtractGeneIDsApp

if __name__ == "__main__":
    sys.exit(ExtractGeneIDsApp().run())
