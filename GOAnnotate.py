#!/usr/bin/env python3
"""
GOAnnotate: functional annotation of gene models using Diamond searches against UniProt.

This script is a thin wrapper around the goannotate package so that
./GOAnnotate.py keeps working from a checkout. Run ./GOAnnotate.py --help
or, after `pip install -e .`, `goannotate --help`.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from goannotate.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
