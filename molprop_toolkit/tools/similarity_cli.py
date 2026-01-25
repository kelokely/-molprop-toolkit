#!/usr/bin/env python
"""
molprop-similarity CLI entry point.

This module provides the command-line interface for similarity searching.
It wraps the core similarity module functionality.
"""

import sys
import os

# Add tools directory to path for the main implementation
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tools.molprop_similarity import main as _main


def main():
    """Entry point for molprop-similarity command."""
    return _main()


if __name__ == "__main__":
    sys.exit(main())
