"""
Entry point for running merPCR as a module: python -m merpcr
"""

from .cli import main

if __name__ == "__main__":
    exit(main())
