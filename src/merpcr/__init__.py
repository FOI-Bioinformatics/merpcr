"""
merPCR - Modern Electronic Rapid PCR

Python reimplementation of the me-PCR (Multithreaded Electronic PCR) program.
"""

__version__ = "1.0.0"
__author__ = "merPCR Contributors"
__license__ = "GPL-3.0"

from .core.engine import MerPCR
from .core.models import FASTARecord, STSHit, STSRecord

__all__ = ["MerPCR", "STSRecord", "FASTARecord", "STSHit"]
