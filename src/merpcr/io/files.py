"""
Opening of plain and compressed input files.

Compression is chosen from the file name extension: .gz (gzip), .bz2 (bzip2)
and .xz or .lzma (xz). Any other name is opened as a plain text file.
"""

import bz2
import gzip
import lzma
import os
from typing import IO, Any

COMPRESSED_SUFFIXES = (".gz", ".bz2", ".xz", ".lzma")


def is_compressed(filename: str) -> bool:
    """Return True if the file name has a recognised compression suffix."""
    return filename.lower().endswith(COMPRESSED_SUFFIXES)


def open_text(filename: str) -> IO[Any]:
    """Open a possibly compressed file for reading as text."""
    name = filename.lower()
    if name.endswith(".gz"):
        return gzip.open(filename, "rt")
    if name.endswith(".bz2"):
        return bz2.open(filename, "rt")
    if name.endswith((".xz", ".lzma")):
        return lzma.open(filename, "rt")
    return open(filename, "r")


def file_size(filename: str) -> int:
    """Return the size of the file on disk in bytes."""
    return os.path.getsize(filename)
