"""
Setup script for merPCR with Cython extensions.
"""

import numpy as np
from setuptools import Extension, setup

# Optional Cython compilation
try:
    from Cython.Build import cythonize

    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False
    print("Cython not available, building without optimized extensions")

extensions = []

if USE_CYTHON:
    extensions = [
        Extension(
            "merpcr.core.engine_fast",
            sources=["src/merpcr/core/engine_fast.pyx"],
            include_dirs=[np.get_include()],
            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
            extra_compile_args=[
                "-O3",  # Highest optimization
                "-march=native",  # Use CPU-specific optimizations
                "-ffast-math",  # Fast math operations
            ],
        )
    ]

    # Cythonize with compiler directives for maximum performance
    extensions = cythonize(
        extensions,
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,  # Disable bounds checking for speed
            "wraparound": False,  # Disable negative indexing
            "cdivision": True,  # Use C division (faster)
            "initializedcheck": False,  # Don't check if memoryviews are initialized
            "nonecheck": False,  # Don't check if variables are None
            "embedsignature": True,  # Embed function signatures in docstrings
        },
        annotate=True,  # Generate HTML annotation files for optimization analysis
    )

setup(
    ext_modules=extensions,
)
