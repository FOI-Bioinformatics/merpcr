#!/usr/bin/env python3
"""
Setup script for merPCR
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="merpcr",
    version="1.0.0",
    author="Andreas Sjödin",
    author_email="andreas.sjodin@gmail.com",
    description="Modern Electronic Rapid PCR - Python reimplementation of me-PCR",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/druvus/merpcr",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.12",
    install_requires=[
        "typing_extensions>=4.0.0",
    ],
    entry_points={
        "console_scripts": [
            "merpcr=merPCR:main",
        ],
    },
)