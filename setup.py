#!/usr/bin/env python3
"""
Setup script for MaMISA
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read version from __init__.py
version = {}
with open("mamisa/__init__.py") as f:
    for line in f:
        if line.startswith("__version__"):
            exec(line, version)

# Read README
readme_file = Path(__file__).parent / "README.md"
long_description = ""
if readme_file.exists():
    long_description = readme_file.read_text()

setup(
    name="mamisa",
    version=version.get("__version__", "1.0.0"),
    description="MaMISA - Manage Misassemblies: Toolkit for metagenomic assembly quality control",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/mamisa",
    license="MIT",
    
    packages=find_packages(),
    
    python_requires=">=3.7",
    
    install_requires=[
        # No external dependencies required for core functionality
    ],
    
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.10',
            'black>=21.0',
            'flake8>=3.8',
        ],
    },
    
    entry_points={
        'console_scripts': [
            'mamisa=mamisa.__main__:main',
        ],
    },
    
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    
    keywords="bioinformatics metagenomics assembly quality-control",
    
    project_urls={
        "Bug Reports": "https://github.com/yourusername/mamisa/issues",
        "Source": "https://github.com/yourusername/mamisa",
    },
)
