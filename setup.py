#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Setup script for LineageQuant."""

from setuptools import setup, find_packages

with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="lineagequant",
    version="1.0.0",
    description="Haplotype-resolved expression quantification for polyploid "
                "genomes: lineage assignment (k-mer/SNP) and global EM read deconvolution",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Yi Chen, Gengrui Zhu",
    url="https://github.com/GengruiZhu/LineageQuant",
    license="Non-Commercial Research License",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "lineagequant": ["core/*.sh"],
    },
    python_requires=">=3.8",
    # NOTE: install_requires only covers the PYTHON deps that pip can install.
    # The pipeline also needs external CONDA tools (bowtie2, samtools, kmc,
    # hisat2, minimap2, gffread, mafft, bcftools, ...). Install those via:
    #     conda env create -f environment.yml
    # or install missing ones individually; `lineagequant run` will tell you
    # exactly which tool is missing and how to install it.
    install_requires=[
        "numpy",
        "pandas",
        "pysam",
        "biopython",
    ],
    entry_points={
        "console_scripts": [
            "lineagequant=lineagequant.cli:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Development Status :: 5 - Production/Stable",
    ],
)
