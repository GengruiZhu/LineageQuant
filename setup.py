from setuptools import setup, find_packages
import os

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="lineagequant",
    version="0.1.0",
    author="Yi Chen, Gengrui Zhu",
    description=(
        "Allele-resolved RNA-seq quantification for polyploid "
        "and hybrid genomes"
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    packages=find_packages(),
    package_data={
        "lineagequant": ["scripts/*.sh"],
    },
    install_requires=[
        "biopython",
        "pysam",
        "pandas",
    ],
    entry_points={
        "console_scripts": [
            "lineagequant=lineagequant.pipeline:main",
        ],
    },
    python_requires=">=3.7",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
)
