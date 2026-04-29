"""
MetaboGraph Setup Configuration
"""

from setuptools import setup, find_packages
import os

# Read the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="metabograph",
    version="1.0.1",
    author="Oluwatosin Daramola",
    author_email="oluwatosin.daramola@ttu.edu",
    description="Metabolite Annotation & Pathway Network Analysis Tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/odaramola92/metabograph",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "metabograph=gui.main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.feather", "*.txt"],
    },
)
