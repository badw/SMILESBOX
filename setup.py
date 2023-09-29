"""
smilesbox: simulation structures from SMILES.
"""

from setuptools import find_packages, setup

from sumo import __version__

with open("README.md") as file:
    long_description = file.read()

setup(
    name="smilesbox",
    version=__version__,
    description="SMILESBOX",
    url="https://github.com/badw/SMILESBOX",
    author="Benjamin A. D. Williamson",
    author_email="benjamin.williamson@ntnu.no",
    long_description=long_description,
    license="MIT",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords="chemistry ase vasp structure",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "ase",
        "tqdm",
        "wget",
        "numpy",
 #       "openbabel"
    ],
    entry_points={
        "console_scripts":[
            "smilesbox-generate = smilesbox.cli:main"
            ]
        }
)
