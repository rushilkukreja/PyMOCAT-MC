"""Setup configuration for PyMOCAT-MC package."""

from setuptools import setup, find_packages
import os

# Read the README file for the long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements from requirements.txt
def read_requirements():
    """Read requirements from requirements.txt file."""
    with open("requirements.txt", "r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]

setup(
    name="pymocat-mc",
    version="1.0.0",
    author="Rushil Kukreja, Edward J. Oughton, Giovanni Lavezzi, Enrico M. Zucchelli, Daniel Jang, Richard Linares",
    author_email="rushil.kukreja@tjhsst.edu",
    description="Python implementation of the MIT Orbital Capacity Assessment Toolbox - Monte Carlo",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rushilkukreja/PyMOCAT-MC-2",
    project_urls={
        "Bug Tracker": "https://github.com/rushilkukreja/PyMOCAT-MC-2/issues",
        "Documentation": "https://github.com/rushilkukreja/PyMOCAT-MC-2#readme",
        "Source": "https://github.com/rushilkukreja/PyMOCAT-MC-2",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "python_implementation"},
    packages=find_packages(where="python_implementation"),
    python_requires=">=3.8",
    install_requires=read_requirements() if os.path.exists("requirements.txt") else [
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "pandas>=1.3.0",
        "matplotlib>=3.5.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
            "mypy>=0.950",
        ],
        "docs": [
            "sphinx>=4.5.0",
            "sphinx-rtd-theme>=1.0.0",
            "sphinx-autodoc-typehints>=1.18.0",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["supporting_data/**/*", "examples/**/*"],
    },
    entry_points={
        "console_scripts": [
            "pymocat-mc=mocat_mc:main",
        ],
    },
    keywords="orbital mechanics, space debris, Monte Carlo simulation, space situational awareness, satellite tracking",
    zip_safe=False,
)