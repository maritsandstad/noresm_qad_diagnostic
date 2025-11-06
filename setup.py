from setuptools import setup, find_packages
from pathlib import Path

README = Path(__file__).parent.joinpath("README.md").read_text(encoding="utf-8")

setup(
    name="noresm-qad-diagnostic",
    version="0.1.0",
    description="Minimal packaging for the noresm_qad_diagnostic module",
    long_description=README,
    long_description_content_type="text/markdown",
    license="MIT",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    python_requires=">=3.12",
    install_requires=[
        "cartopy==0.24.0",
        "cftime==1.6.4",
        "matplotlib==3.10.3",
        "numpy==2.2.6",
        "xarray==2025.4.0",
        "xesmf==0.8.7",
        "netCDF4==1.7.2"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering",
    ],
)