# Pymela - Python Matrix Element Analysis
A Python Application for analyzing correlation functions data to extract matrix elements from Lattice QCD measurements.

The package computes ratios of three- and two-point functions and performs fits	to extract matrix elements that lead to Parton Distribution Functions (PDFs).	Currently, plateau/constant fits on the plain ratio and linear fits on the summed ratio	are supported.
From the extracted matrix elements, the	reduced	Ioffe-time distributions are computed for the types of fits supported and they are stored in HDF5 files.

## Contents
* **Pymela**: Contains modules and class definitions related to the operations supported by the	application. Submodules:
	* **io**: Parse and check JSON input files; file conventions.
	* **fit**: Constant and linear fits.
	* **tools**: Tools and utilities, including a module for Jackknife sampling.

* **Tests**: Tests that parse an input JSON file and perform various operations. Currently supported tests and operations are:
	* `tests/read_2pt_corr.py`: Read two-point correlation functions in ASCII format and write the data in HDF5 format.
	* `tests/read_3pt_corr.py`: Read three-point correlation functions in ASCII format and write the data in HDF5 format.
	* `tests/effective_energy.py`: Read two-point correlation functions in ASCII format, compute Effective Energy and perform constant fits on the Effective Energy; write the data in HDF5 format.
	* `tests/compute_ratio.py`: Read two- and three-point correlation functions in ASCII format, compute ratios of three- and two-point functions and store the data in HDF5 format.
	* `tests/fit_ratio.py`: Read two- and three-point correlation functions in ASCII format, compute ratios of three- and two-point functions, perform fits on the ratio to extract matrix elements and store the data in HDF5 format.
	* `tests/compute_rITD.py`: Most comprehensive test. Read two- and three-point correlation functions in ASCII format, compute ratios of three- and two-point functions, perform fits on the ratio to extract matrix elements and compute reduced Ioffe-time distributions (rITD) from the matrix elements. Store all the data in HDF5 format.

## Dependencies
The following packages are required:
* json
* numpy
* h5py
* scipy (scipy.optimize)

## Author and Contact
* **Christos Kallidonis** - Jefferson Lab
* Web: [https://christoskallidonis.com](https://christoskallidonis.com).
* Copyright (C) 2020. All rights reserved.
