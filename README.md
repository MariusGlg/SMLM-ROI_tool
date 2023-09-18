# SMLM-ROI_tool

Script loads a single picasso DBSCAN or dbcluster file, plots the single-molecule localizations and allows the user to draw a region of interes (ROI) on top of the plot. Saves localizations within the ROI if mask is set to False in the config.ini file or saves localizations outside the ROI (Mask == True, imaged is masked). Modified datasets are saved as "name_Mask" or "name_ROI" in the path provided in the config.ini file. Additionally generates and saves a .xlsm file that contains ROI edge coordinates. 

    Loads HDF5 files (Picasso DBSCAN or dbcluster files 
    Plots the data
    Allows user to draw a ROI on top of the plot
    Saves localizations within or outside the ROI depending on user settings.

Requirements: python 3.7, mpl_point_clicker, mpl_interactions os, configparser, h5py, numpy, matplotlib, yaml, numba, pandas

Input file: Picasso[1] hdf5 (picasso dbscan or dbcluster file)

Execution: main.py

Config file:

[INPUT_FILES] path1: path to picasso dbscan/dbcluster files filename: name of picasso dbscan/dbcluster file (name.hdf5) Mask: True (mask is generated), False (ROI is generated)


links: [1] https://github.com/jungmannlab/picasso [2] https://mpl-point-clicker.readthedocs.io/en/latest/
