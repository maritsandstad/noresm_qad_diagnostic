import numpy as np

import xarray as xr

import matplotlib.pyplot as plt

import os, sys, glob

# import netCDF4 as nc4
import warnings

warnings.filterwarnings("ignore")

import cartopy.crs as ccrs

from .plotting_methods import make_generic_regridder, regrid_se_data, make_bias_plot
from .infrastructure_help_functions import setup_nested_folder_structure_from_dict, read_pam_file#, clean_empty_folders_in_tree
from  .misc_help_functions import get_unit_conversion_and_new_label, make_regridding_target_from_weightfile, get_unit_conversion_from_string, do_light_unit_string_conversion

MONTHS = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]

SEASONS = ["DJF", "MAM", "JJA", "SON"]


class NorESMQADD:
    """
    Class to hold a general NorESM diagnostic object
    with functionality to setup output folders, keep track of
    runs, and make plots
    """

    def __init__(self, main_run, comp_runs=None, ilamb_confs=None, outdir=None, pamfile=None):
        # Probably setup this:
        self.main_run = main_run
        # Probably loop over setting up these:
        self.comp_run_list = comp_runs
        self.ilamb_confs = ilamb_confs
        # Update method to do this: 
        self.var_pams = read_pam_file(pamfile)
        if outdir is None:
            outdir = "figs/"
        self.setup_folder_structure(outdir)

    def setup_folder_structure(self, outdir):
        if not os.path.exists(outdir):
            raise ValueError(f"{outdir} must be an existing directory")
        
        subfolder_structure = {
            f"{self.casename}": {
                "trends_drift": None, 
                "clim_maps": ["ANN", "DJF", "MAM", "JJA", "SON"], 
                "perf_metrics": None, 
            }
        }
                
        setup_nested_folder_structure_from_dict(outdir, subfolder_structure)
        self.outdir = f"{outdir}/{self.casename}"

    def make_all_map_plots(self):
        pass

    def make_single_comparison_map(self, vari_info):
        pass

    def make_all_timeseries_plots(self):
        pass

    def add_to_ts_ax(self, vari_info):
        pass

    def check_drift_criteria(self):
        pass

    def make_perf_heatmaps(self):
        pass