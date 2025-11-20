import copy

import numpy as np

import xarray as xr

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

import os, sys, glob

# import netCDF4 as nc4
import warnings

warnings.filterwarnings("ignore")

from .infrastructure_help_functions import setup_nested_folder_structure_from_dict, read_pam_file
from .noresm_full_model import NorESMFullModel
from .noresm_qad_diagnostic import MONTHS, SEASONS

class NorESMQADPreprocessor:
    """
    Class to preprocess NorESM output from one or more NorESM output runs
    Setting up folders or checking for previously calculated files,
    keeping track of configs for various models etc.
    """
    
    def __init__(self, runs_dict, config_dict, paths_preproc, path_plots = None):
        self.runs_dict = runs_dict
        self.config_dict = config_dict
        self.paths_preproc = paths_preproc
        if path_plots is None:
            self.paths_plot = os.path.join(self.paths_preproc, "mask_plots")
        self.done_list = None
        self.check_if_preproc_files_exits()
        self.var_to_comp_map = {}
        self.make_variable_to_comp_inverse_mapping()
        print(self.done_list)
        self.setup_and_preprocess_missing()
        self.make_mask_plots()

    def check_if_preproc_files_exits(self):
        if self.done_list is None:
            self.done_list = {}
        for run in self.runs_dict.keys():
            if run not in self.done_list:
                self.done_list[run] = []
            elif "All done" in self.done_list[run]:
                continue
            all_done = True
            for component, var_list in self.config_dict["VAR_LIST_MAIN"].items():
                print(component)
                if component in self.done_list[run]:
                    continue
                comp_done = True
                for var in var_list:
                    if var in self.done_list[run]:
                        continue
                    if os.path.exists(f"{self.paths_preproc}/{run}/{var}_{component}_{run}_preproc.nc"):
                        self.done_list[run].append(var)
                    elif var == "AMOC" and os.path.exists(f"{self.paths_preproc}/{run}/amoc_26N_{component}_{run}_preproc.nc"):
                        self.done_list[run].append(var)
                    else:
                        comp_done = False
                        all_done = False
                if comp_done:
                    self.done_list[run].append(component)
            if all_done:
                self.done_list[run].append("All done")

    def setup_and_preprocess_missing(self):
        for run, path in self.runs_dict.items():
            if "All done" in self.done_list[run]:
                continue
            if not os.path.exists(f"{self.paths_preproc}/{run}/"):
                os.mkdir(f"{self.paths_preproc}/{run}/")
            pam_dict_copy = self.make_pams_to_send_to_noresm_object(self.done_list[run])
            noresmfull = NorESMFullModel(path, pamfile=pam_dict_copy, casename=run)
            for component in pam_dict_copy["VAR_LIST_MAIN"].keys():
                noresmfull.components[component].get_area_mean_ts_data_per_variable(self.paths_preproc, varlist_in = pam_dict_copy["VAR_LIST_MAIN"][component])
                noresmfull.components[component].delete_regridder()


    def make_variable_to_comp_inverse_mapping(self):
        for comp, varlist in self.config_dict["VAR_LIST_MAIN"].items():
            for variable in varlist:
                self.var_to_comp_map[variable] = comp
    
    def make_pams_to_send_to_noresm_object(self, run_done_dict):
        pam_dict_copy = copy.deepcopy(self.config_dict)
        for comp in self.config_dict["VAR_LIST_MAIN"].keys():
            if comp in run_done_dict:
                pam_dict_copy["VAR_LIST_MAIN"].pop(comp)
                continue
            to_pop = []
            for variable in self.config_dict["VAR_LIST_MAIN"][comp]:
                if variable in run_done_dict:
                    to_pop.append(variable)
            for variable in to_pop:
                pam_dict_copy["VAR_LIST_MAIN"][comp].remove(variable)
        return pam_dict_copy


    def make_mask_plots(self):
        pass
