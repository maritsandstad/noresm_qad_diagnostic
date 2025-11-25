import copy
import grp
import stat
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
from .noresm_full_model import NorESMFullModel, KEY_COMP_MAPPING
from .noresm_qad_diagnostic import MONTHS, SEASONS
from .noresm_concrete_components import NorESMAtmComponent, NorESMLndComponent, NorESMOcnComponent, NorESMOcnbgcComponent, NorESMGlcComponent, NorESMIceComponent

class NorESMQADPreprocessor:
    """
    Class to preprocess NorESM output from one or more NorESM output runs
    Setting up folders or checking for previously calculated files,
    keeping track of configs for various models etc.
    """
    
    def __init__(self, runs_dict, config_dict, paths_preproc, run_gm_preproc=True, path_plots = None):
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
        if run_gm_preproc:
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

    def check_if_preproc_files_spec_comp_exist(self, comp):
        if self.done_list is None:
            self.done_list = {}
        if comp not in self.done_list:
            self.done_list[comp] = {}
        for run in self.runs_dict.keys():
            if run not in self.done_list[comp]:
                self.done_list[comp][run] = []
            elif "All done" in self.done_list[comp][run]:
                continue
            all_done = True
            if not os.path.exists(f"{self.paths_preproc}/{run}/{comp}_{run}_preproc_climmaps.nc"):
                for var, obslist in self.config_dict["OBS_COMPARE_MAPS"][comp].items():
                    print(var)
                    if var in self.done_list[comp][run]:
                        continue
                    if os.path.exists(f"{self.paths_preproc}/{run}/{var}_{comp}_{run}_preproc_climmap.nc"):
                            self.done_list[comp][run].append(var)
                    else:
                        all_done = False
            if all_done:
                self.done_list[comp][run].append("All done")

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
    
    def make_pams_to_send_to_noresm_object(self, run_done_dict, comp = None):
        if comp is None:
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

        else:
            varlist = []
            for obs_comparison_var in self.config_dict["OBS_COMPARE_MAPS"][comp].keys():
                if "var" in run_done_dict:
                    continue
                varlist.append(obs_comparison_var)
            return varlist
        return pam_dict_copy


    def make_mask_plots(self):
        pass

    def coax_and_regrid_map_plot_component(self, comp="ocn", year_range=None):
        self.check_if_preproc_files_spec_comp_exist(comp)
        for run, path in self.runs_dict.items():
            if "All done" in self.done_list[comp][run]:
                continue
            varlist = self.make_pams_to_send_to_noresm_object(self.done_list[comp][run], comp=comp)
            noresm_comp = KEY_COMP_MAPPING[comp](path, varlist, casename = run)
            outd = noresm_comp.get_annual_mean_map_data(year_range=year_range)

            # If dataset uses y/x dims but lon/lat coordinates exist separately
            # create 2D lon/lat coordinates from available lon/lat variables
            try:
                dims = list(outd.dims)
            except Exception:
                dims = []
            if "y" in dims and "x" in dims:
                # Try to find lon/lat variables in dataset
                lon_var = None
                lat_var = None
                for name in ["lon", "longitude", "LONGXY", "LONGX"]:
                    if name in outd.variables:
                        lon_var = outd[name]
                        lon_name = name
                        break
                for name in ["lat", "latitude", "LATIXY", "LATIXY"]:
                    if name in outd.variables:
                        lat_var = outd[name]
                        lat_name = name
                        break
                # If lon/lat variables found and are 1D or 2D with correct axes,
                # construct 2D coordinate arrays matching y,x
                if lon_var is not None and lat_var is not None:
                    try:
                        lon_arr = lon_var.values
                        lat_arr = lat_var.values
                        # Derive 1D vectors for meshgrid: lon along x, lat along y
                        if lon_arr.ndim == 1:
                            lon_1d = lon_arr
                        elif lon_arr.ndim == 2:
                            lon_1d = lon_arr[0, :]
                        else:
                            lon_1d = None
                        if lat_arr.ndim == 1:
                            lat_1d = lat_arr
                        elif lat_arr.ndim == 2:
                            lat_1d = lat_arr[:, 0]
                        else:
                            lat_1d = None
                        if lon_1d is not None and lat_1d is not None:
                            outd = outd.assign_coords(lon=(("x"), lon_1d), lat=(("y"), lat_1d))
                            # If lon/lat exist as data variables, drop them to avoid duplication
                            data_vars_to_drop = []
                            if lon_name in outd.data_vars:
                                data_vars_to_drop.append(lon_name)
                            if lat_name in outd.data_vars and lat_name != lon_name:
                                data_vars_to_drop.append(lat_name)
                            if len(data_vars_to_drop):
                                outd = outd.drop_vars(data_vars_to_drop)
                            outd = outd.rename({"x": "lon", "y": "lat"})
                    except Exception:
                        pass
            fpath = f"{self.paths_preproc}/{run}/{comp}_{run}_preproc_climmaps.nc"
            outd.to_netcdf(fpath)
            os.chown(fpath, os.getuid(), grp.getgrnam("ns9560k").gr_gid)
            os.chmod(fpath, stat.S_IRWXU | stat.S_IRWXG)


