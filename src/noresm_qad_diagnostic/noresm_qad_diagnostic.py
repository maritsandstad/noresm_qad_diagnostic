import numpy as np

import xarray as xr

import matplotlib.pyplot as plt

import os, sys, glob

# import netCDF4 as nc4
import warnings

warnings.filterwarnings("ignore")

import cartopy.crs as ccrs

from .noresm_full_model import NorESMFullModel
from .plotting_methods import make_generic_regridder, regrid_se_data, make_bias_plot
from .infrastructure_help_functions import setup_nested_folder_structure_from_dict, read_pam_file, make_monthly_taxis_from_years#, clean_empty_folders_in_tree
from  .misc_help_functions import get_unit_conversion_and_new_label, make_regridding_target_from_weightfile, get_unit_conversion_from_string, do_light_unit_string_conversion

MONTHS = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]

SEASONS = ["DJF", "MAM", "JJA", "SON"]


class NorESMQADD:
    """
    Class to hold a general NorESM diagnostic object
    with functionality to setup output folders, keep track of
    runs, and make plots
    """

    def __init__(self, main_run, pamfile, comp_runs=None, ilamb_confs=None, outdir=None, casename=None):
        # Probably setup this:
        self.main_run = NorESMFullModel(main_run, pamfile=pamfile)
        if casename is None:
            self.casename = self.main_run.casename
        else:
            self.casename = casename
        # Probably loop over setting up these:
        self.comp_run_list = comp_runs
        self.ilamb_confs = ilamb_confs
        # Update method to do this: 
        self.varpams = read_pam_file(pamfile)
        self.cut_pamfile_to_existing()
        if outdir is None:
            outdir = "figs/"
        self.setup_folder_structure(outdir)

    def cut_pamfile_to_existing(self):
        print(self.varpams)
        for comp, variables in self.varpams["VAR_LIST_MAIN"].items():
            vars_missing = []
            for item in variables:
                if item not in self.main_run.components[comp].varpams:
                    vars_missing.append(item)
            self.varpams["VAR_LIST_MAIN"][comp] = list(set(self.varpams["VAR_LIST_MAIN"][comp]) - set(vars_missing))
        # TODO: Deal with OBS_COMPARE


    def setup_folder_structure(self, outdir):
        if not os.path.exists(outdir):
            raise ValueError(f"{outdir} must be an existing directory")
        
        subfolder_structure = {
            f"{self.main_run.casename}": {
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
        self.setup_ts_figs_and_axs()

    def setup_ts_figs_and_axs(self, title = None):
        print(self.varpams['VAR_LIST_MAIN'])
        if title is None:
            title = f"Averaged trends {self.casename}"
        fig = plt.figure(constrained_layout=True, figsize=(30,30))
        fig.suptitle(title)
        subfigs = fig.subfigures(nrows=len(self.varpams['VAR_LIST_MAIN'].keys()), ncols=1)
        print(len(self.varpams.keys()))
        print(self.varpams.keys())

        for subfignum, (name, items) in enumerate(self.varpams['VAR_LIST_MAIN'].items()):
            print(subfignum)
            if len(self.varpams['VAR_LIST_MAIN'].keys()) == 1:
                subfignow = subfigs
            else:
                subfignow = subfigs[subfignum]
            subfignow.suptitle(f"Trends from {name}")
            print(items)
            print(f"{name}")
            print(len(items))
            rownums = len(items)//5 + 1
            colnums = np.min((5, len(items)))
            axs = subfignow.subplots(nrows=rownums, ncols=colnums)
            tyaxis = self.main_run.components[name].get_year_range()
            taxis = make_monthly_taxis_from_years(tyaxis)
            outd_year, outd = self.main_run.components[name].get_area_mean_ts_data(varlist=items)
            for subnum, item in enumerate(items):
                if colnums == 1:
                    axnow = axs
                elif rownums == 1:
                    axnow = axs[subnum]
                    print(f"{item}, {subnum}")
                else:
                    axnow = axs[subnum//5, subnum%5]
                    print(f"{item}, {subnum}, {subnum//5}, {subnum%5}")
                # Monthly:
                #axnow.plot(taxis, outd[item].values.flatten())
                # Yearly
                if len(outd_year[item].shape) > 0:
                    axnow.plot(tyaxis, outd_year[item].values.flatten())
                else:
                    axnow.plot(tyaxis, outd_year[item])
                axnow.set_title(f"{item} ({self.main_run.components[name].unit_dict[item]})")
                axnow.set_xlabel("Year")
                axnow.set_ylabel(f"{item} ({self.main_run.components[name].unit_dict[item]})")

        fig.savefig(f"{self.outdir}/trends_drift/trend_overview_{self.casename}_yearly.png")

    def add_to_ts_ax(self, vari_info):
        pass

    def check_drift_criteria(self):
        pass

    def make_perf_heatmaps(self):
        pass