import numpy as np

import xarray as xr

import matplotlib.pyplot as plt
import matplotlib as mpl

import os, sys, glob

# import netCDF4 as nc4
import warnings

warnings.filterwarnings("ignore")

import cartopy.crs as ccrs

from .noresm_full_model import NorESMFullModel
from .plotting_methods import make_bias_plot, make_bias_plot_latixy_longxy
from .infrastructure_help_functions import setup_nested_folder_structure_from_dict, read_pam_file, make_monthly_taxis_from_years#, clean_empty_folders_in_tree
from  .misc_help_functions import get_unit_conversion_and_new_label, make_regridding_target_from_weightfile, get_unit_conversion_from_string, do_light_unit_string_conversion

MONTHS = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]

SEASONS = ["DJF", "MAM", "JJA", "SON"]

DRIFT_YRS = 30
#mpl.rc('font', size=30) # Set default font size to 20

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
        self.var_to_comp_map = {}
        self.make_variable_to_comp_inverse_mapping()

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
                "clim_maps": None, 
                "perf_metrics": None, 
            }
        }
                
        setup_nested_folder_structure_from_dict(outdir, subfolder_structure)
        self.outdir = f"{outdir}/{self.casename}"

    def make_variable_to_comp_inverse_mapping(self):
        for comp, varlist in self.varpams["VAR_LIST_MAIN"].items():
            for variable in varlist:
                self.var_to_comp_map[variable] = comp
        to_pop = []
        for variable in self.varpams['OBS_COMPARE_MAPS'].keys():
            if variable not in self.var_to_comp_map.keys():
                to_pop.append(variable)
        for variable in to_pop:
            self.varpams['OBS_COMPARE_MAPS'].pop(variable)


    def make_all_map_plots(self, title = None):
        if title is None:
            title = f"Last 20 year climatology comparison {self.casename}"
        print(self.varpams["OBS_COMPARE_MAPS"])
        fig = plt.figure(constrained_layout=True, figsize=(30,40))
        fig.suptitle(title, fontsize=20)
        subfigs = fig.subfigures(nrows=len(self.varpams['OBS_COMPARE_MAPS'].keys()), ncols=1)    
        for subfignum, variable in enumerate(self.varpams['OBS_COMPARE_MAPS'].keys()):
            if len(self.varpams['OBS_COMPARE_MAPS'].keys()) == 1:
                subfignow = subfigs
            else:
                subfignow = subfigs[subfignum] 
            self.make_single_comparison_map(variable, subfignow)
        fig.savefig(f"{self.outdir}/clim_maps/climatology_maps_{self.casename}.png")

    def make_comparison_map_no_comparison_data(self, variable, subfig_handle):
        #Put in single map-plot
        print(f"No observational data for {variable}")
        axs = subfig_handle.subplots(
            nrows=1, 
            ncols=1,            
            subplot_kw={"projection": ccrs.Robinson()},
            )
        plot_data = self.main_run.components[self.var_to_comp_map[variable]].get_annual_mean_map_data(varlist_in=variable)
        print(plot_data)
        if variable.startswith("aice"):
            make_bias_plot_latixy_longxy(plot_data, plot_data["TLAT"], plot_data["TLON"], self.casename, ax=axs)
        make_bias_plot(plot_data, self.casename, ax = axs)

    def make_single_comparison_map(self, variable, subfig_handle):
        print(variable)
        subfig_handle.suptitle(f"20 year climatological comparison for {variable}", fontsize=20)
        if len(self.varpams['OBS_COMPARE_MAPS'][variable]) == 1 and self.varpams['OBS_COMPARE_MAPS'][variable][0] == 'False':
            self.make_comparison_map_no_comparison_data(variable, subfig_handle)
        elif self.varpams['OBS_COMPARE_MAPS'][variable][0] not in self.ilamb_confs.configurations[variable].obsdatasets.keys():
            self.make_comparison_map_no_comparison_data(variable, subfig_handle)
        else:
            # Put in map + map + changeplot
            print(self.varpams['OBS_COMPARE_MAPS'][variable])
            axs = subfig_handle.subplots(
                nrows=1, 
                ncols = 3, 
                subplot_kw={"projection": ccrs.Robinson()},
                )
            yminv, ymaxv, diffrange, negdiffrange = self.ilamb_confs.configurations[variable].obs_limits
            # Make this work...
            #unit_conversion_factor, unit_to_print = get_unit_conversion_from_string(self.ilamb_confs.get_variable_plot_unit(variable), self.unit_dict[variable])
            varname_ilamb =  self.ilamb_confs.get_vaname_in_ilamb_cfgs(variable)
            plot_data = self.main_run.components[self.var_to_comp_map[variable]].get_annual_mean_map_data(varlist_in=variable)
            make_bias_plot(
                plot_data,
                f"{self.casename}",
                ax=axs[0],
                yminv = yminv,
                ymaxv = ymaxv,
                xlabel = f"{variable}" #[{unit_to_print}]"
                )
            plot_data_obs = self.ilamb_confs.get_data_for_map_plot(varname_ilamb, self.varpams['OBS_COMPARE_MAPS'][variable][0], self.main_run.components[self.var_to_comp_map[variable]].regrid_target)
            make_bias_plot(
                plot_data_obs,
                f"{self.varpams['OBS_COMPARE_MAPS'][variable][0]}",
                ax=axs[1],
                yminv = yminv,
                ymaxv = ymaxv,
                xlabel = f"{variable}"#[{unit_to_print}]"
                )
            make_bias_plot(
                plot_data_obs,
                f"{self.casename} - {self.varpams['OBS_COMPARE_MAPS'][variable][0]}",
                ax=axs[2],
                yminv = negdiffrange,
                ymaxv = diffrange,
                xlabel = f"{variable}",#[{unit_to_print}]"
                cmap = "RdYlBu_r"
                )
            axs[0].set_title(self.casename, fontsize=20)
            axs[1].set_title(self.varpams['OBS_COMPARE_MAPS'][variable][0], fontsize=20)
            axs[2].set_title("bias", fontsize=20)

    def make_all_timeseries_plots(self):
        self.setup_ts_figs_and_axs()

    def setup_ts_figs_and_axs(self, title = None):
        if title is None:
            title = f"Averaged trends {self.casename}"
        fig = plt.figure(constrained_layout=True, figsize=(40,40))
        fig.suptitle(title, fontsize=20)
        subfigs = fig.subfigures(nrows=len(self.varpams['VAR_LIST_MAIN'].keys()), ncols=1)
        print(self.varpams['VAR_LIST_MAIN'])
        #sys.exit(4)
        for subfignum, (name, items) in enumerate(self.varpams['VAR_LIST_MAIN'].items()):
            if len(self.varpams['VAR_LIST_MAIN'].keys()) == 1:
                subfignow = subfigs
            else:
                subfignow = subfigs[subfignum]
            subfignow.suptitle(f"Trends from {name}", fontsize=20)
            rownums = len(items)//5 + 1
            colnums = np.min((5, len(items)))
            axs = subfignow.subplots(nrows=rownums, ncols=colnums)
            tyaxis = self.main_run.components[name].get_year_range()
            outd_year = self.main_run.components[name].get_area_mean_ts_data()
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
                if item == "AMOC":
                    for data_var in outd_year.data_vars:
                        # Could potentially bring in additional amoc lines here
                        if data_var.startswith("amoc_26N"):
                            axnow.plot(tyaxis, outd_year[data_var].values)
                        outd_ts = outd_year["amoc_26N"].values
                    #axnow.legend(fontsize = 10)
                else:
                    if len(outd_year[item].shape) > 0:
                        outd_ts = outd_year[item].values.flatten()
                    else:
                        outd_ts = outd_year[item].values
                    # Add multiple lines for amoc..., how to...
                    axnow.plot(tyaxis, outd_ts)
                axnow.set_title(f"{item} ({self.main_run.components[name].get_variable_unit(item)})", fontsize=20)
                axnow.set_xlabel("Year", fontsize=20)
                axnow.set_ylabel(f"{item} ({self.main_run.components[name].get_variable_unit(item)})", fontsize=20)
                print(len(outd_ts))
                print(len(tyaxis))
                print(outd_ts)
                self.ilamb_confs.add_target_and_drift(item, axnow, outd_ts[np.max((0, len(tyaxis)-DRIFT_YRS)):], tyaxis)


        fig.savefig(f"{self.outdir}/trends_drift/trend_overview_{self.casename}_yearly.png")

    def add_to_ts_ax(self, vari_info):
        pass

    def check_drift_criteria(self):
        pass

    def make_perf_heatmaps(self):
        pass