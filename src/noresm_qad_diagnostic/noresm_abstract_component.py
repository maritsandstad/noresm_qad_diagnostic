import numpy as np

import xarray as xr

import os, glob

# import netCDF4 as nc4
import warnings

from abc import ABC, abstractmethod

warnings.filterwarnings("ignore")

from .plotting_methods import make_generic_regridder, regrid_se_data, make_bias_plot
from .infrastructure_help_functions import setup_nested_folder_structure_from_dict, read_pam_file#, clean_empty_folders_in_tree
from  .misc_help_functions import get_unit_conversion_and_new_label, make_regridding_target_from_weightfile, get_unit_conversion_from_string, do_light_unit_string_conversion



#def get_minimal_intersecting_year_range(year_range, year_range_other):



class NorESMAbstractComponent(ABC):
    """
    Class that holds and organises info on modelling outputs from some component,
    Methods on getting data on various formats etc.
    """
    def __init__(
        self, datapath, varlists, casename=None,
    ):
        self.datapath = datapath
        #self.var_pams = read_pam_file(pamfile)
        self.filelist = self.get_filelist()
        self.filelist.sort()
        self.varpams = varlists
        if casename is None:
            self.casename = ".".join(self.filelist[0].split("/")[-1].split(".")[:-4])
        else:
            self.casename = casename
        self.unit_dict = {}
        vars_missing = self.wet_and_cut_varlists()
        if len(vars_missing):
            print("Not all requested variables are available in output, ignoring these:")
            print(vars_missing)

    def wet_and_cut_varlists(self):
        # TODO: Make sure all items in SEASONAL is also in var_list_main 
        read = xr.open_dataset(self.filelist[0]).keys()
        lists_check =["VAR_LIST_MAIN", "COMPARE_VARIABLES"]

        vars_missing = []
        for list_n in lists_check:
            for item in self.var_pams[list_n]:
                if item not in read:
                    vars_missing.append(item)
            self.var_pams[list_n] = list(set(self.var_pams[list_n]) - set(vars_missing))
        return vars_missing

    def setup_folder_structure(self, outdir):
        if not os.path.exists(outdir):
            raise ValueError(f"{outdir} must be an existing directory")
        
        subfolder_structure = {
            f"{self.casename}": {
                "trends": None, 
                "clim_maps": ["ANN", "DJF", "MAM", "JJA", "SON"], 
                "seasonal_cycle": None, 
            }
        }
        
        setup_nested_folder_structure_from_dict(outdir, subfolder_structure)
        self.outdir = f"{outdir}/{self.casename}"

    #def clean_out_empty_folders(self):
    #    clean_empty_folders_in_tree(self.outdir)
    @abstractmethod
    def get_filelist(self):
        """
        Get a list of all the files for the case for specific component

        Returns
        -------
        list
            with paths to all the files
        """
        pass

    def get_annual_map_data(self, year_range, varlist=None):
        """
        Get annual mean data for variables in varlist

        Parameters
        ----------
        year_range : range
            Of years to include
        varlist : list
            List of variables to get data for, if not supplied, the objects
            varlist will be used. If the list includes variables not in the 
            outputfiles, they will be 
        """
        outd = None
        if varlist is None:
            varlist = self.var_pams["VAR_LIST_MAIN"]
        for year in year_range:
            for month in range(12):
                mfile = f"{self.datapath}/{self.casename}.clm2.h0.{year:04d}-{month + 1:02d}.nc"
                outd_here = xr.open_dataset(mfile, engine="netcdf4")[varlist]
                # print(outd_here)
                # sys.exit(4)
                if not outd:
                    outd = outd_here
                else:
                    outd = xr.concat([outd, outd_here], dim="time")
        outd = outd.mean(dim="time")
        return outd
    
    def get_annual_mean_ts(self, year_range, varlist=None):
        """
        Get annual mean data for variables in varlist

        Parameters
        ----------
        year_range : range
            Of years to include
        varlist : list
            List of variables to get data for, if not supplied, the objects
            varlist will be used. If the list includes variables not in the 
            outputfiles, they will be 
        """
        outd = None
        if varlist is None:
            varlist = self.var_pams["VAR_LIST_MAIN"]
        for year in year_range:
            outd_yr = None
            for month in range(12):                         
                mfile = f"{self.datapath}/{self.casename}.clm2.h0.{year:04d}-{month + 1:02d}.nc"
                outd_here = xr.open_dataset(mfile, engine="netcdf4")[varlist]
                # print(outd_here)
                # sys.exit(4)
                if not outd_yr:
                    outd_yr = outd_here
                else:
                    outd_yr = xr.concat([outd_yr, outd_here], dim="time")
            outd_yr = outd_yr.mean(dim="time")
            if not outd:
                outd = outd_yr
            else: 
                outd = xr.concat([outd, outd_yr], dim="time")
        return outd
    
    def get_area_mean_ts_data(self, varlist = None, year_range=None, area_def = None):
        outd = None
        if varlist is None:
            varlist = self.var_pams["VAR_LIST_MAIN"]
        if year_range is None:
            year_range = self.get_year_range()
        for year in year_range:
            for month in range(12): 
                mfile = f"{self.datapath}/{self.casename}.clm2.h0.{year:04d}-{month + 1:02d}.nc"
                outd_here = xr.open_dataset(mfile, engine="netcdf4")[varlist]
                outd_here = self.regrid(outd_here)
                if area_def is not None:
                    outd_here = outd_regr.sel(
                    lat=slice(area_def["lat_s"], area_def["lat_n"]),
                    lon=slice(area_def["lon_w"], area_def["lon_e"]),
                )                      
                weights = np.cos(np.deg2rad(outd_here.lat))
                weighted_data = outd_here.weighted(weights)
                ts_data = weighted_data.mean(["lon", "lat"])
                if not outd:
                    outd = ts_data
                else:
                    outd = xr.concat([outd, ts_data], dim="time")
        return outd
    
    
    def regrid(self, data):
        return data

    def get_seasonal_data(self, season, year_range, varlist=None):
        """
        Get climatological mean data for a season

        Parameters
        ----------
        season : int
            The season-number 0 is DJF, 1 is MAM, 2 is JJA
            and 3 is SON
        year_range : np.ndarray
            Range of years to include in climatological mean
        varlist : list
            List of variables to make plots of

        Returns
        -------
        xr.Dataset
            Of climatological seasonal means of the variables in varlist
            for the season requested
        """
        outd = None
        if varlist is None:
            varlist = self.var_pams["VAR_LIST_MAIN"]
        for year in year_range:
            for monthincr in range(3):

                month = monthincr + season * 3
                if month == 0:
                    month = 12
                # print(f"Season: {season}, monthincr: {monthincr}, month: {monthincr}")
                mfile = (
                    f"{self.datapath}/{self.casename}.clm2.h0.{year:04d}-{month:02d}.nc"
                )
                outd_here = xr.open_dataset(mfile, engine="netcdf4")[self.var_pams["VAR_LIST_MAIN"]]
                # print(outd_here)
                # sys.exit(4)
                if not outd:
                    outd = outd_here
                else:
                    outd = xr.concat([outd, outd_here], dim="time")
        outd = outd.mean(dim="time")
        return outd

    def get_monthly_climatology_data(self, year_range, varlist=None):
        outd_months = None
        if varlist is None:
            varlist = self.var_pams["VAR_LIST_MAIN"]
        
        for month in range(12):
            outd = None
            for year in year_range:
                # print(f"Season: {season}, monthincr: {monthincr}, month: {monthincr}")
                mfile = f"{self.datapath}/{self.casename}.clm2.h0.{year:04d}-{month+1:02d}.nc"
                outd_here = xr.open_dataset(mfile, engine="netcdf4")[self.var_pams["VAR_LIST_MAIN"]]
                # print(outd_here)
                # sys.exit(4)
                if not outd:
                    outd = outd_here
                else:
                    outd = xr.concat([outd, outd_here], dim="time")
            outd = outd.mean(dim="time", keepdims=True)
            if outd_months is None:
                outd_months = outd
            else:
                outd_months = xr.concat([outd_months, outd], dim="time")
        return outd_months

    def get_year_range(self, short=False):
        year_start, year_end, files_missing = self.find_case_year_range()
        # TODO: Deal with missing files, also for different year_range
        if not files_missing and short:
            year_range = np.arange(max(year_start, year_end - 10), year_end + 1)
        elif not files_missing:
            year_range = np.arange(year_start, year_end)
        else:
            raise ValueError(f"Files are missing in the year range from {year_start}, {year_end}") 
        return year_range