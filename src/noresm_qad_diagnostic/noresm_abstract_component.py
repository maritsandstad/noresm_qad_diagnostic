import numpy as np

import xarray as xr

import os, glob, sys

# import netCDF4 as nc4
import warnings

from abc import ABC, abstractmethod

warnings.filterwarnings("ignore")

from .plotting_methods import make_generic_regridder, regrid_se_data, make_bias_plot
from .infrastructure_help_functions import setup_nested_folder_structure_from_dict, get_spatial_coordinates#, clean_empty_folders_in_tree
from  .misc_help_functions import get_unit_conversion_and_new_label, make_regridding_target_from_weightfile, get_unit_conversion_from_string, do_light_unit_string_conversion
from .general_util_functions import yearly_avg, DAYS_PER_MONTH


#def get_minimal_intersecting_year_range(year_range, year_range_other):



class NorESMAbstractComponent(ABC):
    """
    Class that holds and organises info on modelling outputs from some component,
    Methods on getting data on various formats etc.
    """
    def __init__(
        self, datapath, varlist, casename=None,
    ):
        self.datapath = datapath
        #self.var_pams = read_pam_file(pamfile)
        self.filelist = self.get_filelist()
        self.filelist.sort()
        self.varpams = varlist
        print(self.varpams)
        #print(self.filelist)
        if casename is None:
            self.casename = ".".join(self.filelist[0].split("/")[-1].split(".")[:-4])
        else:
            self.casename = casename
        self.unit_dict = {}
        self.set_composite_variable_dict()
        vars_missing = self.wet_and_cut_varlists()
        print(self.varpams)

        if len(vars_missing):
            print("Not all requested variables are available in output, ignoring these:")
            print(vars_missing)
        self.add_to_unit_dict(self.varpams)
    
    def set_composite_variable_dict(self):
        self.composite_variable_dict = {}

    def wet_and_cut_varlists(self):
        # TODO: Make sure all items in SEASONAL is also in var_list_main 
        read = xr.open_dataset(self.filelist[0]).keys()
        #lists_check =["VAR_LIST_MAIN", "COMPARE_VARIABLES"]

        vars_missing = []
        vars_needed_for_composites = []
        for item in self.varpams:
            if item not in read:
                if item in self.composite_variable_dict:
                    vars_needed_for_composites.extend(self.composite_variable_dict[item])
                else:
                    vars_missing.append(item)
        self.varpams = list(set(self.varpams) - set(vars_missing))
        self.varpams = list(set(self.varpams).union(set(vars_needed_for_composites)))

        return vars_missing       
        #for list_n in lists_check:
        #    for item in self.varpams[list_n]:
        #        if item not in read:
        #            vars_missing.append(item)
        #    self.varpams[list_n] = list(set(self.varpams[list_n]) - set(vars_missing))
        #return vars_missing

    def add_to_unit_dict(self, varlist):
        missing = list(set(varlist) - set(self.unit_dict.keys()))
        if len(missing) < 1:
            return
        read = xr.open_dataset(self.filelist[0])
        for vrm in missing:
            if vrm in read.keys():
                if "units" in read[vrm].attrs.keys():
                    self.unit_dict[vrm] = do_light_unit_string_conversion(read[vrm].attrs["units"])
                else:
                    self.unit_dict[vrm] = "No unit"

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
    
    def get_file_str_for_regex(self):
        filestr_for_regex = ".".join(self.filelist[0].split(".")[-4:-2])
        return filestr_for_regex
    
    def get_file_str_ending(self):
        fileending = self.filelist[0].split(".")[-2]
        if len(fileending) <= 7:
            return ""
        return fileending[7:]

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

    def get_annual_mean_map_data(self, year_range=None, varlist=None):
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
        if year_range is None:
            year_range = self.get_year_range(short=True)
        outd = None
        if varlist is None:
            varlist = self.varpams
        for year in year_range:
            for month in range(12):
                mfile = f"{self.datapath}/{self.casename}.{self.get_file_str_for_regex()}.{year:04d}-{month + 1:02d}{self.get_file_str_ending()}.nc"
                outd_here = xr.open_dataset(mfile, engine="netcdf4")[varlist]
                # print(outd_here)
                # sys.exit(4)
                if outd is None:
                    outd = outd_here
                else:
                    outd = xr.concat([outd, outd_here], dim="time")
        outd = outd.mean(dim="time")
        return self.regrid(outd)
    
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
        factors = [DAYS_PER_MONTH["noleap"][month+1] / np.sum(DAYS_PER_MONTH["noleap"]) for month in range(12)]
        if varlist is None:
            varlist = self.varpams["VAR_LIST_MAIN"]
        for year in year_range:
            outd_yr = None
            for month in range(12):                         
                mfile = f"{self.datapath}/{self.casename}.{self.get_file_str_for_regex()}.{year:04d}-{month + 1:02d}{self.get_file_str_ending()}.nc"
                outd_here = xr.open_dataset(mfile, engine="netcdf4")[varlist]
                # print(outd_here)
                # sys.exit(4)
                if not outd_yr:
                    outd_yr = outd_here * factors[month]
                else:
                    outd_yr = xr.concat([outd_yr, outd_here * factors[month]], dim="time")
            outd_yr = outd_yr.sum(dim="time")
            if not outd:
                outd = outd_yr
            else: 
                outd = xr.concat([outd, outd_yr], dim="time")
        return outd
    
    def get_area_mean_ts_data(self, varlist = None, year_range=None, area_def = None):
        outd = None
        outd_year = None
        print(varlist)
        print(self.varpams)
        factors = [DAYS_PER_MONTH["noleap"][month+1] / np.sum(DAYS_PER_MONTH["noleap"]) for month in range(12)]
        if varlist is None:
            varlist = self.varpams["VAR_LIST_MAIN"]
        if year_range is None:
            year_range = self.get_year_range()
        weights = None
        spatial_coords = None
        for year in year_range:
            outd_yr_here = None

            mfile = f"{self.datapath}/{self.casename}.{self.get_file_str_for_regex()}.{year:04d}-*.nc"
            outd_here = xr.open_mfdataset(mfile, engine="netcdf4")[varlist]
            #print(outd_here)
            #sys.exit(4)
            #for month in range(12):

             #   mfile = f"{self.datapath}/{self.casename}.{self.get_file_str_for_regex()}.{year:04d}-{month + 1:02d}{self.get_file_str_ending()}.nc"
            #outd_here = xr.open_mfdataset(mfile, engine="netcdf4")[varlist]
            outd_here_t_weighted = outd_here.weighted(xr.DataArray(factors, dims=("time",)))
            outd_here = outd_here_t_weighted.sum(dim="time")
            
            outd_here = self.regrid(outd_here)
            #print(outd_here)
            if weights is None:
                weights = self.get_weights(outd_here)
            if spatial_coords is None:
                get_spatial_coordinates(outd_here)
            if area_def is not None:
                outd_here = outd_here.sel(
                lat=slice(area_def["lat_s"], area_def["lat_n"]),
                lon=slice(area_def["lon_w"], area_def["lon_e"]),
            )                      
                
            weighted_data = outd_here.weighted(weights)
            ts_data = weighted_data.mean(spatial_coords)
            #print(spatial_coords)
            #print(len(weighted_data.mean(spatial_coords)))
            #print(len(factors))
            #print(ts_data.coords)
            #ts_weight = ts_data.weighted(xr.DataArray(factors, dims=("time",), coords = ts_data.coords))
            #if not outd_yr_here:
                #outd_yr_here = ts_data * factors
            #else:
                #outd_yr_here = xr.concat([outd_yr_here, ts_data*factors[month]], dim="time")
            #ts_data = ts_weight.sum(dim="time")
            if not outd_year:
                outd_year = ts_data
            else: 
                outd_year = xr.concat([outd_year, ts_data], dim="time")
        return outd_year
    
    def get_weights(self, outd_here):
        return np.cos(np.deg2rad(outd_here.lat))
    
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
            varlist = self.varpams["VAR_LIST_MAIN"]
        for year in year_range:
            for monthincr in range(3):

                month = monthincr + season * 3
                if month == 0:
                    month = 12
                # print(f"Season: {season}, monthincr: {monthincr}, month: {monthincr}")
                mfile = (
                    f"{self.datapath}/{self.casename}.{self.get_file_str_for_regex()}.{year:04d}-{month:02d}{self.get_file_str_ending()}.nc"
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
                mfile = f"{self.datapath}/{self.casename}.{self.get_file_str_for_regex()}.{year:04d}-{month+1:02d}{self.get_file_str_ending()}.nc"
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
            year_range = np.arange(max(year_start, year_end - 20), year_end + 1)
        elif not files_missing:
            year_range = np.arange(year_start, year_end)
        else:
            raise ValueError(f"Files are missing in the year range from {year_start}, {year_end}") 
        return year_range
    
    def find_case_year_range(self):
        year_start = int(self.filelist[0].split(".")[-2].split("-")[0])
        year_end = int(self.filelist[-1].split(".")[-2].split("-")[0])
        files_missing = False
        if len(self.filelist) < (year_end - year_start) * 12:
            files_missing = True
        return year_start, year_end, files_missing