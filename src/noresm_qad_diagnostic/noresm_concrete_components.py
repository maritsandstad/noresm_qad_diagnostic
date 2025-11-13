import glob, sys

# import netCDF4 as nc4
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import xarray as xr

from .setup_loging import get_logger
from .noresm_abstract_component import NorESMAbstractComponent
from .plotting_methods import make_generic_regridder, regrid_se_data
#from .regrid_functions import make_regridder, regrid_file
from .misc_help_functions import make_regridding_target_from_weightfile, get_unit_conversion_from_string
from .general_util_functions import amoc
from .regrid_functions import make_regridder, make_regular_outgrid
# setup logging
logger = get_logger(__name__)

MONTHS = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]

SEASONS = ["DJF", "MAM", "JJA", "SON"]


#def get_minimal_intersecting_year_range(year_range, year_range_other):

_weight_file_dict = {
    48600 : "/nird/datalake/NS9560K/diagnostics/land_xesmf_diag_data/map_ne30pg3_to_0.5x0.5_nomask_aave_da_c180515.nc",
    13824: "/nird/datalake/NS9560K/diagnostics/land_xesmf_diag_data/map_ne16pg3_to_1.9x2.5_nomask_scripgrids_c250425.nc"
}

class NorESMLndComponent(NorESMAbstractComponent):
    """
    Class that holds and organises info on modelling outputs,
    regridding variables to plot up in diagnostics etc
    """
    def __init__(
        self, datapath,  varlist, weightfile=None, casename=None
    ):
        super().__init__(f"{datapath}/lnd/hist/", varlist=varlist, casename=casename)
        self._set_correct_weightfile()
        self.regridder = make_generic_regridder(self.weightfile, self.filelist[0])
        self.regrid_target = make_regridding_target_from_weightfile(self.weightfile, self.filelist[0])


    #def clean_out_empty_folders(self):
    #    clean_empty_folders_in_tree(self.outdir)

    def _set_correct_weightfile(self):
        first_file = self.filelist[0]
        data_check = xr.open_dataset(first_file)
        if "lon" in data_check.dims and "lat" in data_check.dims:
            self.weightfile = None
            return
        if "lndgrid" in data_check.dims:
            print(len(data_check["lndgrid"]))
            if len(data_check["lndgrid"]) in _weight_file_dict:
                self.weightfile = _weight_file_dict[len(data_check["lndgrid"])]
                self.regrid_target = make_regridding_target_from_weightfile(None, self.filelist[0])
                return

        print("No weightfile found for land data, but one might be needed")
        self.weightfile = None

    def get_filelist(self):
        """
        Get a list of all the files for the case

        Returns
        -------
        list
            with paths to all the clm2.h0 files
        """
        #(f"{self.datapath}*.clm2.h0.*.nc")
        return glob.glob(f"{self.datapath}*.clm2.h0.*.nc")

    def regrid(self, data):
        return regrid_se_data(self.regridder, data)

    def get_weights_regrid(self, outd_here):
        outd_here = self.regrid(outd_here)
        return outd_here, np.cos(np.deg2rad(outd_here.lat))



class NorESMAtmComponent(NorESMAbstractComponent):
    """
    Class that holds and organises info on modelling outputs,
    regridding variables to plot up in diagnostics etc
    """
    def __init__(
        self, datapath, varlist, casename=None
    ):
        super().__init__(f"{datapath}/atm/hist/", varlist=varlist, casename=casename)
        if "toa" in self.varpams:
            self.unit_dict["toa"] = self.unit_dict["FSNT"]


    #def clean_out_empty_folders(self):
    #    clean_empty_folders_in_tree(self.outdir)

    def get_filelist(self):
        """
        Get a list of all the files for the case

        Returns
        -------
        list
            with paths to all the clm2.h0 files
        """
        #(f"{self.datapath}*.clm2.h0.*.nc")
        if len(glob.glob(f"{self.datapath}*.cam.h0.*.nc")) > 0:
            return glob.glob(f"{self.datapath}*.cam.h0.*.nc")
        else:
            return glob.glob(f"{self.datapath}*.cam.h0a.*.nc")

    def set_composite_variable_dict(self):
        self.composite_variable_dict = {"toa": ["FSNT", "FLNT"]}

    def make_spatial_means_do_unit_fixes_etc(self, weighted_data, spatial_coords):
        mean_ts = weighted_data.mean(spatial_coords)
        if "toa" in self.varpams:
            mean_ts["toa"] = mean_ts["FSNT"] - mean_ts["FLNT"]
        return mean_ts

class NorESMOcnComponent(NorESMAbstractComponent):
    """
    Class that holds and organises info on modelling outputs,
    regridding variables to plot up in diagnostics etc
    """
    def __init__(
        self, datapath, varlist, casename=None, weightfile = None
    ):
        super().__init__(f"{datapath}/ocn/hist/", varlist=varlist, casename=casename)

        if weightfile is None:
            self.weightfile = ("/nird/datalake/NS9560K/diagnostics/land_xesmf_diag_data/grid_tnx1v4_20170622.nc")
        self.regrid_target = make_regular_outgrid()
        self.unit_dict["AMOC"] = "Sv"

        self.regridder = make_regridder(self.weightfile, self.regrid_target)
    #def clean_out_empty_folders(self):
    #    clean_empty_folders_in_tree(self.outdir)

    def get_filelist(self):
        """
        Get a list of all the files for the case

        Returns
        -------
        list
            with paths to all the clm2.h0 files
        """
        #(f"{self.datapath}*.clm2.h0.*.nc")
        return glob.glob(f"{self.datapath}*.blom.hm.*.nc")

    def set_composite_variable_dict(self):
        self.composite_variable_dict = {"AMOC": ["mmflxd"]}

    def get_weights_regrid(self, outd_here):
        # TODO: Get areacella weights here
        grid = xr.open_dataset(self.weightfile)
        pweight = grid.parea * grid.pmask.where(grid.pmask > 0)
        pweight = pweight.fillna(0)
        # add latitude and longitude info good to have in case of e.g. regridding
        pweight = pweight.assign_coords(lat=grid.plat)
        pweight = pweight.assign_coords(lon=grid.plon)
        return outd_here, pweight

    def add_variables_pre_regridding(self, outd, varlist_in=None):
        if varlist_in is None:
            varlist_in = self.varpams
        if "AMOC" in varlist_in:
            if len(varlist_in) == 1:
                outd["AMOC"] = amoc(outd)
            else:
                outd = xr.merge([outd, amoc(outd)])
        return outd

    def regrid(self, data):
        return self.regridder(data)

    #def regrid(self, data):
    #    return self.regridder(data)

class NorESMOcnbgcComponent(NorESMOcnComponent):

    def __init__(self, datapath, varlist, casename=None, weightfile=None):
        super().__init__(datapath, varlist, casename, weightfile)
        if "fgco2" in self.varpams:
            self.unit_dict["fgco2"] = self.unit_dict["co2fxu"]
        self.expected_unit = {
            "fgco2": "PgC y-1",
            "ppint": "PgC y-1",
            "epc100": "PgC y-1",
        }
        #self._run_unit_conversion_test()
        #sys.exit(4)

    def _run_unit_conversion_test(self):
        for variable in self.varpams:
            if variable in self.expected_unit.keys():
                orig_unit = self.unit_dict[variable]
                expected_unit = self.expected_unit[variable]
                #print(f"Attempting unit conversion for {variable} with orig_unit {orig_unit} and expected unit {expected_unit}")
                #print(f"Getting conversion factor{get_unit_conversion_from_string(expected_unit, orig_unit, areasummed=True)}")

    def get_variable_unit(self, variable):
        if variable in self.expected_unit.keys():
            return self.expected_unit[variable]
        return super().get_variable_unit(variable)

    def get_compname(self):
        return "ocnbgc"

    def get_filelist(self):
        """
        Get a list of all the files for the case

        Returns
        -------
        list
            with paths to all the clm2.h0 files
        """
        #(f"{self.datapath}*.clm2.h0.*.nc")
        return glob.glob(f"{self.datapath}*.blom.hbgcm.*.nc")


    def set_composite_variable_dict(self):
        self.composite_variable_dict = {"fgco2": ["co2fxu", "co2fxd"]}

    def make_spatial_means_do_unit_fixes_etc(self, weighted_data, spatial_coords):
        mean_ts = weighted_data.sum(spatial_coords)
        if "fgco2" in self.varpams:
            mean_ts["fgco2"] = mean_ts["co2fxd"] - mean_ts["co2fxu"]
        for variable in mean_ts.data_vars:
            if variable in self.expected_unit.keys():
                conv_factor, new_unit =  get_unit_conversion_from_string(self.expected_unit[variable], self.unit_dict[variable], areasummed=True)
                #print(variable)
                #print(mean_ts[variable])
                mean_ts[variable] = mean_ts[variable] * conv_factor

        return mean_ts


class NorESMIceComponent(NorESMAbstractComponent):
    """
    Class that holds and organises info on modelling outputs,
    regridding variables to plot up in diagnostics etc
    """
    def __init__(
        self, datapath, varlist, casename=None, weightfile=None
    ):
        super().__init__(f"{datapath}/ice/hist/", varlist=varlist, casename=casename)
        self.weight_dict = {}
        self.regrid_target = make_regular_outgrid()
        if weightfile is None:
            self.weightfile = ("/nird/datalake/NS9560K/diagnostics/land_xesmf_diag_data/grid_tnx1v4_20170622.nc")
        self.precalc_weight_dict()
        if "aice-antarctic" in self.varpams:
            self.unit_dict["aice-antarctic"] = self.unit_dict["aice"]
        if "aice-antarctic" in self.varpams:
            self.unit_dict["aice-arctic"] = self.unit_dict["aice"]
        self.regridder = make_regridder(self.weightfile, self.regrid_target)
        #self.regridder = make_regridder(self.get_filelist[0], _weight_file_dict[13824], seaice=True)
    #def clean_out_empty_folders(self):
    #    clean_empty_folders_in_tree(self.outdir)

    def get_filelist(self):
        """
        Get a list of all the files for the case

        Returns
        -------
        list
            with paths to all the clm2.h0 files
        """
        #(f"{self.datapath}*.clm2.h0.*.nc")
        return glob.glob(f"{self.datapath}*.cice.h.*.nc")

    def set_composite_variable_dict(self):
        self.composite_variable_dict = {
            "aice-arctic": ["aice", "TLAT", "TLON"],
            "aice-antarctic": ["aice", "TLAT", "TLON"]
            }

    def calc_composite_components(self, outd, varlist_in):
        return self.get_hemispheric_split_variables(outd, varlist_in)

    def set_hemisplit(self, datause):
        full_len = len(datause.nj)
        test = datause.where(datause > -1e-5, drop= True)
        len_nh = len(test.nj)
        self.hemisplit_nj = full_len -len_nh

    def add_variables_pre_regridding(self, outd, varlist_in=None):
        return self.get_hemispheric_split_variables(outd, varlist_in=varlist_in)


    def precalc_weight_dict(self):
        vars_for_weights = list(set(self.varpams.copy()) - set(self.composite_variable_dict.keys()))

        ext_varlist = vars_for_weights.copy()
        for areatype in ["tarea", "uarea", "narea","earea", "TLAT"]:
            ext_varlist.append(areatype)
        datause = xr.open_dataset(self.filelist[0])[ext_varlist]
        self.set_hemisplit(datause["TLAT"])
        for var in vars_for_weights:
            if "cell_measures" in datause[var].attrs.keys():
                #print(datause[var].attrs["cell_measures"].split(" ")[-1])
                self.weight_dict[var] = datause[datause[var].attrs["cell_measures"].split(" ")[-1]]

    def get_weights_regrid(self, outd_here):
        for key, weight in self.weight_dict.items():
            if key in outd_here.keys():
                return outd_here, weight
        return outd_here, super().get_weights(outd_here)

    ##def regrid(self, data):
        return self.regridder(data)

    def get_hemispheric_split_variables(self, orig_data, varlist_in = None):
        if varlist_in is None:
            varlist_in = self.varpams
        if "aice-antarctic" in varlist_in:
            orig_data["aice-antarctic"] = orig_data["aice"].where(orig_data.nj<self.hemisplit_nj, 0)
        if "aice-arctic" in varlist_in:
            orig_data["aice-arctic"] = orig_data["aice"].where(orig_data.nj >= self.hemisplit_nj, 0)
        return orig_data[varlist_in]

    def regrid(self, data):
        regridded = super().regrid(data)
        #egridded = self.get_hemispheric_split_variables(regridded)
        return regridded

    def make_spatial_means_do_unit_fixes_etc(self, weighted_data, spatial_coords):
        mean_ts = weighted_data.sum(spatial_coords)
        return mean_ts

    def make_comparison_map_no_comparison_data(self, variable, subfig_handle):
        #Put in single map-plot
        logger.warning(f"No observational data for {variable}")
        axs = subfig_handle.subplots(
            nrows=1,
            ncols=1,
            subplot_kw={"projection": ccrs.Robinson()},
            )
        plot_data = self.main_run.components[self.var_to_comp_map[variable]].get_annual_mean_map_data(varlist_in=variable)
        make_bias_plot(plot_data, self.casename, ax = axs)


class NorESMGlcComponent(NorESMAbstractComponent):
    """
    Class that holds and organises info on modelling outputs,
    regridding variables to plot up in diagnostics etc
    """
    def __init__(
        self, datapath, varlist, casename=None
    ):
        super().__init__(f"{datapath}/glc/hist/", varlist=varlist, casename=casename)


    #def clean_out_empty_folders(self):
    #    clean_empty_folders_in_tree(self.outdir)

    def get_filelist(self):
        """
        Get a list of all the files for the case

        Returns
        -------
        list
            with paths to all the clm2.h0 files
        """
        #(f"{self.datapath}*.clm2.h0.*.nc")
        return glob.glob(f"{self.datapath}*.glc.h.*.nc")