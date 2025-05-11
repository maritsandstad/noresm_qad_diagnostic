import glob

# import netCDF4 as nc4
import warnings

warnings.filterwarnings("ignore")

import xarray as xr

from .noresm_abstract_component import NorESMAbstractComponent
from .plotting_methods import make_generic_regridder, regrid_se_data
from  .misc_help_functions import make_regridding_target_from_weightfile

MONTHS = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]

SEASONS = ["DJF", "MAM", "JJA", "SON"]


#def get_minimal_intersecting_year_range(year_range, year_range_other):

_weight_file_dict = {
    48600 : "/datalake/NS9560K/diagnostics/land_xesmf_diag_data/map_ne30pg3_to_0.5x0.5_nomask_aave_da_c180515.nc",
    13824: "/datalake/NS9560K/diagnostics/land_xesmf_diag_data/map_ne16pg3_to_1.9x2.5_nomask_scripgrids_c250425.nc"
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
    

class NorESMAtmComponent(NorESMAbstractComponent):
    """
    Class that holds and organises info on modelling outputs,
    regridding variables to plot up in diagnostics etc
    """
    def __init__(
        self, datapath, varlist, casename=None
    ):
        super().__init__(f"{datapath}/atm/hist/", varlist=varlist, casename=casename)  


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
        
class NorESMOcnComponent(NorESMAbstractComponent):
    """
    Class that holds and organises info on modelling outputs,
    regridding variables to plot up in diagnostics etc
    """
    def __init__(
        self, datapath, varlist, casename=None
    ):
        super().__init__(f"{datapath}/ocn/hist/", varlist=varlist, casename=casename)  


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
        return glob.glob(f"{self.datapath}*.blom.hd.*.nc")
    
    def get_weights():
        # TODO: Get areacella weights here
        pass
    
class NorESMIceComponent(NorESMAbstractComponent):
    """
    Class that holds and organises info on modelling outputs,
    regridding variables to plot up in diagnostics etc
    """
    def __init__(
        self, datapath, varlist, casename=None
    ):
        super().__init__(f"{datapath}/ice/hist/", varlist=varlist, casename=casename)

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