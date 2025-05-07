import glob

# import netCDF4 as nc4
import warnings

warnings.filterwarnings("ignore")

import cartopy.crs as ccrs

from .noresm_abstract_component import NorESMAbstractComponent
from .plotting_methods import make_generic_regridder, regrid_se_data
from  .misc_help_functions import make_regridding_target_from_weightfile

MONTHS = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]

SEASONS = ["DJF", "MAM", "JJA", "SON"]


#def get_minimal_intersecting_year_range(year_range, year_range_other):



class NorESMLndComponent(NorESMAbstractComponent):
    """
    Class that holds and organises info on modelling outputs,
    regridding variables to plot up in diagnostics etc
    """
    def __init__(
        self, datapath,  varlist, weightfile, casename=None
    ):
        super.__init__(f"{datapath}/lnd/hist/", varlist=varlist, casename=casename)
        self.weightfile = weightfile
        self.regridder = make_generic_regridder(self.weightfile, self.filelist[0])
        self.regrid_target = make_regridding_target_from_weightfile(self.weightfile, self.filelist[0])    


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
        super.__init__(f"{datapath}/atm/hist/", varlist=varlist, casename=casename)  


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
            return glob.glob(f"{self.datapath}*.cam.h0i.*.nc")
        
class NorESMOcnComponent(NorESMAbstractComponent):
    """
    Class that holds and organises info on modelling outputs,
    regridding variables to plot up in diagnostics etc
    """
    def __init__(
        self, datapath, varlist, casename=None
    ):
        super.__init__(f"{datapath}/ocn/hist/", varlist=varlist, casename=casename)  


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
        return glob.glob(f"{self.datapath}*.blom.h0.*.nc")
    
class NorESMIceComponent(NorESMAbstractComponent):
    """
    Class that holds and organises info on modelling outputs,
    regridding variables to plot up in diagnostics etc
    """
    def __init__(
        self, datapath, varlist, casename=None
    ):
        super.__init__(f"{datapath}/ice/hist", varlist=varlist, casename=casename)  


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
        super.__init__(f"{datapath}/glc/hist/", varlist=varlist, casename=casename)  


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