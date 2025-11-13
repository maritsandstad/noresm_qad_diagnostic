import numpy as np

import xarray as xr

import matplotlib.pyplot as plt

import json
import os, sys, glob

# import netCDF4 as nc4
import warnings

warnings.filterwarnings("ignore")

from .setup_logging import get_logger
from .noresm_concrete_components import NorESMAtmComponent, NorESMLndComponent, NorESMOcnComponent, NorESMOcnbgcComponent, NorESMGlcComponent, NorESMIceComponent
# setup logging
logger = get_logger(__name__)

# Ocean might need to be multiple components...
DEFAULT_VAR_LIST_MAIN = {
    "atm": ["TREFHT", "DOD550"],
    "ocn": ["sst", "sss"],
    "ice": ["aice"],
    "lnd": ["nbp", "TOTSOILC", "TOTVEGC", "FATES_GPP", "FATES_LAI"]
}

KEY_COMP_MAPPING = {
    "atm" : NorESMAtmComponent,
    "lnd" : NorESMLndComponent,
    "ocn" : NorESMOcnComponent,
    "ocn-bgc": NorESMOcnbgcComponent,
    "ice": NorESMIceComponent,
    "glc": NorESMGlcComponent
}

class NorESMFullModel:

    def __init__(self, datapath, pamfile=None, casename=None):
        self.datapath = datapath
        self.casename = casename
        self._set_components(pamfile=pamfile)
        if self.casename is None:
            for compname, comp in self.components.items():
                if comp.casename is not None:
                    self.casename = comp.casename
                    break

    def _set_components(self, pamfile):
        self.components = {}
        if pamfile is None:
            pams = {'VAR_LIST_MAIN': DEFAULT_VAR_LIST_MAIN}
        elif isinstance(pamfile, dict):
            pams = pamfile
        else:
            with open(pamfile, 'r') as jsonfile:
                pams = json.load(jsonfile)
        if "VAR_LIST_MAIN" not in pams:
                raise ValueError("pamfile dictionary must contain 'VAR_LIST_MAIN' key.")
        for key, values in pams["VAR_LIST_MAIN"].items():
            logger.info(f"Now initialising {key}")
            #Deal with weighting for land:
            self.components[key] = KEY_COMP_MAPPING[key](self.datapath, values, casename = self.casename)