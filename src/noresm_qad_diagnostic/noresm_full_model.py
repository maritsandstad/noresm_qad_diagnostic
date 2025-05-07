import numpy as np

import xarray as xr

import matplotlib.pyplot as plt

import json
import os, sys, glob

# import netCDF4 as nc4
import warnings

warnings.filterwarnings("ignore")

from .noresm_concrete_components import NorESMAtmComponent, NorESMLndComponent, NorESMOcnComponent, NorESMGlcComponent, NorESMIceComponent

DEFAULT_PAMS = {
    "atm": ["TREFHT", "DOD550"],
    "ocn": ["sst", "sss"],
    "ice": ["aice"], 
    "lnd": ["nbp", "TOTSOILC", "TOTVEGC", "FATES_GPP", "FATES_LAI"]
}

KEY_COMP_MAPPING = {
    "atm" : NorESMAtmComponent,
    "lnd" : NorESMLndComponent,
    "ocn" : NorESMOcnComponent,
    "ice": NorESMIceComponent,
    "glc": NorESMGlcComponent
}

class NorESMFullModel:

    def __init__(self, datapath, pamfile=None):
        self.datapath = datapath
        self._set_components()

    def _set_components(self, pamfile):
        self.components = {}
        if pamfile is None:
            pams = DEFAULT_PAMS
        elif isinstance(pamfile, dict):
            pams = pamfile
        else:
            with open(pamfile, 'r') as jsonfile:
                pams = json.load(jsonfile)
        for key, values in pams.items():
            #Deal with weighting for land:
            self.components[key] = KEY_COMP_MAPPING[key](self.datapath, values)

        