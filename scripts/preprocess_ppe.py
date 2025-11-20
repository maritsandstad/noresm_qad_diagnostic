import os
import sys
import glob
import json

import numpy as np
import matplotlib as mpl
sys.path.append(os.path.join(os.path.dirname(__file__), "../", "src"))

from noresm_qad_diagnostic import noresm_qad_diagnostic, ilamb_configurations, noresm_qad_preprocessor
#mpl.rc('font', size=30) # Set default font size to 20

standard_run_dict = {
    "weight" : "/datalake/NS9560K/diagnostics/land_xesmf_diag_data/map_ne16pg3_to_1.9x2.5_nomask_scripgrids_c250425.nc",
    "outpath" : "figs/ppe_runs/",
    "pamfile" : f"{os.path.dirname(__file__)}/pam_test.json",
    "compare": None,
    "year_range_compare": None, 
}

ppe_path = "/datalake/NS9560K/noresm3/cases/coupled_ppe.20251108/"
runs_dict = {}
for member_path in glob.glob(f"{ppe_path}ensemble*"):
    run_name = member_path.split("/")[-1]
    print(run_name)
    runs_dict[run_name] = member_path

#runs_dict = {"n1850.ne16pg3_tn14.ppe_base_run.nor3_b04._202511105": "/nird/datalake/NS9560K/noresm3/cases/n1850.ne16pg3_tn14.ppe_base_run.nor3_b04._202511105/"}

print(runs_dict)
with open("pam_test_preproc.json", "r") as cfile:
    config_dict = json.load(cfile)

preprocessor = noresm_qad_preprocessor.NorESMQADPreprocessor(runs_dict=runs_dict, config_dict= config_dict,  paths_preproc=standard_run_dict["outpath"])