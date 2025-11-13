import os, sys, glob, grp
import json
import stat

from .setup_loging import get_logger

# Set up logger for this module
logger = get_logger(__name__)

def setup_folder_if_not_exists(path):
    """
    Create a folder if it does not exist already

    Parameters
    ----------
    path : str
        Path to folder that should be created

    """
    # TODO: Throw error
    if not os.path.exists(path):
        os.mkdir(path)
        os.chown(path, os.getuid(), grp.getgrnam("ns9560k").gr_gid)
        os.chmod(path, stat.S_IRWXU | stat.S_IRWXG)

def setup_nested_folder_structure_from_dict(root, subfolder_dict):
    setup_folder_if_not_exists(root)
    if isinstance(subfolder_dict, dict):
        for key, value in subfolder_dict.items():
            setup_nested_folder_structure_from_dict(f"{root}/{key}", value)
    elif isinstance(subfolder_dict, list):
        for value in subfolder_dict:
            setup_folder_if_not_exists(f"{root}/{value}")
    elif isinstance(subfolder_dict, str):
        setup_folder_if_not_exists(f"{root}/{subfolder_dict}")
    return

def clean_empty_folders_in_tree(root):
    empty_below = 0
    for sub in glob.glob(root):
        if os.path.isdir(sub):
            empty_below = empty_below + clean_empty_folders_in_tree(sub)
        else:
            empty_below = empty_below + 1
    if empty_below == 0:
        os.rmdir(root)
        return 0
    return empty_below

def read_pam_file(pam_file_path):
    required = {"VAR_LIST_MAIN": dict,"OBS_COMPARE_MAPS":dict}
    try:
        with open(pam_file_path, "r") as jsonfile:
            data = json.load(jsonfile)
    except json.decoder.JSONDecodeError as err:
        logger.error(f"{pam_file_path} must be a valid json-file, error decoding: {err}")
        sys.exit(4)
    if not isinstance(data, dict):
        raise ValueError(f"{pam_file_path} must evaluate to dict")
    for elem, e_type in required.items():
        if elem not in data.keys():
            if (elem == "OBS_COMPARE_MAPS"):
                data["OBS_COMPARE_MAPS"] = None
                continue
            else:
                raise ValueError(f"{pam_file_path} must include {elem}")
        if not isinstance(data[elem], e_type):
            raise TypeError(f"{pam_file_path} element {elem} must be a {e_type}, but is {type(data[elem])}")
    return data

def make_monthly_taxis_from_years(years):
    return [years[0] + tstep/12. for tstep in range(12*len(years))]

def get_spatial_coordinates(spatial_data):
    if "lat" and "lon" in spatial_data.keys():
        return ["lat", "lon"]
    if "x" and "y" in spatial_data.keys():
        return ["x", "y"]
    if "nj" and "ni" in spatial_data.keys():
        return ["nj", "ni"]




