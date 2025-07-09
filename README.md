# noresm_qad_diagnostic

A quick and dirty diagnostic tool for NorESM to produce various overview plots and tables

Currently very much a work in progress based on xesmf land diagnostic,

NB: This README is currently not correct for this diagnostic at all, will be updated

## Prerequisites

In order to use the tool you need to load and ESMF module and build an xesmf-containing conda environment on top of it. On Nird running, navigating to the folder called `scripts` and running:

```
. setup.sh
```
will automatically do this for you.

Otherwise you can install an environment using the file `conda_env.yaml` to install such an environment yourself using conda (or Miniforge). However, be aware that to do so you also need ESMF installed and loaded before building, and you need to loaded the same ESMF module when loading the environment. We recommend editing the `setup.sh` file to reflect the setup needed on your machine if you build your own environment this way.

## Usage

To run navigate to the folder called scripts (or extend paths for run-scripts to include the full path to that folder in the following commands) and run: 
```
. setup.sh
```
Then run 
```
python test_setup.py path_1 weight=weight_path outpath=opt_out_path pamfile=pamfile_path
```
where `path_1` is the path to the lnd/hist folder containing your output.

The other arguments are optional:

`weight_path` is a path to a weight-file if the standard one is not to be used. For regular lat, lon runs, the weight-file is not used, and hence you can send a dumy argument, use the standard, or if not working on nird, the tool will try to add a dummy and run without if you don't send a weight_path. The standard weight path is one for an ne30 grid, so for this you also won't need to provide a weight_path. For an ne16 grid, this path will work: /datalake/NS9560K/diagnostics/land_xesmf_diag_data/map_ne16pg3_to_1.9x2.5_nomask_scripgrids_c250425.nc

`outpath` is the path of where you want the output diagnostic figures filetree to go. If not sent the figures will be expected to go in a folder called figs situated in whatever directory you ran the command from.
If you want this to be viewable by web, choose a web-facing directory. For instance if you have access to the NS9560K account, make a subdirectory with the same name as your username in /datalake/NS9560K/www/diagnostics/noresm/ and make that your outpath, i.e. `outpath=/datalake/NS9560K/www/diagnostics/noresm/username`

`pamfile_path` is the path to a parameterfile in which you can specify which variables to plot in the various plots. For most runs we recommend just running with the default file (i.e. no need to send this as an argument, the standard pams_tes.json will be used), to get standard plots. The
This file should be a json-file containing the three (four) keyword arguments:
* VAR_LIST_MAIN - and this should be followed by the list of main variables to plot on maps and for trends
* OBS_COMPARE_MAPS - this should be followed by a dictionary of named variable sets for which to plot map comparison to observations as keys and items should be lists with observation names that exist in the ilamb_qad-diags.cfg file
* COMPARE_RUNS - is optional and where you might possibly want to make changes from the default parameterfile to include or omit runs. It should be formatted as a dictionary with a top folder path as key, and a list of the casenames for the runs to be used from that folder as the dictionary item.