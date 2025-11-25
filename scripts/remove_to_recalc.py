import sys, os, glob

keep = ["E-P", "PRECC","QFLX", "PRECL", "PRECT", "co2fxu", "epc100", "co2fxd", "ppint", "fgco2"]
for folder in glob.glob("figs/ppe_runs/ensemble_member.???"):
    for file in glob.glob(f"{folder}/*.nc"):
        varname = file.split("/")[-1].split("_")[0]
        if varname in keep:
            continue
        #print(f"Now deleting {file}")
        os.remove(file)
    #sys.exit(4)

