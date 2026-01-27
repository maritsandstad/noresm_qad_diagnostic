import os, sys, glob, shutil, grp, stat


def move_file_from_folder_ens_member(
        root_target = "/datalake/NS9560K/diagnostics/noresm_qad_data/ppe_runs_partII", 
        root_orig = "figs/ppe_runs_partII",
        ens_member = 0,
        subfolder_name = "netcdf_dumps",
        fname_start = "",
        copy = True
    ):
    ens_member_name = f"ensemble_member.{ens_member:03d}"
    if not os.path.exists(f"{root_target}/{ens_member_name}/{subfolder_name}"):
        os.makedirs(f"{root_target}/{ens_member_name}/{subfolder_name}")
    print(f"{root_orig}/{ens_member_name}/{subfolder_name}/{fname_start}*")
    for file in glob.glob(f"{root_orig}/{ens_member_name}/{subfolder_name}/{fname_start}*"):
        fname = file.split("/")[-1]
        target_path = f"{root_target}/{ens_member_name}/{subfolder_name}/{fname}"
        print(f"Now moving {file} to {target_path}")

        if copy:
            shutil.copy(file, target_path)
        else: 
            shutil.move(file, target_path)
        os.chown(target_path, os.getuid(), grp.getgrnam("ns9560k").gr_gid)
        os.chmod(target_path, stat.S_IRWXU | stat.S_IRWXG)

ensmembers = [i for i in range(61)]#[ 0, 3, 5, 9, 19, 21, 24, 30, 32, 39, 41, 42, 48, 49, 58, 59, 64, 65, 69, 72, 73]
#seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
#ensmembers = [19, 21, 32, 65, 69]
for ens_member in ensmembers:
    print(ens_member)
    move_file_from_folder_ens_member(ens_member=ens_member, subfolder_name=f"")#, fname_start="FLNS")