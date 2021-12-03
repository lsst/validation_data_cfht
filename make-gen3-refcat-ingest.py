"""Generate an astropy-readable .ecsv files for `butler ingest-files`, to ingest an existing gen2 refcat.

The `refcat_dir` variable needs to be modified for each refcat being converted.
"""
import os
import glob
import astropy.table

refcat_dir = "ref_cats/gaia_dr2_20200414"
out_dir = "."

out_file = f"{out_dir}/{os.path.basename(refcat_dir)}.ecsv"

table = astropy.table.Table(names=("filename", "htm7"), dtype=("str", "int"))
files = glob.glob(f"{refcat_dir}/[0-9]*.fits")

for i, file in enumerate(files):
    # running status, overwriting each print statement as it proceeds
    print(f"{i}/{len(files)} ({100*i/len(files):0.1f}%)", end="\r")

    # extract file index; add row to table
    file_index = int(os.path.basename(os.path.splitext(file)[0]))
    table.add_row((file, file_index))

table.write(out_file)
print(f"Saving to: {out_file}")
