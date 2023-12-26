import json
import subprocess
import os

output_dir = '/home/tdpham/work/research/MOF_databases/Anionic_MOFs_project'
input_dir = '/home/tdpham/work/research/MOF_databases/CSD_nondisorder'

with open('../data/Anionic_MOFs_formula.json', 'r') as f:
    data = json.load(f)

for mof in data:
    mof_path = os.path.join(input_dir, mof + '.Non-disordered_MOF_subset.cif')
    if os.path.isfile(mof_path):
        subprocess.run("cp {} {}".format(mof_path, output_dir), shell=True)
