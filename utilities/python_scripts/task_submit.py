# read json files

from ase.io.jsonio import encode, decode
from ase.io import write
import json
import os

project_name = input('Enter the file name of your current project with JSON files and the code: ')
os.chdir("/ocean/projects/mat210008p/jhuanga/project/%s" % project_name)
os.system('mkdir -p task_all')
os.system('mkdir -p task_pretest')

json_filename = input("Please enter the file name of your JSON file: ")
with open(json_filename, "r") as f:
    out = json.load(f)

ase_list = []
for values in out.values():
    ase_struc = decode(encode(values))
    ase_list.append(ase_struc)
# Selection
action = input("pretest or get_all: ").lower()
if action == "pretest":
    ase_struc_pre = decode(encode(list(out.values())[0]))
    os.chdir('/ocean/projects/mat210008p/jhuanga/project/%s/task_pretest'%project_name)
    write("geometry.in", ase_struc_pre, format="aims")

    dir_lo_pre = f'/ocean/projects/mat210008p/jhuanga/project/{project_name}/task_pretest/'
    os.chdir('..')
    os.system(f'cp -i submit.sh control.in {dir_lo_pre}')
    command = input('Submit the task now? (yes/no): ').lower()
    if command == 'yes':
        os.chdir(dir_lo_pre)
        os.system('sbatch submit.sh')
    elif command == 'no':
        print('Submit it later')

elif action == "get_all":
    os.chdir('/ocean/projects/mat210008p/jhuanga/project/%s/task_all'%project_name)
    for i, ase_struc in enumerate(ase_list):
        dir_name = 'structure%i' % i
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        dir_lo=f'/ocean/projects/mat210008p/jhuanga/project/{project_name}/task_all/{dir_name}/'
        os.chdir(f'/ocean/projects/mat210008p/jhuanga/project/{project_name}')
         # Copy the necessary files to the destination directory
        os.system(f'cp -i submit.sh control.in {dir_lo}')
        os.chdir(dir_lo)
        write("geometry.in", ase_struc, format="aims")
        os.chdir('..')  # Move back to the parent directory (task_all)

# Run the task
    command = input('Submit the task now? (yes/no): ').lower()
    if command == 'yes':
        for i, ase_struc in enumerate(ase_list):
            dir_name = 'structure%i' % i
            os.chdir(dir_lo)
            os.system('sbatch submit.sh')
    elif command == 'no':
        print('Submit it later')
    

        # End of the script
