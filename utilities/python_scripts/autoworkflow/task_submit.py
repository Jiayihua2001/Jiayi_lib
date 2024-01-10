# read json files
from ase.io.jsonio import encode, decode
from ase.io import write
import json
import os

def task_submission(struc_json,pretest=True):
    """
    pretest or run all tasks by inputing json file ,automatically submit the task.

    Arguments
    ---------
    struct json: str
        json_file path
    

    pretest :bool
        defalut= true
        after the pretest, if you want to run all the task, change the pretest=False
    
    """
    #prepare for the path
    tp=os.getcwd()
    os.mkdir('calc')
    os.mkdir('pretest')
    control_p = os.path.abspath('control.in')
    submit_p = os.path.abspath('submit.sh')
    calc_p = os.path.abspath('calc')
    pre_p = os.path.abspath('pretest')
    

    #get the struct
    with open(struc_json, "r") as f:
        out = json.load(f)

    ase_dict = {k:decode(encode(v)) for k,v in out.items()}

    # Selection pretest or not
    if pretest == True:
        ase_struc_pre =list(ase_dict.values())[0]
        os.chdir(pre_p)
        write("geometry.in", ase_struc_pre, format="aims")
        os.system('cp {control_p} .')
        os.system('cp {submit_p} .')
        os.system('sbatch submit.sh')
        os.chdir(tp)

    elif pretest == False:
        os.chdir(calc_p)
        for key,values in ase_dict.items():
            temp_path = os.path(calc_p,key)
            if not os.path.isdir(temp_path):
                os.mkdir(temp_path)
            os.chdir(temp_path)
            write("geometry.in", values, format="aims")
            os.system('cp {control_p} .')
            os.system('cp {submit_p} .')
            os.system('sbatch submit.sh')
            os.chdir(calc_p)


        # End of the script
if __name__ == "__main__":
    struc_json=''
    task_submission(struc_json,pretest=True)
    