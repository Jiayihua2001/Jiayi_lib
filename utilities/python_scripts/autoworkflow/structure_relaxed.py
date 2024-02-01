
#!/usr/bin/python

import os
import json
from ase.io import read, write
from ase.io.jsonio import encode, decode

def is_converge(aims_path):
    """
    check if converged
    """
    converge = False
    for line in reversed(list(open(aims_path))):
        if 'Have a nice day.' in line:  
            converge = True
            break
    return converge
            
def find_energy(aims_path):
    """
    Get the total energy
    """
    for line in reversed(list(open(aims_path))):
        if "| Total energy of the DFT / Hartree-Fock s.c.f. calculation      :" in line:
            energy = line.split(":")[1]
            energy = energy.split()[0]
            break
    return float(energy)

def write_json(out_dict, filename):
    """
    Write a Json file 
    out_dic should be ase_dic
    """

    str_list = [f'"{k}": {encode(v)},\n' for k, v in out_dict.items()]
    str_list[-1] = str_list[-1][:-2]
    with open(filename, "w") as wfile:
        wfile.write("{\n")
        wfile.writelines(str_list)
        wfile.write("\n}")

    return
    
def main(json_path,relaxed= True):
    """
    post-processing 
    """
    with open(json_path,'r') as f:
        struct_dict=json.load(f)

    energy_dic={}
    relax_dic={}
    num_unconverge=0
    unconverge_dic={}

    for geometry_id,values in struct_dict.items():
        temp_path=os.path.join(calc_p,geometry_id)
        os.chdir(temp_path)
        aims_path=os.path.abspath('aims.out')
        if is_converge(aims_path):
            energy_dic[geometry_id]=find_energy(aims_path)
            if relaxed:
                ase_struc_relaxed = read("geometry.in")
                relax_dic[geometry_id]=ase_struc_relaxed
                os.chdir(calc_p)
            else:
                os.chdir(calc_p)
           

        else:
            if os.path.exists('geometry.in.next_step'):
                os.system("mv geometry.in.next_step geometry.in")
                os.system('sbatch submit.sh')
            else:
                pass
            num_unconverge += 1 
            unconverge_dic[geometry_id]=values
            os.chdir(calc_p)
            continue
    
    sorted_energy_dic = {k: v for k, v in sorted(energy_dic.items(), key=lambda item: item[1])}
    sorted_relax_dic={}
    for k,v in sorted_energy_dic.items():
        if k in relax_dic:
            sorted_relax_dic[v]=relax_dic[k]
    os.chdir(tp)
    write_json(sorted_relax_dic,'sructure_relaxed.json')
    write_json(sorted_energy_dic,'sructure_energy.json')
    
    with open("unconverged_dic",'w',encoding='utf8') as fe:
        json.dump(unconverge_dic,fe,indent=6) 

    return num_unconverge

tp=os.getcwd()
calc_p=os.path.abspath('calc')
json_path='/ocean/projects/mat210008p/jhuanga/projects/target11/ase_aims/relaxation_dup/structures_removed.json'
main(json_path)
