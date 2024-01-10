
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

#func of extracting info
def struc_relax(temp_path):
        return ase_struc_relaxed
    
def main(json_path):
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
        temp_path=os.path(calc_p,geometry_id)
        os.chdir(temp_path)
        aims_path=os.path.abspath('aims.out')
        if is_converge(aims_path):
            energy_dic[geometry_id]=find_energy(aims_path)
            if os.path.exists('geometry.in.next_step'):
                os.system("mv geometry.in.next_step geometry.in")
                ase_struc_relaxed = read("geometry.in")
                relax_dic[geometry_id]=json.loads(encode(ase_struc_relaxed))
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

    os.chdir(tp)
    with open("Total_energy.json",'w',encoding='utf8') as fe:
        json.dump(energy_dic,fe,indent=6) 
    with open("sructure_relaxed.json",'w',encoding='utf8') as fr:
        json.dump(relax_dic,fr,indent=6) 
    with open("unconverged_dic",'w',encoding='utf8') as fe:
        json.dump(unconverge_dic,fe,indent=6) 

    return num_unconverge

tp=os.getcwd()
calc_p=os.path.abspath('calc')
json_path=''
main(json_path)
