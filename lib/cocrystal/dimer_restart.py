
import os 
from ase.io import read,write
import json




def dimer_setting(user_settings,num_1,num_2):

    out_dir=user_settings
    mol_path_l=out_dir["rtm_set"]["molecule_path"]
    geo_p=mol_path_l[0]
    head, tail = os.path.split(geo_p)
    mol1_p=os.path.join(head,'mol1.in')
    mol2_p=os.path.join(head,'mol2.in')
    mol_path_l.remove(geo_p)
    mol_path_l.append(mol1_p)
    mol_path_l.append(mol2_p)

    out_dir["rtm_set"]["molecule_path"]=mol_path_l
    n_mol_l=out_dir["rtm_set"]["n_atoms_in_mol"]
    n_dimer=n_mol_l[0]
    try:
        if num_1 + num_2 == n_dimer:
            n_mol_l.remove(n_dimer)
            n_mol_l.append(num_1)
            n_mol_l.append(num_2)
    except:
        print(f'error: total n of atoms {n_dimer} != that sum of mol1 and mol2: {num_1}+{num_2}')
    out_dir["rtm_set"]["n_atoms_in_mol"]=n_mol_l
    
    out_dir["usr_set"]["generation"]["generation_type"]='cocrystal'
    out_dir["usr_set"]["generation"]["stoichiometry"]= [1,1]
    out_dir["usr_set"]["master"]["molecule_path"] = ["./mol1.in","./mol2.in"]

    return out_dir

def main():
    os.system('cp ./mol1.in ./tmp/molecule')
    os.system('cp ./mol2.in ./tmp/molecule')
    mol_1=read('./mol1.in')
    mol_2=read('./mol2.in')
    os.system('rm -f ./tmp/molecule/geometry_0.in')
    restart_path = os.path.join('./tmp', "restart.json")
    with open(restart_path,'r') as jf:
        user_settings=json.load(jf)
    num_1=len(mol_1)
    num_2=len(mol_2)
    new_user_settings=dimer_setting(user_settings,num_1,num_2)
    os.system('rm -f ./tmp/restart.json')
    with open(restart_path,'w') as f:
        json.dump(new_user_settings,f,indent=6)
    os.system('sbatch submit.sh')


if __name__ =="__main__":
    cp=os.getcwd()
    for i in os.listdir(cp):
        try:
            tmp=os.path.join(cp,i)
            os.chdir(tmp)
            os.chdir('40000')
            os.chdir('relax')
            os.chdir('rigid')
            main()
            os.chdir(cp)
        except:
            print('attempt fail')

