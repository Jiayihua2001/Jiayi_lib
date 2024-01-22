import os
import ibslib
from ibslib.io import read, write
import numpy as np

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

def task_submit(pp):
    tp=os.path.join(pp,'test_dir')
    os.chdir(tp)
    temp_list=[os.path.join(tp,f) for f in os.listdir(tp)]
    for temp_path in temp_list:
        os.chdir(temp_path)
        aims_path=os.path.join(temp_path,'aims.out')
        os.system(f'cp /ocean/projects/mat210008p/jhuanga/projects/target11/submit.sh {temp_path}')
        os.remove('aims.out')
        os.system('sbatch submit.sh')
        os.chdir(tp) 
    os.chdir(pp)
    return 

def sort_energy(pp):
    tp=os.path.join(pp,'test_dir')
    temp_list=[os.path.join(tp,f) for f in os.listdir(tp)]
    energy_dic={} 
    os.chdir(tp)
    for temp_path in temp_list:  
        os.chdir(temp_path)
        aims_path=os.path.join(temp_path,'aims.out')
        if  is_converge(aims_path):
            energy=find_energy(aims_path)
            energy_dic[temp_path]=energy
        else:
            print(f'task not convergent:{temp_path}')               
    os.chdir(tp) 
    #sort energy
    min_value = min(energy_dic.values())  # Find the minimum value once
    # List comprehension to get all key-value pairs with the minimum value
    min_temp_path = [(key, value) for key, value in energy_dic.items() if value == min_value]
    print(f'minimum energy:{min_temp_path[0][1]}')
    os.chdir(min_temp_path[0][0])
    dimer_struct = read('geometry.in',file_format='geo')
    os.chdir(tp) 
    return dimer_struct
    

pp = os.getcwd()
final_path=os.path.join(pp,'final_dimer')
final_dimer = sort_energy(pp)
write(final_path,final_dimer,file_format='geo')
