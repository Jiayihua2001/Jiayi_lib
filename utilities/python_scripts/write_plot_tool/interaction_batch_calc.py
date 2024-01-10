from ase.io import read, write
from copy import deepcopy
import math
import numpy as np 
import os

inter_path = '/ocean/projects/mat210008p/shared/Software/FHIaims_2023_1/fhi-aims.221103_1/species_defaults/defaults_2020/light'
dft = ['ts', 'mbd', 'pbe0']

arjuna_arguments = {'-J':'CSP', '-N': 1, '-n': 16, '-o': 'j_%j.out', '-p': 'RM', '-A': 'mat210008p', '-t':'48:00:00',
                                         'pre-command': ['ulimit -s unlimited', 'module load anaconda', 'module load openmpi/4.0.2-intel20.4', 'conda activate gnrs_env'],
                                          'command': 'mpirun -np 16 /ocean/projects/mat210008p/shared/bin/aims.210716_1.scalapack.mpi.x > aims.out'}

SPE_settings = {'xc': 'pbe',
                                        'spin': 'none',
                                        'relativistic': 'atomic_zora scalar',
                                        'charge': 0,
                                        'occupation_type': 'gaussian 0.01',
                                        'mixer': 'pulay',
                                        'n_max_pulay': 8,
                                        'charge_mix_param': 0.02, 
                                        'sc_accuracy_rho': 0.0001, 
                                        'sc_accuracy_eev': 0.01,
                                        'sc_accuracy_etot': 1e-06, 
                                        'sc_iter_limit': 10000, 
                                        'KS_method': 'parallel', 
                                        'empty_states': 6, 
                                        'basis_threshold': 1e-05,
                                        'k_grid': '3 3 3', 
                                        'sc_accuracy_forces': 0.0001,  
                                         'many_body_dispersion': ''}

def unique_preserve_order(array):
        _, idx = np.unique(array, return_index=True)
        array = array[np.sort(idx)]
        return array

def write_control(struct, species_path, settings):
        #get the unique elements in the structure
        symbols = np.array(struct.get_chemical_symbols())
        symbols = unique_preserve_order(symbols)
        #get the corresponding atomic numbers
        ase_element = struct.get_atomic_numbers()
        element = []
        for ele in ase_element:
                if ele < 10:
                        new_ele = str(0) + str(ele)
                else:
                        new_ele = str(ele)
                element.append(new_ele)
        element = np.array(element)
        element = unique_preserve_order(element)
        control_file = 'control.in'
        with open(control_file, 'w') as f:
                for setting, value in settings.items():
                        f.write('{}    {} \n'.format(setting, value))
        f.close()
        #create a list (species_list) of "species paths" from the elements and symbols in the array
        element = list(element)
        symbols = list(symbols)
        #TODO: fix this so that the elements correspond with the right symbols, the list is getting sorted somewhere... maybe in unique call
        species_list = []
        for ele, symbol in zip(element, symbols):
                species_list.append('{}_{}_default'.format(ele, symbol))
        species_list = [os.path.join(species_path, specie) for specie in species_list]
        for specie in species_list:
                with open(specie, 'r') as first_file, open(control_file, 'a') as second_file:
                        for line in first_file:
                                second_file.write(line)

def k_grid_25(struct):
        #can pass in the pymatgen lattice vectors and get a clear representation of the lattice vectors
        lv = list(struct.cell.lengths())
        k_grid = [math.ceil(25/vector) for vector in lv]
        k_grid = '{} {} {}'.format(k_grid[0], k_grid[1], k_grid[2])
        return k_grid

def write_batch_script(filename, arguments):
        with open(filename, 'w') as f:
                f.write('#!/bin/sh\n')
                for arg, key in arguments.items():
                        if arg != 'pre-command' and arg != 'command' and arg != '--mem-per-cpu':
                                f.write('#SBATCH {} {} \n'.format(arg, key))
                        elif arg == '--mem-per-cpu':
                                f.write('#SBATCH {}={}\n'.format(arg, key))
                        elif arg == 'pre-command':
                                for command in key:
                                        f.write(command + '\n')
                        elif arg == 'command':
                                f.write(key)
        f.close()


cwd = os.getcwd()

s = read('geometry.in')
SPE_settings['k_grid'] = k_grid_25(s)

for method in dft:
    method_path = os.path.join(cwd, method)
    os.mkdir(method_path)
    os.chdir(method_path)
    if method == 'ts':
        TS_settings = deepcopy(SPE_settings)
        del(TS_settings['many_body_dispersion'])
        TS_settings['vdw_correction_hirshfeld'] = ''
        write_control(s, inter_path, TS_settings)
    elif method == 'mbd':
        MBD_settings = SPE_settings
        print(SPE_settings)
        write_control(s, inter_path, MBD_settings)
    elif method == 'pbe0':
        PBE0_settings = deepcopy(SPE_settings)
        PBE0_settings['xc'] = 'pbe0'
        write_control(s, inter_path, PBE0_settings)
    write_batch_script('submit.sh', arjuna_arguments)
    write('geometry.in', s, format='aims')
    os.chdir(cwd)