# automatically write control.in file 
from ibslib import read,write
import ase.io
import numpy as np 
import os
import json
from ase.io.jsonio import decode, encode


light_path = '/ocean/projects/mat210008p/shared/Software/FHIaims_2023_1/fhi-aims.221103_1/species_defaults/defaults_2020/light'

relaxed_settings = {'xc': 'pbe',
					'spin': 'none',
					'relativistic': 'none',
					'charge': 0,     #physical model

                                        #SCF convergence
					'occupation_type': 'gaussian 0.01',
					'mixer': 'pulay',
					'n_max_pulay': 8,
					'charge_mix_param': 0.2, 
					'sc_accuracy_rho': 0.0001, 
					'sc_accuracy_eev': 0.01,
					'sc_accuracy_etot': 1e-06, 
					'sc_iter_limit': 10000, 
                                        'check_cpu_consistency': '.true.',
                                        #  Eigenvalue solution
					'KS_method': 'parallel', 
					'empty_states': 6, 
					'basis_threshold': 1e-05,
                                        # For periodic boundary conditions
					'k_grid': '3 3 3', 
                                        #  Relaxation   if just calculate energy ,please uncommon 
					'sc_accuracy_forces': 0.0001, 
					'relax_geometry': 'trm 1e-2', 
					'relax_unit_cell': 'full',
					 'hessian_to_restart_geometry': '.false.', 
                                         # relax_unit_cell fixed_angles #uncomment for unit cell relaxation with fixed angles
					 # restart_relaxations .true.
                                         'harmonic_length_scale': 0.01, #uncomment in case of relaxation errors
					 'energy_tolerance': 0.0005,  #uncomment in case of relaxation errors
					 'many_body_dispersion': ''}

energy_settings = {'xc': 'pbe',
					'spin': 'none',
					'relativistic': 'none',
					'charge': 0,     #physical model

                                        #SCF convergence
					'occupation_type': 'gaussian 0.01',
					'mixer': 'pulay',
					'n_max_pulay': 8,
					'charge_mix_param': 0.2, 
					'sc_accuracy_rho': 0.0001, 
					'sc_accuracy_eev': 0.01,
					'sc_accuracy_etot': 1e-06, 
					'sc_iter_limit': 10000, 
                                        'check_cpu_consistency': '.true.',
                                        #  Eigenvalue solution
					'KS_method': 'parallel', 
					'empty_states': 6, 
					'basis_threshold': 1e-05,
                                        # For periodic boundary conditions
					'k_grid': '3 3 3', 
                                         # relax_unit_cell fixed_angles #uncomment for unit cell relaxation with fixed angles
					 # restart_relaxations .true.
                                         'harmonic_length_scale': 0.01, #uncomment in case of relaxation errors
					 'energy_tolerance': 0.0005,  #uncomment in case of relaxation errors
					 'many_body_dispersion': ''}
#  Vdw Corrections  choose vch or mbd 
# del(relaxed_settings['many_body_dispersion'])
# relaxed_settings['vdw_correction_hirshfeld'] = ''   

def unique_preserve_order(array):
        _, idx = np.unique(array, return_index=True)
        array = array[np.sort(idx)]
        return array

def write_control(struct, species_path, settings):
        #get the unique elements in the structure
        symbols = np.array(struct.get_chemical_symbols())
        print('symbols:',symbols)
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
        print('element:',ase_element)
        element = unique_preserve_order(element)
        control_file = 'control.in'
        with open(control_file, 'w') as f:
                for setting, value in settings.items():
                        f.write('{}   {}\n'.format(setting,value))
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




# struc = read('/ocean/projects/mat210008p/jhuanga/projects/target11/unique_dimer/test.json',file_format='json')
# ase_struc = struc.get_ase_atoms()   #the same 
# write('geometry',struc,file_format='geo')

struc=read('mol.in',file_format= 'geo')
# for relaxation
write_control(struc, light_path, relaxed_settings)

# for energy
#write_control(struc, light_path, energy_settings)