# automatically write control.in file 

from ase.io import read, write
import numpy as np 
import os

light_path = '/ocean/projects/mat210008p/shared/Software/FHIaims_2023_1/fhi-aims.221103_1/species_defaults/defaults_2020/light'

relaxed_settings = {'xc': 'pbe',
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
					'relax_geometry': 'trm 1e-2', 
					'relax_unit_cell': 'full',
					 'hessian_to_restart_geometry': '.false.', 
					 'harmonic_length_scale': 0.01, 
					 'energy_tolerance': 0.0005, 
					 'many_body_dispersion': ''}

del(relaxed_settings['many_body_dispersion'])
relaxed_settings['vdw_correction_hirshfeld'] = ''

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

s = read('/ocean/projects/mat210008p/jhuanga/python_scripts/geometry.in')
write_control(s, light_path, relaxed_settings)