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
        

def main(object='xtal',task='relaxation'):
        """
        object: str
                "xtal" or "cluster" (non-period system)
        task: str
                "energy" or "relaxation"

        """
        s = read('geometry.in')
        if object =='xtal':
                if task =="energy":         
                        write_control(s, light_path, energy_settings)
                elif task =="relaxation":
                        write_control(s, light_path, relaxed_settings)
                else:
                        print('incorrect input of task,please choose between "energy" and "relaxation "')
        else:
                if task =="energy":         
                        del(energy_settings['k_grid'])
                        write_control(s, light_path, energy_settings)
                elif task =="relaxation":
                        del(relaxed_settings['k_grid'])
                        del(relaxed_settings['relax_unit_cell'])
                        write_control(s, light_path, relaxed_settings)
                else:
                        print('incorrect input of task, please choose between "energy" and "relaxation "')

        return

if __name__ =="__main__":
        main()
