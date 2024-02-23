#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import time
import os
import numpy as np
import ibslib
from ibslib import Structure
from ibslib.io import read,write
from ase.io.jsonio import encode,decode
from ibslib.driver import BaseDriver_
from ibslib.molecules import FindMolecules
from ibslib.molecules import com
from scipy.spatial.distance import cdist
from ase import Atoms
from find_molecules import FindMolecules
from scipy.spatial.distance import cdist
import cutoff_matrix
from mpi4py import MPI

def read_json(jsonfile_dir):  #read -to lib struc
    """
    Load Json file to return a dictionary storing ase structures
    """
    struc_dict = read(jsonfile_dir,file_format='json')
    return struc_dict

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

def scatter_list(input_list, num_cores):
    """
    Slicing the data based on the number of cores.
    """
    ave, res = divmod(len(input_list), num_cores)
    counts = [ave + 1 if p < res else ave for p in range(num_cores)]
    # determine the starting and ending indices of each sub-task
    starts = [sum(counts[:p]) for p in range(num_cores)]
    ends = [sum(counts[: p + 1]) for p in range(num_cores)]
    # save the starting and ending indices in data
    output_list = [input_list[starts[p] : ends[p]] for p in range(num_cores)]
    return output_list

def find_molecule(struc):
    find = FindMolecules(output_rstruct=True)
    unique = find.calc(struc)
    return unique

def standardize_mol(mol):
    """
    Move the center of geometry to origin and align
    the molecule with the principal axes of rotation.
    """
    
    an = mol.get_atomic_numbers()
    mol.set_atomic_numbers([1 for _ in range(len(mol))])
    com = mol.get_center_of_mass()
    pos = mol.get_positions()
    pos -= com
    mol.set_positions(pos)
    mol.set_atomic_numbers(an)
    ref_mol = mol
    return ref_mol


def get_mol_len(ref_mol):
    """
    Returns the diameter of the smallest sphere that can
    enclose a molecule.
    """

    pos = ref_mol.positions
    cog = np.mean(pos, axis=0)
    dist = cdist(pos, [cog])
    return 2 * dist.max()


def is_duplicate(cog1, cog2):
    """
    Check if the difference in a close to an integer.
    """

    intval = np.round(cog1 - cog2)
    is_dup = np.allclose(intval, cog1 - cog2, atol=0.001, rtol=0)
    return is_dup

class Dimer_RigidPress(BaseDriver_):
    """
    Finds the optimal location for the COM distances between dimers. This is
    done by increasing the distance between the dimers along the vector of their
    center of mass. 
    
    Arguments
    ---------
    fm: ibslib.molecules.FindMolecules
        Driver that will obtain molecules from the dimer structure
    offset_list: iterable
        List of com_vector distances to evaluate for grid based minimization. 
        If the length of this list is zero, the minimization using scipy
        will be performed. 
    minimize_kw: dictionry
        Arguments for minimizer that follows the scipy.optimize.minimize API. 
        Optional input as a standard arguments for dimer relaxations will be 
        initialized if none are provided. 
    
    """
    def __init__(self,int_scale,sr,
                 offset_list=np.arange(-0.5,2+0.05,0.05), 
                 ):
        # initialize
        self.fm = FindMolecules()
        self.offset_list = offset_list
        self.int_scale = int_scale
        self.sr=sr
        self.com_list = []
        self.com_vector = []
        self.dimer_list = []

        # unpack cluster. Get mol properties



    def process_cluster(self,dimer):
        self.fm.calc_struct(dimer)
        mol1=self.fm.molecules[0]
        mol2=self.fm.molecules[1]
        cutoff_m = cutoff_matrix.get_cutoff_matrix(mol1,mol2)
        self.mol1=mol1.get_ase_atoms()
        self.mol2=mol2.get_ase_atoms()
        
        self.ref_mol = standardize_mol(self.mol1)
        mol_length = get_mol_len(self.ref_mol)
        self.natoms = len(self.ref_mol)
        self.radius = cutoff_m[0 : self.natoms, 0 : self.natoms].ravel()
        # Interaction distance = mol length + 2x max cutoff
        self.D = mol_length + 2 * max(self.radius)

    

    def interaction_kernel(self, dist, radius):
        """
        Defines how atoms interact based on distance.

        Returns the total energy of interaction.
        """
        wt = self.int_scale / (self.natoms * self.natoms)
        energy = np.zeros(len(dist))
        energy[dist < self.sr*radius] = np.inf
        energy[dist > self.D] = 0
        energy[(dist > self.sr*radius) & (dist < self.D)] = (wt * (self.D - dist) / (dist - self.sr*radius))[
            (dist > self.sr*radius) & (dist < self.D)
        ]
        return energy.sum()

    
    def total_energy(self,dimer):
        """
        Computes the total energy of the cluster based on positions.

        Parameters:
        - positions: A flat array of positions that will be reshaped into (N, 3).

        Returns the total energy.
        """

        # Reshape positions to (N, 3) for mol1 and mol2
        self.process_cluster(dimer)

        mol1_positions = self.mol1.positions.reshape((-1,3))
        mol2_positions = self.mol2.positions.reshape((-1,3))
        natoms_mol1 = np.shape(mol1_positions)[0]
        dist = cdist(mol1_positions, mol2_positions).T.reshape(-1)
        stretched_radius = np.resize(self.radius, dist.shape)
        return self.interaction_kernel(dist, stretched_radius)

    def calc_struct(self, struct):
        """
            Main route ,only offset method can run
        """
        self.fm.calc_struct(struct)
        self.com_list = [com(x) for x in self.fm.molecules]
        self.com_vector = self.com_list[1] - self.com_list[0]
        
        if len(self.offset_list) > 0:
            return self._calc_struct_offset_range(struct)
        else:
            return self._calc_struct_minimize(struct)
    
    def generate_offset(self, molecules, offset):

        if len(self.com_list) == 0:
            self.com_list = [com(x) for x in molecules]
            self.com_vector = self.com_list[1] - self.com_list[0]
        
        geo0 = molecules[0].get_geo_array()
        geo1 = molecules[1].get_geo_array()
        ele0 = molecules[0].geometry["element"]
        ele1 = molecules[1].geometry["element"]
        
        dot_value = offset / np.linalg.norm(self.com_vector)
        offset_vector = self.com_vector*dot_value
        offset_total_vector = self.com_vector + offset_vector

        ### Fraction by which the com_vector must be multiplied to get the 
        ### given offset
        # offset_ratio = offset / np.linalg.norm(self.com_vector)
        # offset_vector = self.com_vector*dot_value
        # offset_total_vector = self.com_vector + offset_vector
        
        diff = np.linalg.norm(offset_total_vector - self.com_vector)
        if np.abs(np.abs(offset) - diff) > 1e-4:
            raise Exception("New COM not calculated correctly. {} {}"
                    .format(np.abs(offset - diff), 
                            offset,
                            ))
        
        geo1 = geo1 + offset_vector
        
        dimer = Structure.from_geo(np.vstack([geo0, geo1]),
                                   np.hstack([ele0, ele1]))
        
        #### Visualization purposes
        # dimer.append(self.com_list[0][0],
        #              self.com_list[0][1],
        #              self.com_list[0][2],
        #              "Cl")
        
        # dimer.append(self.com_list[1][0] + offset_vector[0],
        #              self.com_list[1][1] + offset_vector[1],
        #              self.com_list[1][2] + offset_vector[2],
        #              "Br")
        
        dimer.properties["offset"] = offset
        dimer.properties["original_com_vector"] = list(self.com_vector)
        dimer.properties["offset_vector"] = list(offset_vector)
        dimer.properties["com_vector"] = list(offset_total_vector)
        dimer.properties["energy"] = 0
        dimer.struct_id = "Relaxation_{:e}".format(offset)
        
        return dimer
    
    def _calc_struct_offset_range(self, struct):
        """
        Calculate the relaxed dimer structure by calculating over the 
        given offset range provided during initialization of the class. 
        """

        ##### Generate all dimers based on offset values
        self.dimer_dict = {}
        for entry in self.offset_list:
            
            temp_dimer = self.generate_offset(self.fm.molecules, entry) 
            self.dimer_dict[temp_dimer.struct_id] = temp_dimer
        
        ##### Calculate all dimer energies

        ##### Store all results
        offset_list = []
        energy_list = []
        dimer_keys = []
        for temp_struct_id,temp_struct in self.dimer_dict.items():
            offset = temp_struct.properties["offset"]
            energy=self.total_energy(temp_struct)
            offset_list.append(offset)
            energy_list.append(energy)
            dimer_keys.append(temp_struct_id)
        
        sort_idx = np.argsort(offset_list)
        offset_list_sorted = [offset_list[x] for x in sort_idx]
        energy_list_sorted = [energy_list[x] for x in sort_idx]
        struct.properties["dimer_rigidpress"] = [x for x in 
                    zip(offset_list_sorted, energy_list_sorted)]
        
        ##### Set geometry as minimum energy dimer
        geometry_idx = np.argmin(energy_list)
        min_dimer_key = dimer_keys[geometry_idx]
        min_dimer = self.dimer_dict[min_dimer_key]
        struct.geometry = min_dimer.geometry
        for key,value in min_dimer.properties.items():
            struct.properties[key] = value

        return struct

def get_pressed_dic(struc_dict, id_list):
    """
    Runs pymatgen duplicate checks on a distributed struct dictionary.
    target : target structure to be matched
    Returns: List of matching structure names/ids
    """
    if len(id_list) != 0:
        pressed_dic = {}
        dimer_rigid=Dimer_RigidPress(sr=0.80,int_scale=0.1)   #exp -0.82/0.84
        for key in id_list:
            dimer = struc_dict[key]
            pressed_dimer=dimer_rigid.calc_struct(dimer)
            pressed_dic[key]=pressed_dimer.get_ase_atoms()

    else:
        pressed_dic = {}

    return pressed_dic

def main():
    # Initialise MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

     # Only the root process reads the entire JSON
    dimer_dict = None
    if rank == 0:
        start_time = time.perf_counter()
        json_dir = '/ocean/projects/mat210008p/jhuanga/projects/target11/dimer_generation_0.83_0.85'
        dimer_dict = read_json(json_dir)
        struct_id_list = list(dimer_dict.keys())
        struct_id_list = scatter_list(struct_id_list, size)
    else:
        struct_id_list = None

    dimer_dict = comm.bcast(dimer_dict,root=0)
    # dispatches struct_id_list from main process across all processes
    struct_id_list = comm.scatter(struct_id_list, root=0)

    # tasks for each rank
    pressed_dic = get_pressed_dic(dimer_dict,struct_id_list)
    #gather to  rank 0
    all_pressed_dic = comm.gather(pressed_dic, root=0)
    
    if rank == 0:
        # flatten
        pressed_dic = {k:v for d in all_pressed_dic for k,v in d.items()}
        write_json(pressed_dic,'rigid_press.json')
        end_time  = time.perf_counter()
        elapsed_time = end_time - start_time
        print(f'Time used is: {elapsed_time}')
    

    return
   


if __name__ =='__main__':
    main()

