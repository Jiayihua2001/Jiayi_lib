
import os
import numpy as np
import ibslib
from ibslib import Structure,SDS
from ibslib.io import read,write,check_dir
from ibslib.driver import BaseDriver_
from ibslib.molecules import check_molecule,FindMolecules,align
from ibslib.molecules.symmetry import get_symmetry
from ibslib.molecules import com,moit,orientation,show_axes,align
from ibslib.molecules.cutoff_matrix import CutoffMatrix
from ibslib.molecules.orientations import get_unique_angle_grid
from ibslib.crystals import align_structs,calc_rmsd,rot_struct
from ibslib.calculators.aims_driver import AimsDriver
from ibslib.dimers.cutoff_matrix import get_cutoff_matrix
from scipy.spatial.distance import cdist
from scipy.optimize import minimize
from ase import Atoms
from ase.io import read,write
from find_molecules import FindMolecules
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize
import cutoff_matrix



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

class DimerRelaxation(BaseDriver_):
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
    def __init__(self,cluster,int_scale,
                 offset_list=np.arange(-0.5,2+0.05,0.05), 
                 ):
        # initialize
        self.fm = FindMolecules()
        self.offset_list = offset_list
        self.cluster = cluster
        self.int_scale = int_scale

        self.com_list = []
        self.com_vector = []
        self.dimer_list = []

        # unpack cluster. Get mol properties
        self.process_cluster()



    def process_cluster(self):
        self.fm.calc_struct(cluster)
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
        energy[dist < radius] = np.inf
        energy[dist > self.D] = 0
        energy[(dist > radius) & (dist < self.D)] = (wt * (self.D - dist) / (dist - radius))[
            (dist > radius) & (dist < self.D)
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
        
    
   


if __name__ =='__main__':

    cluster = read('/ocean/projects/mat210008p/jhuanga/projects/target11/rigid_press/test.in')
    cluster = Structure.from_ase(cluster)
    dimer_rigid=DimerRelaxation(cluster,int_scale=0.1)
    pressed_struc=dimer_rigid.calc_struct(cluster)
    ibslib.io.write('pressed_geo.in',pressed_struc,file_format='geo')