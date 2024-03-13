
from ibslib.molecules.symmetry import get_symmetry
import numpy as np
from ibslib.molecules import com,orientation
from ibslib.crystals import calc_rmsd
from mpi4py import MPI
from ibslib.libmpi.dup_check import DupCheck
from ibslib.molecules import FindMolecules
from ibslib import Structure
from ibslib.molecules.orientations import get_unique_angle_grid
import json
from ibslib.io import read,write
from ase.io.jsonio import encode,decode

def center_dimer(struct, choose_min=True, fm=None):
    """
    Centers the system for one the molecules in the dimer such that the COM is 
    at the origin and the axes defined by the moment of inertia tensor are  
    oriented with the origin. This is done by translating and rotation the 
    entire dimer system.
    
    Arguments
    ---------
    struct: Structure
        Structure to adjust
    choose_min: bool
        If True, will center on molecule closest to the origin.
        If False, will center on molecule furtherst from origin. 
    fm: FindMolecules
        FindMolecules object that has been calculated using Structure. This will
        save time. Otherwise, default will calculate the molecules in the 
        Structure. 
    
    """
    if fm == None:
        fm = FindMolecules()
        fm.calc(struct)
        
    # check_dimer(struct, fm=fm)
    
    com_list = np.vstack([com(x) for x in fm.molecules])
    com_list = np.linalg.norm(com_list, axis=-1)
    
    if choose_min:
        chosen_idx = np.argmin(com_list)
        other_idx = np.argmax(com_list)
    else:
        chosen_idx = np.argmax(com_list)
        other_idx = np.argmin(com_list)
        
    trans = com(fm.molecules[chosen_idx])
    rot = orientation(fm.molecules[chosen_idx])
    
    geo = struct.get_geo_array()
    geo = geo - trans
    geo = np.dot(geo, rot.T)
    struct.from_geo_array(geo, struct.geometry["element"])
    
#    ## fm should also be modified
#    align(fm.unique_molecule_list[0])
#    for idx,molecule in enumerate(fm.molecule_struct_list):
#        geo = molecule.get_geo_array()
#        geo = geo - trans
#        geo = np.dot(geo, rot.T)
#        molecule.from_geo_array(geo, molecule.geometry["element"])    
#        fm.molecule_struct_list[idx] = molecule


def align_dimers(struct1, struct2, mol_rot=[], tol=0.1):
    """
    Attempts to align dimers using the symmetry given. 
        
    """
    if len(mol_rot) == 0:
        raise Exception()
    
    ### For compatability layer with libmpi.DupCheck
    dimer1 = struct1
    dimer2 = struct2
    
    geo1 = dimer1.get_geo_array()
    geo2 = dimer2.get_geo_array()
    ele = dimer1.geometry["element"]
    
    for rot in mol_rot:
        # Rot = R.from_matrix(rot)
        # temp_geo = Rot.apply(geo2)
        temp_geo = np.dot(geo2,rot)
        
        temp_orientation = Structure.from_geo(temp_geo, ele)
        temp_rmsd = calc_rmsd(dimer1, temp_orientation)
        
        # temp_rmsd = np.mean(np.linalg.norm(geo1 - temp_geo,axis=-1))
        if temp_rmsd < tol:
            return True

    #### It's also an allowed symmetry operation that the indices of the dimer 
    #### switch such that the other dimer is centered at the origin.
        
    # center_dimer(dimer2, choose_min=False)
    # geo1 = dimer1.get_geo_array()
    # geo2 = dimer2.get_geo_array()
    # for rot in mol_rot:
    #     # Rot = R.from_matrix(rot)
    #     # temp_geo = Rot.apply(geo2)
    #     temp_geo = np.dot(geo2,rot)
        
    #     temp_orientation = Structure.from_geo(temp_geo, ele)
    #     temp_rmsd = calc_rmsd(dimer1, temp_orientation)
        
    #     ##### Shouldn't have to reorder indices at all because changing the 
    #     ##### centered dimer is explicitly considered as a symmetry operation. 
    #     # combined = align_structs(dimer1,temp_orientation)
    #     # print("NAIVE,COMBINED: {}, {}".format(temp_rmsd,combined.properties["RMSD"]))
    #     # if combined.properties["RMSD"] < 0.7:
    #     #     test_stream.update(combined)
        
    #     # temp_rmsd = np.mean(np.linalg.norm(geo1 - temp_geo,axis=-1))
    #     if temp_rmsd < tol:
    #         return True
    
    return False

def remove_symmetric_orientations(self, mol):
    """
    New implementation utilizes the molecule.orientations.get_unique_angle_grid
    such that finding orientations takes advantage of molecular symmetry
    without having to explicity compute RMSD as in previous method. 
    
    """
    ## Build orientation dict on rank 0. 
    ## Will be communicated to other ranks in DupCheck
    if self.rank == 0:
        return get_unique_angle_grid(
                mol, 
                angle_spacing=self.angle_spacing, 
                max_angle=360,
                max_rot=10,
                tol=self.tol)
    else:
        return []


def decide_mol(mol,mol_symmetry=True):
    ### Choose the central molecule based on which choice will produce 
    ### the smaller number of unique angles
    fm = FindMolecules()
    fm.calc(mol)
    mol1=fm.molecules[0]
    mol2=fm.molecules[1]

    if mol_symmetry:
        angle_list_1 = remove_symmetric_orientations(mol1)
        angle_list_2 = remove_symmetric_orientations(mol2)
        if len(angle_list_1) > len(angle_list_2):
            mol1 = mol1
            mol2 = mol2
            angle_list = angle_list_2
        else:
            mol1 = mol2
            mol2 = mol1
            angle_list = angle_list_1
    else:
        ### Otherwise just have the molecule with the smaller number
        ### of atoms at the center           
        if len(mol1.geometry) > len(mol2.geometry):
            mol1 = mol1
            mol2 = mol2
        else:
            mol1 = mol2
            mol2 = mol1
    return mol1,mol2

def read_json(jsonfile):
    """
    Load Json file to return a dictionary storing ase structures
    """

    with open(jsonfile, "r") as f:
        stru_dict = json.load(f)

    ase_dict = {k: decode(encode(v)) for k, v in stru_dict.items()}
    return ase_dict

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

def main():

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        mol=read('mol.in')
        mol1,mol2 = decide_mol(mol)
        # mol_rot should be the mol that is actually going to be rotate
        mol_rot = get_symmetry(mol2)
        #load json as lib struc
        dimer_dict = read('./dimer_generation') 
        #### Performing duplicate check
    

    mol_rot=comm.bcast(mol_rot,root=0)
    dimer_dict=comm.bcast(dimer_dict,root=0)
    dc = DupCheck(comm=comm)
    unique_dict = dc.calc(align_dimers, dimer_dict,
                            compare_fn_kw={"mol_rot": mol_rot})
    if rank==0:
        tol=0.5  # change the tol in align_dimers
        ase_unique_dict ={k:v.get_ase_atoms for k,v in unique_dict.items()}
        write_json(ase_unique_dict,f'dup {tol} _generation.json')
    
    
   
    