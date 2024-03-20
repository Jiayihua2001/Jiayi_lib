
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

    return False


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
        dimer_dict = read('./dimer_generation') 
        #### Performing duplicate check
    mol_rot=np.eye(3)  #identity matrix if there is no symmetry
    mol_rot=comm.bcast(mol_rot,root=0)
    dimer_dict=comm.bcast(dimer_dict,root=0)
    dc = DupCheck(comm=comm)
    unique_dict = dc.calc(align_dimers, dimer_dict,
                            compare_fn_kw={"mol_rot": mol_rot})
    if rank==0:
        tol=0.5  # change the tol in align_dimers
        ase_unique_dict ={k:v.get_ase_atoms() for k,v in unique_dict.items()}
        write_json(ase_unique_dict,f'dup {tol} _generation.json')
    
    
   
    