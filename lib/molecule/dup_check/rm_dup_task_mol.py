#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import gc
import os
import json
import time
import numpy as np
import copy
from mpi4py import MPI
from ase.io import read
from ase.io.jsonio import encode, decode
from ibslib.structure import Structure
from scipy.optimize import linear_sum_assignment
from scipy.spatial.transform import Rotation 
from scipy.spatial.distance import cdist
from pymatgen.analysis.molecule_matcher import MoleculeMatcher, IsomorphismMolAtomMapper
from ibslib.molecules.align import fast_align
from ibslib.molecules.utils import rot_mol

def calc_rmsd(struct1,struct2):
    temp1 = copy.deepcopy(struct1)
    temp2 = copy.deepcopy(struct2)
    ## Align the inertial tensor  and COM with origin
    temp1 = fast_align(temp1)
    temp2 = fast_align(temp2)

    geo1 = temp1.get_geo_array()
    geo2 = temp2.get_geo_array()
    ele1 = temp1.geometry["element"]
    ele2 = temp2.geometry["element"]

    dist = cdist(geo1,geo2)
    idx1,idx2 = linear_sum_assignment(dist)
    
    geo1 = geo1[idx1,:]
    geo2 = geo2[idx2,:]
    ele1 = ele1[idx1]
    ele2 = ele2[idx2]

    Rot,rmsd = Rotation.align_vectors(geo1,geo2)
    rot = Rot.as_matrix()

    geo2 = Rot.apply(geo2)
    temp2.from_geo_array(geo2, ele2)

    # temp2=rot_mol(rot,geo2)

    geo1 = temp1.get_geo_array()
    geo2 = temp2.get_geo_array()
    
    dist = cdist(geo1,geo2)
    
    idx1,idx2 = linear_sum_assignment(dist)
    
    geo1 = geo1[idx1]
    geo2 = geo2[idx2]
    ele1 = ele1[idx1]
    ele2 = ele2[idx2]

    result = np.mean(np.linalg.norm(geo1 - geo2,axis=-1))
    return result


def read_json(jsonfile):
    """
    Load Json file to return a dictionary storing ase structures
    """

    with open(jsonfile, 'r') as f:
        out = json.load(f)

    out_dict = {}
    for i in out.keys():
        ase_structs = decode(encode(out[i]))
        out_dict[i] = ase_structs

    return out_dict


def write_json(out_dict, filename):
    """
    Write a Json file
    """

    str_list = [f'"{k}": {encode(v)},\n' for k, v in out_dict.items()]
    str_list[-1] = str_list[-1][:-2]
    with open(filename, "w") as wfile:
        wfile.write("{\n")
        wfile.writelines(str_list)
        wfile.write("\n}")

    return


comm = None
rank = None
size = None
origin = None
structs_dict = None
matcher = None
names = None
new_structs_dict = None


def init_param(comm_in):
    """
    Initial all the parameters
    """
    global comm
    global rank
    global size
    global origin
    global structs_dict
    global matcher
    global names
    global new_structs_dict
    global tol

    comm = comm_in
    size = comm.Get_size()
    rank = comm.Get_rank()
    origin = os.getcwd()
    structs_dict = read_json("./generation.json")
    new_structs_dict = structs_dict.copy()
    names = list(structs_dict.keys())
    tol=0.5



def combine_dict(out_list):
    """
    Return a dictionary that combines all ranks.
    """
    dic = {}
    for i in out_list:
        if i is not None:
            dic.update(i)
        else:
            continue
    return dic


def check_go_dir(dir_name):
    """
    Create a new directory and change to this directory
    """
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
        path = os.path.abspath(dir_name)
        os.chdir(path)


def scatter_list(input_list, num_cores):
    """
    Slicing the data based on the number of cores.
    """
    ave, res = divmod(len(input_list), num_cores)
    counts = [ave + 1 if p < res else ave for p in range(num_cores)]
    # determine the starting and ending indices of each sub-task
    starts = [sum(counts[:p]) for p in range(num_cores)]
    ends = [sum(counts[:p + 1]) for p in range(num_cores)]
    # save the starting and ending indices in data
    output_list = [input_list[starts[p]:ends[p]] for p in range(num_cores)]
    return output_list


def get_match_dict(new_structs_dict, name_list, ibs_target):
    """
    Runs pymatgen duplicate checks on a distributed struct dictionary.
    target : target structure to be matched
    Returns: Dictionary of matching structure names and their dftb+ energy
    """
    if len(name_list) != 0:
        match_list = []
        for name in name_list:
            if name in new_structs_dict:
                mol = new_structs_dict[name]
                ibs_mol = Structure.from_ase(mol)
                rmsd=calc_rmsd(ibs_target,ibs_mol)
                if rmsd <=tol:
                    is_match = True
                else:
                    is_match = False
                if is_match:
                    match_list.append(name)
            else:
                continue

        return match_list
    else:
        match_list = []
        return match_list


def remove_match(match_list, stru_dict):
    """
    Remove duplicates and keep only the lowest energy structure
    """
    
    for name in match_list:
        if name in stru_dict:
            del stru_dict[name]
        else:
            continue
    return stru_dict


def main():
    comm = MPI.COMM_WORLD
    init_param(comm)
    global new_structs_dict
    for i in range(len(structs_dict)):
        if i == len(structs_dict) - 1:
            continue
        if names[i] not in new_structs_dict:
            continue
        target = structs_dict[names[i]]
        ibs_target = Structure.from_ase(target)
        # get match list on each process
        if rank == 0:
            # determine the size of each sub-task
            name_list = list(structs_dict.keys())[i+1:]
            name_list = scatter_list(name_list, size)
        else:
            name_list = None
        name_list = comm.scatter(name_list, root=0)
        comm.Barrier()
        match_list = get_match_dict(new_structs_dict, name_list, ibs_target)
        match_list = comm.gather(match_list, root=0)
        if rank == 0:
            match_list = [item for sublist in match_list for item in sublist]
            if len(match_list) != 0:
                new_structs_dict = remove_match(match_list, new_structs_dict)
            else:
                new_structs_dict = new_structs_dict
        else:
            new_structs_dict = None
        comm.Barrier()
        new_structs_dict = comm.bcast(new_structs_dict, root=0)
        del match_list
        del name_list
        del target
        del ibs_target
        gc.collect()
    comm.Barrier()
    if rank == 0:
        print(f"total {len(new_structs_dict)} unique structures.")
        os.chdir(origin)
        write_json(new_structs_dict, "structures_removed.json")


if __name__ == '__main__':

    start = time.time()
    main()
    end = time.time()
    if rank==0:
        print(f"Completed in {end-start:.2f} seconds.")
    


