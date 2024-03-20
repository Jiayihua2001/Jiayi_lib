#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import json
import time
from mpi4py import MPI
from ase.io import read
from ase.io.jsonio import encode, decode
from compare import compare
from ibslib import Structure

def read_json(jsonfile):
    """
    Load Json file to return a dictionary storing ase structures
    """

    with open(jsonfile, "r") as f:
        stru_dict = json.load(f)

    struc_dict = {k: decode(encode(v)) for k, v in stru_dict.items()}
    return struc_dict

def compare_mol(mol1,mol2,tol=0.6):
    
    result=compare(mol1,mol2,rmsd_tol=tol)
    ret_dup=result[0]
    min_rmsd=result[1]
    if ret_dup == True:
        return min_rmsd, mol2
    elif ret_dup == False:
        return None, None

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


def get_match_list(structs_dict,exp_dimer, id_list):
    """
    Runs pymatgen duplicate checks on a distributed struct dictionary.
    target : target structure to be matched
    Returns: List of matching structure names/ids
    """
    if len(id_list) != 0:
        match_dic ={}
        for key in id_list:
            dimer = Structure.from_ase(structs_dict[key]) 
            min_rmsd,cop=compare_mol(exp_dimer,dimer)
            if cop is not None:
                match_dic[min_rmsd]=cop
                print(f'Matched structure: energy: {key} ,min_rmsd: {min_rmsd} < 0.6')
    
    else:
        match_dic ={}

    return match_dic

def write_json(out_dict, filename):
    """
    Write a Json file 
    out_dic should be ase_dic
    """

    str_list = [f'"{k}": {encode(v.get_ase_atoms())},\n' for k, v in out_dict.items()]
    str_list[-1] = str_list[-1][:-2]
    with open(filename, "w") as wfile:
        wfile.write("{\n")
        wfile.writelines(str_list)
        wfile.write("\n}")

    return



def main():
    # Initialise MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    
    dimer_dict = None
    if rank == 0:
        exp_dimer = Structure.from_ase(read('exp_geometry.in'))   #relaxed
        json_dir = './generation.json'
        dimer_dict = read_json(json_dir)
        struct_id_list = list(dimer_dict.keys())
        struct_id_list = scatter_list(struct_id_list, size)
    else:
        struct_id_list = None

    dimer_dict = comm.bcast(dimer_dict,root=0)
    exp_dimer = comm.bcast(exp_dimer,root=0)
    # dispatches struct_id_list from main process across all processes
    struct_id_list = comm.scatter(struct_id_list, root=0)
    # or comm.barrier(). Not necessary.

    match_dic= get_match_list(dimer_dict,exp_dimer, struct_id_list)
    match_dic = comm.gather(match_dic, root=0)
    if rank == 0:
        # flatten
        match_dic = {k:v for d in match_dic for k,v in d.items()}
        sorted_match_dic = {k:match_dic[k] for k in sorted(match_dic)}
        if len(match_dic) !=0:
            write_json(sorted_match_dic,'sorted_match.json')
    


if __name__ == "__main__":
    main()
