#!/usr/bin/env python
# -*- coding: utf-8 -*-
import ibslib 
import json 
from mpi4py import MPI
from compare import compare
from ase.io.jsonio import encode,decode
from ibslib import Structure

def compare_mol(mol1,mol2,tol):
    
    result=compare(mol1,mol2,rmsd_tol=tol)
    ret_dup=result[0]
    return ret_dup
    

def read_json(json_path):
    with open(json_path,'r')as f:
        out_dic=json.load(f)
    struc_dic={k:Structure.from_ase(decode(encode(v))) for k,v in out_dic.items()}
    return struc_dic


def dup_removal(structs_dict,id_list,tol):
    if len(id_list) != 0:
        struc_dic={k:structs_dict[k] for k in id_list}
        to_remove = set()  # A set to keep track of keys to remove
        keys = id_list # List of keys to iterate over
        for i in range(len(keys)):
            k1 = keys[i]
            v1 = struc_dic[k1]
            for j in range(i + 1, len(keys)):
                k2 = keys[j]
                v2 = struc_dic[k2]
                if compare_mol(v1,v2,tol):
                    to_remove.add(k2)  # Mark for removal

        # Remove marked keys
        for k in to_remove:
            id_list.remove(k)

        unique_key=id_list

    else:
        unique_key = []

    return unique_key


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


def write_json(out_dict, filename):
    """
    Write a Json file 
    out_dic should be ase_dic
    """
    ase_dict = {k:v.get_ase_atoms() for k,v in out_dict.items()}
    str_dict = {k:encode(v) for k, v in ase_dict.items()}
    with open(filename,'w') as wf:
        json.dump(str_dict,wf,indent=6)
    return


def main():
    # Initialise MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    structs_dict = read_json("./generation.json")

    if rank == 0: 
        # determine the size of each sub-task
        struct_id_list = list(structs_dict.keys())
        struct_id_list = scatter_list(struct_id_list, size)
    else:
        struct_id_list = None

    # dispatches struct_id_list from main process across all processes
    struct_id_list = comm.scatter(struct_id_list, root=0)
    # or comm.barrier(). Not necessary.
    comm.Barrier()
    unique_list = dup_removal(structs_dict, struct_id_list,tol=0.6)
    unique_list = comm.gather(unique_list, root=0)
    if rank == 0:
        # flatten
        tol=0.6   # remember to change the params of dup_removal() accordingly
        unique_list = [item for sublist in unique_list for item in sublist]
        # unique_list_final=dup_removal(structs_dict,unique_list)
        unique_dict={k:structs_dict[k] for k in unique_list}
        write_json(unique_dict,f'dup_tol{tol}.json')
        

if __name__ == "__main__":
    main()



    


