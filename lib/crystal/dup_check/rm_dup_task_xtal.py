import gc
import os
import json
import time
import numpy as np

from mpi4py import MPI
from ase.io import read
from ase.io.jsonio import encode, decode
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher

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

    comm = comm_in
    size = comm.Get_size()
    rank = comm.Get_rank()
    origin = os.getcwd()
    structs_dict = read_json("./structures.json")
    new_structs_dict = structs_dict.copy()
    settings = {"stol": 0.5, "ltol": 0.5, "angle_tol": 10}
    matcher = StructureMatcher(**settings)
    names = list(structs_dict.keys())


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


def get_match_dict(new_structs_dict, name_list, pmg_target, matcher):
    """
    Runs pymatgen duplicate checks on a distributed struct dictionary.
    target : target structure to be matched
    Returns: Dictionary of matching structure names and their dftb+ energy
    """
    if len(name_list) != 0:
        match_list = []
        for name in name_list:
            if name in new_structs_dict:
                xtal = new_structs_dict[name]
                pmg_xtal = AseAtomsAdaptor.get_structure(xtal)
                is_match = matcher.fit(pmg_target, pmg_xtal)
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
        pmg_target = AseAtomsAdaptor.get_structure(target)
        # get match list on each process
        if rank == 0:
            # determine the size of each sub-task
            name_list = list(structs_dict.keys())[i+1:]
            name_list = scatter_list(name_list, size)
        else:
            name_list = None
        name_list = comm.scatter(name_list, root=0)
        comm.Barrier()
        match_list = get_match_dict(new_structs_dict, name_list, pmg_target, matcher)
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
        del pmg_target
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
    print(f"Completed in {end-start:.2f} seconds.")
    


