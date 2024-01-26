import os
from ase.io import read, write
import numpy as np
from ase.calculators.aims import Aims
import json
from mpi4py import MPI
from ase.io.jsonio import encode, decode
from ase.optimize import BFGS
import time

### this file is just a pseudocode for scattering task on diff ranks use mpi


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
       json.dump(str_list,wfile)
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


def energy_calc(xtal):
    with open('aims_settings.json', 'r') as f:
        aims_set = json.load(f)
    calc = Aims(parallel=False, **aims_set)
    xtal.calc = calc
    energy = xtal.get_potential_energy()
    return energy

def optimizer(xtal):
    with open('aims_settings.json', 'r') as f:
        aims_set = json.load(f)
    calc = Aims(parallel=False, **aims_set)
    xtal.calc = calc
    opt = BFGS(xtal)
    opt.run(fmax=0.2)
    energy = xtal.get_potential_energy()
    return xtal,energy

def sort_energy(energy_dic):
    min_energy = min(energy_dic.keys())
    min_energy_struc = energy_dic[min_energy]
    write('min_energy_geometry.in',min_energy_struc)
    print(f'min_energy:{min_energy}')
    return 

def main():
    # Initialise MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    start_time = time.perf_counter()
    structs_dict = read_json("./generation_test.json")

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

    dimer_list=[]
    for key in struct_id_list:
        xtal = structs_dict[key]
        # run you task 
        run_task()
        dimer_list.append()


    comm.Barrier()
    dimer_list = comm.gather(dimer_list, root=0)
    if rank == 0:
        # flatten
        dimer_list = [k for sublist in dimer_list for k in sublist]
        sort_energy(dimer_list) 
    end_time  = time.perf_counter()
    elapsed_time = end_time-start_time
    print(f'Time used in processing 20 structures is: {elapsed_time}')
if __name__ == "__main__":
    main()


