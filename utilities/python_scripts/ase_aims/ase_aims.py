#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from ase.io import read, write
import numpy as np
from ase.calculators.aims import Aims
import json
# from mpi4py import MPI
from ase.io.jsonio import encode, decode
from ase.optimize import BFGS
import time

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


def energy_calc(xtal):
    with open('aims_settings.json', 'r') as f:
        aims_set = json.load(f)
    calc = Aims(parallel=True, **aims_set)
    xtal.calc = calc
    energy = xtal.get_potential_energy()
    return energy

def optimizer(xtal,fmax_v):
    """
    by defalut fmax_v =0.02
    for pre -relax to select the xtal ,can use 0.2 
    """
    with open('aims_settings.json', 'r') as f:
        aims_set = json.load(f)
    calc = Aims(parallel=False, **aims_set)
    xtal.calc = calc
    opt = BFGS(xtal)
    opt.run(fmax=fmax_v)  #pre-relax
    energy = xtal.get_potential_energy()
    return xtal,energy

def sort_energy(energy_dic):
    sorted_dic={k:energy_dic[k] for k in sorted(energy_dic)[:10]}
    write_json(sorted_dic,'sorted.json')
    min_energy = min(energy_dic.keys())
    min_energy_struc = energy_dic[min_energy]
    write('min_energy_geometry.in',min_energy_struc)
    print(f'min_energy:{min_energy}')
    return 
    
def main():
    start_time = time.perf_counter()
    ase_dict = read_json("./generation.json")
    energy_dic={}
    # for key,value in ase_dict.items():
    #     xtal,energy = optimizer(value,fmax_v=0.2)
    #     energy_dic[energy] = xtal           (too time consuming!!!)
    for key,value in ase_dict.items():
        energy = energy_calc(value)
        energy_dic[energy] = value

    sort_energy(energy_dic)
    end_time  = time.perf_counter()
    elapsed_time = end_time-start_time

    print(f'Time used is: {elapsed_time}')

if __name__ == "__main__":
    main()


