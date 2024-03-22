import json
from ase.io.jsonio import decode,encode
import copy
import numpy as np
import os
from ase.io import read,write
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher



def matcher(exp,struct):
    """
    exp struc must be ase struc
    """
    settings = {"stol": 0.5, "ltol": 0.5, "angle_tol": 10, 'scale': False}
    matcher = StructureMatcher(**settings)
    pmg_target = AseAtomsAdaptor.get_structure(exp)
    pmg_xtal = AseAtomsAdaptor.get_structure(struct)
    if matcher.fit(pmg_target, pmg_xtal):
        rmsd = matcher.get_rms_dist(pmg_target, pmg_xtal)
    return rmsd


def get_dict_from_json(struc_path,hash_index,exp):
    with open(struc_path,'r') as f:
        out=json.load(f)

    for keys,values in out.items():
        if keys == hash_index:
            ase_struc = decode(encode(values))
            Rmsd=matcher(exp,ase_struc)
            os.chdir('bestmatch')
            write(f'bestmatch{keys}_RMSD:{Rmsd}.in',ase_struc)
            os.chdir(cp)
    return 

cp = os.getcwd()
path_json=os.path.join(cp,'structures/symm_rigid_press/structures.json')
match_list=['c82a00d70b7fa18', 'f6d916352c3aa35', '5dd3408edb201a3']
exp = read("exp_geometry.in")
os.mkdir('bestmatch')
for i in match_list:
    get_dict_from_json(path_json,i,exp)