import json

from mpi4py import MPI
from ase.io import read
from ase.io.jsonio import encode, decode
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher



def main():
   
    exp = read("exp_geometry.in")
    struct = read('geometry.in')
    matcher = StructureMatcher(ltol=2,stol=3,angle_tol=30)
    pmg_target = AseAtomsAdaptor.get_structure(exp)
    pmg_xtal = AseAtomsAdaptor.get_structure(struct)
    if matcher.fit(pmg_target, pmg_xtal):
        rmsd = matcher.get_rms_dist(pmg_target, pmg_xtal)
        print(rmsd)

if __name__ == "__main__":
    main()
