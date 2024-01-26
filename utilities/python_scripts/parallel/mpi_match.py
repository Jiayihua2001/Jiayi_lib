import json

from mpi4py import MPI
from ase.io import read
from ase.io.jsonio import encode, decode
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher


def read_json(jsonfile):
    """
    Load Json file to return a dictionary storing ase structures
    """

    with open(jsonfile, "r") as f:
        stru_dict = json.load(f)

    ase_dict = {k: decode(encode(v)) for k, v in stru_dict.items()}
    return ase_dict



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


def get_match_list(structs_dict, id_list, pmg_target, matcher):
    """
    Runs pymatgen duplicate checks on a distributed struct dictionary.
    target : target structure to be matched
    Returns: List of matching structure names/ids
    """
    if len(id_list) != 0:
        match_list = []
        for key in id_list:
            xtal = structs_dict[key]
            pmg_xtal = AseAtomsAdaptor.get_structure(xtal)
            is_match = matcher.fit(pmg_target, pmg_xtal)
            if is_match:
                match_list.append(key)
    else:
        match_list = []

    return match_list


def main():
    # Initialise MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    exp = read("exp_geometry.in")
    structs_dict = read_json("./structures.json")
    matcher = StructureMatcher()
    pmg_target = AseAtomsAdaptor.get_structure(exp)

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
    match_list = get_match_list(structs_dict, struct_id_list, pmg_target, matcher)
    match_list = comm.gather(match_list, root=0)
    if rank == 0:
        # flatten
        match_list = [item for sublist in match_list for item in sublist]
        print("# of unique structrues", len(match_list))
        print("unique structrues", match_list)


if __name__ == "__main__":
    main()
