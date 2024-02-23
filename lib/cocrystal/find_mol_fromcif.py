from ibslib.molecules.find_molecules import  FindMolecules
from pymatgen.io.cif import CifParser
from ibslib import Structure
from ibslib.io import read,write

def read_cif(file_path):

    # Initialize the parser with the CIF file
    parser = CifParser(file_path)

    # Get the structure
    structure = parser.get_structures()[0]
    return structure

struct_p = read_cif('/ocean/projects/mat210008p/jhuanga/projects/darlug/1531947.cif')
struct = Structure.from_pymatgen(struct_p)
fm=FindMolecules()
fm.calc_struct(struct)
write('mol1',fm.molecules[0],file_format='geo',overwrite=True)
write('mol2',fm.molecules[1],file_format='geo',overwrite=True)
