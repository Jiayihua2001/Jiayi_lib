from pymatgen.core.structure import Structure, Molecule
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifParser
from ibslib.io import read,write
from ibslib.molecules import com,orientation
import numpy as np
from ibslib.molecules import FindMolecules
import ibslib
import ase

# still need further improvement
### the problem: the molecule not complete:
### possible method :   must take distance into account.

def extract_asymmetric_unit(structure):
    """
    Extracts the asymmetric unit from a given crystal structure.

    Args:
        structure (Structure): A pymatgen Structure object representing the crystal unit cell.

    Returns:
        Structure: A pymatgen Structure object containing the asymmetric unit.
    """
    
    # Analyze the structure to get the symmetrized structure, which contains the asymmetric unit
    structure.make_supercell(scaling_matrix=[3, 3, 3])
    finder = SpacegroupAnalyzer(structure)
    symmetrized_structure = finder.get_symmetrized_structure()
    
    # Extract asymmetric unit sites from the symmetrized structure
    asymmetric_sites = symmetrized_structure.equivalent_sites

    # Assemble the asymmetric unit from unique sites
    unique_sites = []
    for site_list in asymmetric_sites:
        unique_sites.append(site_list[0])  # Take the first site from each equivalent site list
    
    # Create a new structure from the unique sites
    asymmetric_unit_structure = Structure.from_sites(unique_sites)

    return asymmetric_unit_structure

def read_cif(file_path):
    """
    Reads a CIF file and returns the structure.

    Args:
        file_path (str): Path to the CIF file.

    Returns:
        Structure: A pymatgen Structure object of the CIF file.
    """
    parser = CifParser(file_path)
    structure = parser.get_structures()[0]
    return structure

def center_dimer(struct, choose_min=True, fm=None):
    """
    Centers the system for one the molecules in the dimer such that the COM is 
    at the origin and the axes defined by the moment of inertia tensor are  
    oriented with the origin. This is done by translating and rotation the 
    entire dimer system.
    
    Arguments
    ---------
    struct: Structure
        Structure to adjust
    choose_min: bool
        If True, will center on molecule closest to the origin.
        If False, will center on molecule furtherst from origin. 
    fm: FindMolecules
        FindMolecules object that has been calculated using Structure. This will
        save time. Otherwise, default will calculate the molecules in the 
        Structure. 
    
    """
    if fm == None:
        fm = FindMolecules()
        fm.calc(struct)
        
    # check_dimer(struct, fm=fm)
    
    com_list = np.vstack([com(x) for x in fm.molecules])
    com_list = np.linalg.norm(com_list, axis=-1)
    
    if choose_min:
        chosen_idx = np.argmin(com_list)
        other_idx = np.argmax(com_list)
    else:
        chosen_idx = np.argmax(com_list)
        other_idx = np.argmin(com_list)
        
    trans = com(fm.molecules[chosen_idx])
    rot = orientation(fm.molecules[chosen_idx])
    
    geo = struct.get_geo_array()
    geo = geo - trans
    geo = np.dot(geo, rot.T)
    struct.from_geo_array(geo, struct.geometry["element"])
    return struct

def to_mol(struc):
    species = [site.species_string for site in struc]
    coords = [site.coords for site in struc]
    # Create the Molecule object
    molecule = Molecule(species, coords)
    return molecule

# Example usage
file_path = './2306641.cif'
structure = read_cif(file_path)
asymmetric_unit = extract_asymmetric_unit(structure)
mol=to_mol(asymmetric_unit)
struc=ibslib.structure.Structure.from_pymatgen(mol)
center_mol= center_dimer(struc)
write('asy.in',center_mol,file_format='geo')


