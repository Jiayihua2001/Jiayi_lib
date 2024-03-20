import os
from ibslib.molecules.find_molecules import  FindMolecules
from ibslib.io import read,write

def find_mol(residues=0):   
    """
    only apply for asymmetric unit packed with 2 molecules.
    residues: int
        the number of unique molecules in the asymmetric unit, sometimes it need to be identified.
    """
    struct = read('mol.in',file_format='geo')
    n_mol=struct.get_n_atoms()
    fm=FindMolecules(residues)
    fm.calc_struct(struct)
    num_1=fm.molecules[0].get_n_atoms()
    num_2=fm.molecules[1].get_n_atoms()
    is_found= False
    if n_mol == num_1 + num_2:
        write('mol1',fm.molecules[0],file_format='geo',overwrite=True)
        write('mol2',fm.molecules[1],file_format='geo',overwrite=True)
        is_found =True

    else:
        print(f'n_mol{n_mol} != num_1 {num_1}+ num_2{num_2},you need to check the unique molecule and constrain residues')

    return is_found


def main():
    """
    please modify according to your task.
    """
    cp=os.getcwd()
    for i in os.listdir(cp):
        target_p=os.path.join(cp,i)
        try:
            os.chdir(target_p)
            os.chdir('4000')
            # is_found=find_mol()
            # if is_found :
            #     os.system('sbatch submit.sh')
            # else:
            #     pass
            os.system('sbatch submit.sh')
            os.chdir(target_p)
            os.chdir('10000')
            # is_found=find_mol()
            # if is_found :
            #     os.system('sbatch submit.sh')
            # else:
            #     pass
            os.system('sbatch submit.sh')
            os.chdir(cp)
        except:
            print('This target fail')
            os.chdir(cp)
        
if __name__ == "__main__":
    main()
    