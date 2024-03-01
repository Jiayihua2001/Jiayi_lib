from pymatgen.io.cif import CifParser
from ibslib import Structure
from ibslib.io import read,write
import os
import write_control



def read_cif(file_path):

    # Initialize the parser with the CIF file
    parser = CifParser(file_path)

    # Get the structure
    structure = parser.get_structures()[0]
    return structure

def relax_xtal(submit_p):
    write('geometry.in',struct,file_format='geo' )
    write_control.main()
    os.system(f'cp {submit_p} .')
    os.system('mv submit_relax.sh submit.sh')
    os.system('sbatch submit.sh')
    return 

def is_converge(aims_path):
    """
    check if converged
    """
    converge = False
    for line in reversed(list(open(aims_path))):
        if 'Have a nice day.' in line:  
            converge = True
            break
    return converge

cp=os.getcwd()
submit_r_p=os.path.abspath('submit_relax.sh')
submit_gnrs_p = os.path.abspath('submit_gnrs.sh')
conf_p =os.path.abspath('ui.conf')
for s in os.listdir(cp):
    tmp=os.path.join(cp,s)
    os.chdir(tmp)
    os.system(f'cp {submit_gnrs_p} .')
    os.system('mv submit_gnrs.sh submit.sh')
    os.system(f'cp {conf_p} .')
    for c in os.listdir(tmp):
        try:
            struct_p = read_cif(c)
            struct = Structure.from_pymatgen(struct_p)
            if not os.isdir('relax_xtal'):
                os.mkdir('relax_xtal') 
            os.chdir('relax_xtal')
            if 'aims.out' in os.listdir(os.path.join(tmp,'relax_xtal')):
                aims_path=os.path.abspath('aims.out')
                if is_converge(aims_path):
                    os.system(f'cp geometry.in.next_step {tmp}')   
                os.chdir(tmp)
                os.system('mv geometry.in.next_step exp_geometry.in')
            else:
                relax_xtal(submit_r_p)
                os.chdir(tmp)
                print('The relaxation task of exp_xtal has been submitted')
        except:
            print(f"{c} : attempt not succeeded")
        


