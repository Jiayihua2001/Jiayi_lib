import ibslib 
import json 
from compare import compare
from ibslib.io import read,write
from ase.io.jsonio import encode,decode
from ibslib import Structure

def compare_mol(mol1,mol2):
    
    result=compare(mol1,mol2,rmsd_tol=0.8)
    ret_dup=result[0]
    min_rmsd=result[1]
    if ret_dup == True:
        print(f'The two dimer structures are silimar, min_rmsd= {min_rmsd} < 0.8')
        return mol2
    elif ret_dup == False:
        return  None
    

def read_json(json_path):
    with open(json_path,'r')as f:
        out_dic=json.load(f)
    struc_dic={k:Structure.from_ase(decode(encode(v))) for k,v in out_dic.items()}
    return struc_dic

json_path='./structure_relaxed.json'
struc_dic=read_json(json_path)

mol1=read('./exp_geometry.in')   #lib stru

similar_dic={}
for e,s in struc_dic.items():
    cop=compare_mol(mol1,s)
    if cop is not None:
        similar_dic[e]=cop

for rex in similar_dic.values():
    write('mol_relax',rex,file_format='geo')
        




    
        

