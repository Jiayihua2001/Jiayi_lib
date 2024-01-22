import os
from ibslib.io import read,write 
from ase.io.jsonio import encode,decode
import ase.io
import json

def write_json(out_dict, filename):
    """
    Write a Json file
    """

    str_list = [f'"{k}": {encode(v)},\n' for k, v in out_dict.items()]
    str_list[-1] = str_list[-1][:-2]
    with open(filename, "w") as wfile:
        wfile.write("{\n")
        wfile.writelines(str_list)
        wfile.write("\n}")

    return


pp=os.getcwd()
dimer_path = os.path.join(pp,'dimer_generation')
struc_dic = read(dimer_path,file_format='json')
ase_dic={}


for key, struc in struc_dic.items():
    ase_struc = struc.get_aims()
    print(ase_struc)
    ase_dic[key] = ase_struc



write_json(ase_dic,'generation_1.json')

# with open('generation_1.json','r') as f:
#     out=json.load(f)
# ase_dic={}
# for key ,value in out.items():
#     ase_struc= decode(encode(value))
#     ase_dic[key]= ase_struc
# print(list(ase_dic.items())[0])
# print(type(list(ase_dic.values())[0]))

