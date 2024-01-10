from ibslib.io  import read

from ibslib.structure import Structure
from ase.io.jsonio import encode,decode
import json
import os
from ibslib.analysis import get
from ibslib.plot.genarris import plot_volume_hist,plot_spg_hist
import os


current_path = os.getcwd()

def get_dict_from_json(struc_path):
    with open(struc_path,'r') as f:
        out=json.load(f)
    struct_dict={}
    volume_dict={}
    sg_dict={}
    for keys,values in out.items():
        ase_struc = decode(encode(values))
        struc=Structure.from_ase(ase_struc)
        volume_value=struc.get_unit_cell_volume()
        sg=struc.get_space_group()
        struct_dict[keys]=struc
        volume_dict[keys]=volume_value
        sg_dict[keys]=sg
    return volume_dict,sg_dict


def save_dict_to_json_file(data_dict, file_path):
    try:
        with open(file_path, 'w') as file:
            # Serialize dict and write to the file
            json.dump(data_dict, file, indent=4)
    except IOError:
        print(f"Error: Unable to write to the file {file_path}.")


def read_dict_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            # Read the file content
            file_content = file.read()
            # Convert the JSON string to a dictionary
            data_dict = json.loads(file_content)
            return data_dict
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
    except json.JSONDecodeError:
        print("Error: The file content is not in valid JSON format.")

### need to be fixed later
# def recreate_plot(arguments):
#     fig = plt.figure()
#     ax = fig.add_subplot(111)

#     # Call the original function with the saved arguments
#     plot_volume_hist(ax=ax, **arguments)



#Usage
#step1:get dict from json
path_json=os.path.join(current_path,'structures/generation/structures.json')
volume_dict,sg_dict=get_dict_from_json(path_json)

# # save the dict to json (if needed)
# file_path_v = 'volume.json' # Replace with your desired file path
# file_path_sg='sg.json'
# save_dict_to_json_file(volume_dict, file_path_v)
# save_dict_to_json_file(sg_dict, file_path_sg)

# #read from json (if needed)
# volume_dict=read_dict_from_file(file_path_v)
# sg_dict=read_dict_from_file(file_path_sg)

#step2:dict to list,plot the graph
volume_list=list(volume_dict.values())   ### !!! attention:plot input must be a list instead of a dictionary
sg_list=list(sg_dict.values())
plot_volume_hist(volume_list,exp_volume=753.299,figname='volume_plot')  #only if you give the filename ,it can save the plot  , exp_volume in A^3

plot_spg_hist(sg_list,exp_spg_arrow=14)



