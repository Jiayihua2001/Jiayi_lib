
#after successful running
#write new_json and get energy

#!/usr/bin/python
# Prompt the user for the file name
import re
import pandas as pd
import os
import json

# input info
project_name=input('project_name: ')
key_list=[]
status=[]
data_s=[]
dic_energy={}
#read files enter keys

json_filename = input("Please enter the file name of your JSON file: ")
path=f'/ocean/projects/mat210008p/jhuanga/project/{project_name}/{json_filename}'
with open(path, "r") as f:
    out = json.load(f)
for keys in out.keys():
    key_list.append(keys)
dir_energy={}
#func of extracting info
def extract_info():
    for i in range(0,len(key_list)):#change back len(key_list)
        dir_name='structure%i'%i
        dir_lo = f'/ocean/projects/mat210008p/jhuanga/project/{project_name}/task_all/{dir_name}'
        os.chdir(dir_lo)
        file_name = 'aims.out'
        task_success='Have a nice day.'
    with open(file_name, 'r') as file:
        content = file.read()
        file.seek(0)
        lines=file.readlines()
        target_sentence='Final output of selected total energy values:'
    
    #judge and extract
        if re.search(task_success,content):
            for line_number, line in enumerate(lines, start=-1):
                if target_sentence in line:
                    start_index = max(0, line_number+6)
                    end_index = min(len(lines), line_number+7)
                    # Extract and print the surrounding context
                    context = lines[start_index:end_index]
                    context_string=''.join(context)
                    energy=float(re.split(r'\s',context_string)[-3])
                    data_s.append(energy)
                else:
                    pass
            if data_s==[]:
                print(f'Total Energy  not found in the {i} file.')
                status.append('not found')
            else:
                status.extend(data_s)

        else:
            print(f"the {key_list(i)} task didn't completed")
            status.append('unconverge')
        
  
        dic_energy[key_list[i]]=status[i]
    return dic_energy

Dic_energy = extract_info()

# Create a Pandas DataFrame from the extracted information
os.chdir('..')

with open("Total_energy.json",'w',encoding='utf8') as f:
    json.dump(Dic_energy,f,indent=6) 

df = pd.DataFrame({"index:":Dic_energy.keys(),'energy:': Dic_energy.values()})
print(df)   
