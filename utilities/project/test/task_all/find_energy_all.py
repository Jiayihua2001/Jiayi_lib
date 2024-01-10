#after successful running
#write new_json and get energy

#!/usr/bin/python
# Prompt the user for the file name
import re
import pandas as pd
import os
import json
# Function to extract words behind 'hello' in a file

n_all=input('total number of structures: ')
project_name=input('project_name: ')
data=[]
status=[]
data_s=[]
json_filename = input("Please enter the file name of your JSON file: ")
with open(json_filename, "r") as f:
    out = json.load(f)

dir_energy={}
def extract_info(n_all):
    for i in range(0,int(n_all)+1):
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
                pass
            else:
                status.extend(data_s)
        else:
            print(f"the {out.keys()[i]} task didn't completed")
            status.append('unconverge')
            pass
        dic_energy[out.keys()[i]]=status[i]
    return dic_energy

Dic_energy=extract_info(n_all)
# Create a Pandas DataFrame from the extracted information
os.chdir('..')

with open("Total_energy.json",'w',encoding='utf8') as f:
    json.dump(Dic_energy,f,indent=6) 

df = pd.DataFrame({"index:":Dic_energy.keys(),'energy:': Dic_energy.values()})
print(df)   
