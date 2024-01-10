import matplotlib.pyplot as plt 
from ibslib.analysis.pltlib import format_ticks
from ibslib.plot import labels_and_ticks
from ibslib.plot.ml import pca
from matplotlib.ticker import FormatStrFormatter
import copy
import numpy as np

def plot_spg_hist(
    spg_values,
    general_spg_values=[],
    special_spg_values=[],
    ax=None,
    figname="",
    general_bar_kw = 
        {
            "width": 0.8,
            "color": "tab:blue",
            "edgecolor": "k"
        },
    special_bar_kw = 
        {
            "width": 0.8,
            "color": "tab:orange",
            "edgecolor": "k",
        },
    exp_spg_arrow =
        {
        
        },
    xlabel_kw = 
        {
            "xlabel": "Space Group",
            "fontsize": 16,
            "labelpad": 0,
        },
    ylabel_kw = 
        {
            "ylabel": "Number of Structures",
            "fontsize": 16,
            "labelpad": 10, 
        },
    xticks = 
        {
            "xlim": [],
            "xticks_kw":
                {
                    "ticks": [],
                },
            "xticklabels_kw": 
                {
                    "labels": [],
                    "fontsize": 6,
                    "rotation": 90,
                },
        },
    yticks = 
        {
            "ylim": [],
            "yticks_kw":
                {
                    "ticks": [],
                },
            "yticklabels_kw": 
                {
                    "labels": [],
                    "fontsize": 12,
                },
            "FormatStrFormatter": "%.0f"
        
        },    
    ):
    """
    Plots a Genarris space group histogram
    
    Arguments
    ---------
    spg_values: list
        List or array of volume values from the pool of structures. It's easy 
        to obtain this by reading in the directory using ibslib.io.read and
        then getting values using ibslib.analysis.get.
    general_spg_values: list
        List or array of space group numbers for which the molecule would sit 
        on a general position. If neither general_spg_values or 
        special_spg_values are provided, all spg_values are assumed to 
        be general positions.
    special_spg_values: list
        List or array of space group numbers for which the molecule would sit 
        on a specical position.  
    ax: matplotlib.pyplot.Axes 
        Provide if this plot is supposed to be drawn on an already existing
        Axes. If ax is None, a new figure is created. This gives the user
        the most flexibility. If you would like the figure to be a certain
        size, please initialize the figure first yourself and feed in the 
        Axes object. For example:
            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111)
    exp_spg_arrow: dict
        Dictionary of arguments to be use to draw an arrow pointing towards
        the experimental space group.
        
    """
    arguments = locals()
    arguments_copy = {}
    
    for key,value in arguments.items():
        if key == "spg_values":
            arguments_copy[key] = value
        elif key == "ax":
            arguments_copy[key] = value
        else:
            arguments_copy[key] = copy.deepcopy(value)
    arguments = arguments_copy
    arguments["spg_values"] = list(arguments["spg_values"])
    
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    if len(general_spg_values) > 0 or len(special_spg_values) > 0:
        original_bins = general_spg_values + special_spg_values
        min_spg = min(original_bins)
        max_spg = max(original_bins)
        original_bins_final = [min_spg-1]
        original_bins_final += original_bins
        original_bins_final.append(max_spg+1)
        original_bins_final.sort()
    else:
        min_spg = min(spg_values)
        max_spg = max(spg_values)
        original_bins_final = np.arange(min_spg-1, max_spg+2, 1).tolist()
        
    ## Run hist the first time to get histogram setting values
    temp_fig = plt.figure()
    temp_ax = temp_fig.add_subplot(111)
    hist = temp_ax.hist(spg_values, bins=original_bins_final, edgecolor="k",
                   align="mid")
    plt.close()
    
    #### Create maping from spg values to consecutive integers so that the 
    ## plot only shows all allowed spg values consecutively.
    all_spg = hist[1][0:-1]
    int_map = np.arange(0,len(all_spg))
    
    ### Now we're going to make a bar graph of the values we obtained, but 
    ## using different formatting for general and special positions
    if len(special_spg_values) > 0:
        special_height = []
        for value in special_spg_values:
            idx = np.where(hist[1] == value)[0]
            special_height.append(hist[0][idx][0])
        
        ## Refer to int_map to get x values
        special_int_idx = np.searchsorted(all_spg, special_spg_values)
        special_int_map = int_map[special_int_idx]
        
        ax.bar(special_int_map, special_height, **special_bar_kw)

    
    if len(general_spg_values) > 0:
        general_height = []
        for value in general_spg_values:
            idx = np.where(hist[1] == value)[0]
            general_height.append(hist[0][idx][0])
            
        ## Refer to int_map to get x values
        general_int_idx = np.searchsorted(all_spg, general_spg_values)
        general_int_map = int_map[general_int_idx]
        
        ax.bar(general_int_map, general_height, **general_bar_kw)
    
    # Handle behavior when both values are zero
    if len(special_spg_values) == 0 and len(general_spg_values) == 0:
        ax.bar(all_spg, hist[0], **general_bar_kw)
    
    # Default for xtick values is over the observed range of tick values
    if len(arguments["xticks"]["xticks_kw"]["ticks"]) == 0 and \
        len(general_spg_values) == 0 and \
        len(special_spg_values) == 0:
        min_spg = min(spg_values)
        max_spg = max(spg_values)
        arguments["xticks"]["xticks_kw"]["ticks"] = np.arange(min_spg,max_spg+1,1)
    elif len(arguments["xticks"]["xticks_kw"]["ticks"]) == 0:
        ### Set xticks to int map
        arguments["xticks"]["xticklabels_kw"]["labels"] = [str(x) for x in all_spg] 
        arguments["xticks"]["xticks_kw"]["ticks"] = int_map.tolist()
        
    format_ticks(ax)
    labels_and_ticks(ax, arguments["xlabel_kw"], arguments["ylabel_kw"],
                     arguments["xticks"],arguments["yticks"])
    
    if len(general_spg_values)>0 and len(special_spg_values) > 0:
        for i,xtick in enumerate(ax.get_xticklabels()):
            if i in general_int_idx:
                xtick.set_color("tab:blue")
            elif i in special_int_idx:
                xtick.set_color("tab:orange")
    else:
        pass


    if len(figname) > 0:
        fig.savefig(figname) 
    ### Cannot return ax object as an argument
    del(arguments["ax"])
    return arguments

import matplotlib.pyplot as plt
import numpy as np
import json
import os
# Other necessary imports

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

file_path_sg='sg.json'
#cp= os.getcwd()
cp='/ocean/projects/mat210008p/jhuanga/gnrsAIMNet/examples/target11_0.84_1.0'
path_json=os.path.join(cp,file_path_sg)

sg_dict=read_dict_from_file(path_json)
# Prepare your data
spg_values = list(sg_dict.values()) # Example data

# Call the function
# Since 'plot_spg_hist' returns a dictionary, we're not capturing its return value here
plot_spg_hist(spg_values,figname='myplot')

# To display the plot
plt.show()

# Or to save the plot to a file, call the function with a 'figname'
#plot_spg_hist(spg_values,figname="my_plot.png")
