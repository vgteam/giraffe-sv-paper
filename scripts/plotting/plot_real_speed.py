import json
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sklearn.linear_model import LinearRegression
import argparse

#Plot speed, memory, and runtime for 1kg, hgsvc, and yeast
#Input file (speed_report_all.tsv) must be formatted:
#graph   algorithm   reads   pairing speed   total_cpu_sec   total_wall_clock_time   max_resident_set_size(kbytes)

matplotlib.rcParams["font.family"]="sans-serif"
matplotlib.rcParams['font.sans-serif']="Arial"
matplotlib.rcParams['font.size']=12.0
#matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['axes.titlesize']="x-large"
matplotlib.rcParams['axes.labelsize']="large"
matplotlib.rcParams['axes.labelcolor'] = "black"

colors_single = { "bowtie2":"#cccc22", "bwa_mem":"#ffee99", "giraffe":"#55ccaa", "giraffe_fast":"#99ddff", "giraffe_primary":"#88bbcc", "graphaligner":"#ccdd44", "hisat2":"#ff9977", "vg_map":"#ffbbcc", "minimap2":"#eaeaea"}
colors_original={ "bowtie2":"#aaaa00", "bwa_mem":"#eedd88", "giraffe":"#44bb99", "giraffe_fast":"#99ddff", "giraffe_primary":"#77aadd", "graphaligner":"#bbcc33", "hisat2":"#ee8866", "vg_map":"#ffaabb", "minimap2":"#dddddd"}
colors_paired = { "bowtie2":"#aaaa00", "bwa_mem":"#ddcc77", "giraffe":"#33aa88", "giraffe_fast":"#77bbdd", "giraffe_primary":"#6699cc", "graphaligner":"#aabb22", "hisat2":"#dd7755", "vg_map":"#ee99aa", "minimap2":"#cccccc"}

fix_names = {"giraffe"         : "vg giraffe", 
             "giraffe_fast"    : " giraffe fast", 
             "giraffe_primary" : "giraffe primary",
             "vg_map"          : "vg map",
             "minimap2"        : "minimap2", 
             "bwa_mem"         : "bwa mem",
             "bowtie2"         : "bowtie2",
             "graphaligner"    : "graphaligner",
             "hisat2"          : "hisat2"}

def plot_histogram(stat_name, graph_name, title_name, data, units):
    #data maps {mapper_name : (single val, paired val)}
    
    fig_width = 14
    fig_height = 10
    plt.figure(figsize = (fig_width, fig_height))
    
    panel_width = .8
    panel_height = .8
    
    panel = plt.axes([0.15, 0.1, panel_width, panel_height])
    max_val = 0
    ylabels = []

    mapper_names = list(data.keys())

    #[ (value, mapper, is_paired)]
    sorted_vals = sorted(filter(lambda x : x[0] != 0 , [( data[mapper_names[i]][0],mapper_names[i], False) for i in range(len(mapper_names))] + 
                         [( data[mapper_names[i]][1],mapper_names[i], True) for i in range(len(mapper_names))])) 


    for i in range(len(sorted_vals)) :
        (val,mapper,  is_paired) = sorted_vals[i]
        y_start = i 
        max_val = max(max_val, val)
        color = colors_paired[mapper] if is_paired else colors_single[mapper]
        paired_str = " paired" if is_paired else " single"
    
        panel.add_patch(mplpatches.Rectangle([0, y_start+0.05], val, 0.9, linewidth=0, color=color, label=fix_names.get(mapper, mapper)+paired_str))
        #panel.text(val_single[0], y_start+0.5, "%.2f" % val_single[0], verticalalignment="center") 
        ylabels.append(fix_names.get(mapper, mapper)+paired_str)
    
    panel.legend(loc="lower right")
    
    
    panel.tick_params(axis='both', which='both',\
                        bottom='on',labelbottom='on',\
                        left='off', labelleft='on',\
                        right='off', labelright='off',\
                        top='off', labeltop='off')   
    panel.set_ylim(-0.05,  len(sorted_vals)+0.05)
    panel.set_xlim(0, max_val*1.02)
    panel.set_yticks([x+0.5 for x in range(len(sorted_vals))])
    panel.set_yticklabels(ylabels)
    
    #panel.set_yticks([(x*3)+1 for x in range(len(mapper_names))])
    #panel.set_yticklabels(mapper_names)
    panel.set_xlabel(stat_name + " (" + units + ")")
    panel.set_title(title_name + " " + stat_name)
    
    #plt.show()
    plt.savefig(graph_name + "_" + stat_name+"_plot.svg", dpi=400)

def main():
    with open("speed_report_all.tsv", "r") as infile:
        labels = infile.readline().split()
    
        #{graph : {mapper name : (single, paired)}} 
        all_speeds = {"1kg":{}, "hgsvc":{}, "yeast":{}}
        all_times = {"1kg":{}, "hgsvc":{}, "yeast":{}}
        all_memory = {"1kg":{}, "hgsvc":{}, "yeast":{}}
    
    
        for line in infile:
            l = line.split()
            graph=l[0] if l[0] == "1kg" or l[0] == "hgsvc" else "yeast"
            mapper=l[1]
            pairing = l[3]
            speed = float(l[4])
            total_sec = float(l[5])
            time_str = l[6].split(":")
            time = 0 if len(time_str) == 1 else float(time_str[0]) + (float(time_str[1])/60.0) + (float(time_str[2])/60.0) 
            memory = float(l[7]) / 1000000.0
        
            old_speeds = all_speeds.get(graph, {}).get(mapper, (0,0))
            old_times =   all_times.get(graph, {}).get(mapper, (0,0))
            old_memory = all_memory.get(graph, {}).get(mapper, (0,0))
            new_speeds = (speed, old_speeds[1])  if pairing == "single" else (old_speeds[0], speed)
            new_times =  (time, old_times[1])    if pairing == "single" else (old_times[0], time)
            new_memory = (memory, old_memory[1]) if pairing == "single" else (old_memory[0], memory)
            all_speeds[graph][mapper] = new_speeds
            all_times [graph][mapper]  = new_times
            all_memory[graph][mapper] = new_memory
    

    plot_histogram("Speed",  "1kg", "1KGP/GRCh37", all_speeds["1kg"], "reads/second/thread")
    plot_histogram("Runtime","1kg", "1KGP/GRCh37", all_times["1kg"],  "hours")
    plot_histogram("Memory", "1kg", "1KGP/GRCh37", all_memory["1kg"], "gbytes")
    plot_histogram("Speed",  "hgsvc", "HGSVC/GRCh38", all_speeds["hgsvc"], "reads/second/thread")
    plot_histogram("Runtime","hgsvc", "HGSVC/GRCh38", all_times ["hgsvc"],  "hours")
    plot_histogram("Memory", "hgsvc", "HGSVC/GRCh38", all_memory["hgsvc"], "gbytes")
    plot_histogram("Speed",  "yeast", "yeast", all_speeds["yeast"], "reads/second/thread")
main()
