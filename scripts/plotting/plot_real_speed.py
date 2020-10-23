import json
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sklearn.linear_model import LinearRegression
import argparse

#Plot speed, memory, and runtime
#Input file is must be formatted:
#graph   algorithm   reads   pairing speed   total_cpu_sec   total_wall_clock_time   max_resident_set_size(kbytes)

matplotlib.rcParams["font.family"]="sans-serif"
matplotlib.rcParams['font.sans-serif']="Arial"
matplotlib.rcParams['font.size']=12.0
#matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['axes.titlesize']="x-large"
matplotlib.rcParams['axes.labelsize']="large"
matplotlib.rcParams['axes.labelcolor'] = "black"

#(light, original, dark)
colors={ "bowtie2"         : ("#cccc22", "#aaaa00", "#aaaa00"), 
         "bwa_mem"         : ("#ffee99", "#eedd88", "#ddcc77"),  
         "minimap2"        : ("#eaeaea", "#dddddd", "#cccccc"),
         "giraffe"         : ("#66ddbb", "#44bb99", "#44bb99"), 
         "giraffe_cover"   : ("#66ddbb", "#44bb99", "#44bb99"), 
         "giraffe_sampled" : ("#66ddbb", "#44bb99", "#44bb99"), 
         "giraffe_full"    : ("#007755", "#118866", "#229977"),
         "giraffe_fast"    : ("#99ddff", "#99ddff", "#77bbdd"), 
         "giraffe_primary" : ("#88bbcc", "#77aadd", "#6699cc"),
         "graphaligner"    : ("#ccdd44", "#bbcc33", "#aabb22"), 
         "vg_map"          : ("#ffbbcc", "#ffaabb", "#ee99aa"),
         "hisat2"          : ("#ff9988", "#ee8877", "#dd7766"),
         "hisat2_sens"     : ("#cc6655", "#bb5544", "#aa4433"),
         "hisat2_vsens"    : ("#993322", "#882211", "#771100")}

fix_names = {"giraffe"         : "vg giraffe", 
             "giraffe_fast"    : "giraffe fast", 
             "giraffe_primary" : "giraffe primary",
             "giraffe_sampled" : "giraffe sampled",
             "giraffe_cover"   : "giraffe cover",
             "giraffe_full"    : "giraffe full",
             "vg_map"          : "vg map",
             "minimap2"        : "minimap2", 
             "bwa_mem"         : "bwa mem",
             "bowtie2"         : "bowtie2",
             "graphaligner"    : "graphaligner",
             "hisat2"          : "hisat2",
             "hisat2_sens"     : "hisat2 sensitive",
             "hisat2_vsens"    : "hisat2 very sensitive"}

def plot_histogram(stat_name, graph_name, title_name, data, units, read_name):
    #data maps {mapper_name : (single val, paired val)}
    
    fig_width = 14
    fig_height = 15
    plt.figure(figsize = (fig_width, fig_height))
    
    panel_width = .8
    panel_height = .8
    
    panel = plt.axes([0.15, 0.1, panel_width, panel_height])
    max_val = 0
    ylabels = []

    mapper_names = list(data.keys())

    #[ (value, mapper, is_paired)]
    sorted_vals = sorted(filter(lambda x : x[0] != 0 , 
        [( data[mapper_names[i]][0],data[mapper_names[i]][2],mapper_names[i], False) for i in range(len(mapper_names))] + 
        [( data[mapper_names[i]][1],data[mapper_names[i]][3],mapper_names[i], True) for i in range(len(mapper_names))])) 


    for i in range(len(sorted_vals)) :
        (val,extra_val,mapper,  is_paired) = sorted_vals[i]
        y_start = i 
        #if extra_val != 0:
        #    panel.add_patch(mplpatches.Rectangle([val, y_start+0.05], extra_val, 0.9, linewidth=0, color='#bb5533'))

        if (stat_name == "Memory") and (val > 240):
            #Hacky way of cutting off the memory when we ran out

            panel.text(0.5, y_start+0.5, "Out of Memory", verticalalignment="center", fontsize=15) 
        else :
            max_val = max(max_val, val)
            color = colors[mapper][2] if is_paired else colors[mapper][0]
            paired_str = " paired" if is_paired else " single"
    
            panel.add_patch(mplpatches.Rectangle([0, y_start+0.05], val, 0.9, linewidth=0, color=color, label=fix_names.get(mapper, mapper)+paired_str))
        ylabels.append(fix_names.get(mapper, mapper)+paired_str)
    
    #panel.legend(loc="lower right")
    
    
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
    plt.savefig(graph_name + "_" + stat_name+"_"+read_name+"_plot.svg", dpi=400)

def main(read_name):
    with open("speed_report_"+read_name+".tsv", "r") as infile:
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
            if len(time_str) == 2:
                time_str = ["0"] + time_str
            time = 0 if len(time_str) == 1 else float(time_str[0]) + (float(time_str[1])/60.0) + (float(time_str[2])/3600.0) 
            memory = float(l[7]) / 1000000.0
            if len(l) == 9:
                time_str = l[8].split(":")
                extra_time = 0 if len(time_str) == 1 else float(time_str[0]) + (float(time_str[1])/60.0) + (float(time_str[2])/3600.0) 
            else:
                extra_time=0

        
            #(single val, paired val, extra val single, extra val paired)
            old_speeds = all_speeds.get(graph, {}).get(mapper, (0,0,0,0))
            old_times =   all_times.get(graph, {}).get(mapper, (0,0,0,0))
            old_memory = all_memory.get(graph, {}).get(mapper, (0,0,0,0))
            new_speeds = (speed, old_speeds[1],0,0)  if pairing == "single" else (old_speeds[0], speed,0,0)
            new_times =  (time, old_times[1],extra_time, old_times[3])  if pairing == "single" else (old_times[0], time,old_times[2], extra_time)
            new_memory = (memory, old_memory[1],0,0) if pairing == "single" else (old_memory[0], memory,0,0)
            all_speeds[graph][mapper] = new_speeds
            all_times [graph][mapper]  = new_times
            all_memory[graph][mapper] = new_memory
    

    plot_histogram("Speed",  "1kg", "1KGP/GRCh37", all_speeds["1kg"], "reads/second/thread", read_name)
    plot_histogram("Runtime","1kg", "1KGP/GRCh37", all_times["1kg"],  "hours", read_name)
    plot_histogram("Memory", "1kg", "1KGP/GRCh37", all_memory["1kg"], "gbytes", read_name)
    plot_histogram("Speed",  "hgsvc", "HGSVC/GRCh38", all_speeds["hgsvc"], "reads/second/thread", read_name)
    plot_histogram("Runtime","hgsvc", "HGSVC/GRCh38", all_times ["hgsvc"],  "hours", read_name)
    plot_histogram("Memory", "hgsvc", "HGSVC/GRCh38", all_memory["hgsvc"], "gbytes", read_name)
    plot_histogram("Speed",  "yeast", "yeast", all_speeds["yeast"], "reads/second/thread", read_name)

main("novaseq600")
