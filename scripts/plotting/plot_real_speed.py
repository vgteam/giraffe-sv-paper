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
#graph   algorithm   reads   pairing speed   total_wall_clock_time   max_resident_set_size(kbytes)



matplotlib.rcParams["font.family"]="sans-serif"
matplotlib.rcParams['font.sans-serif']="Arial"
matplotlib.rcParams['font.size']=12.0
#matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['axes.titlesize']="x-large"
matplotlib.rcParams['axes.labelsize']="large"
matplotlib.rcParams['axes.labelcolor'] = "black"

#(light, original, dark)
colors={ "bowtie2_single"         : "#cccc22", "bowtie2"         : "#aaaa00", "bowtie2_paired"         : "#aaaa00", 
         "bwa_mem_single"         : "#ffee99", "bwa_mem"         : "#eedd88", "bwa_mem_paired"         : "#ddcc77",  
         "minimap2_single"        : "#eaeaea", "minimap2"        : "#dddddd", "minimap2_paired"         : "#cccccc",
         "giraffe_single"         : "#66ddbb", "giraffe"         : "#44bb99", "giraffe_paired"         : "#44bb99", 
         "giraffe_cover_single"   : "#66ddbb", "giraffe_cover"   : "#44bb99", "giraffe_cover_paired"   : "#44bb99", 
         "giraffe_sampled_single" : "#66ddbb", "giraffe_sampled" : "#44bb99", "giraffe_sampled_paired" : "#44bb99", 
         "giraffe_full_single"    : "#007755", "giraffe_full"    : "#118866", "giraffe_full_paired"    : "#229977",
         "giraffe_fast_single"    : "#99ddff", "giraffe_fast"    : "#99ddff", "giraffe_fast_paired"    : "#77bbdd", 
         "giraffe_primary_single" : "#88bbcc", "giraffe_primary" : "#77aadd", "giraffe_primary_paired" : "#6699cc",
         "graphaligner_single"    : "#ccdd44", "graphaligner"    : "#bbcc33", "graphaligner_paired"    : "#aabb22", 
         "vg_map_single"          : "#ffbbcc", "vg_map"          : "#ffaabb", "vg_map_paired"          : "#ee99aa",
         "vg_map_primary_single"  : "#bb6677", "vg_map_primary"  : "#aa5566", "vg_map_primary_paired"  : "#994455",
         "hisat2_single"          : "#ff9988", "hisat2"          : "#ee8877", "hisat2_paired"          : "#dd7766",
         "hisat2_sens_single"     : "#cc6655", "hisat2_sens"     : "#bb5544", "hisat2_sens_paired"     : "#aa4433",
         "hisat2_vsens_single"    : "#993322", "hisat2_vsens"    : "#882211", "hisat2_vsens_paired"    : "#771100"}

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

def plot_histogram(stat_name, graph_name, data, read_name):
    #data maps {mapper_name : (single val, paired val)}

    
    fig_width = 14
    fig_height = 16
    plt.figure(figsize = (fig_width, fig_height))
    
    panel_width = .8
    panel_height = .8
    
    panel = plt.axes([0.15, 0.1, panel_width, panel_height])
    max_val = 0
    ylabels = []
    panel.set_axisbelow(True)
    plt.grid()

    mapper_names = list(data.keys())

    #[ (value, mapper)]
    sorted_vals = sorted(filter(lambda x : x[0] != 0 , 
        [( data[mapper_names[i]], mapper_names[i]) for i in range(len(mapper_names))] )) 


    for i in range(len(sorted_vals)) :
        (val,mapper) = sorted_vals[i]
        y_start = i 
        #if extra_val != 0:
        #    panel.add_patch(mplpatches.Rectangle([val, y_start+0.05], extra_val, 0.9, linewidth=0, color='#bb5533'))

        if (stat_name == "memory") and (val > 240):
            #Hacky way of cutting off the memory when we ran out

            panel.text(0.5, y_start+0.5, "Out of Memory", verticalalignment="center", fontsize=15) 
        else :
            max_val = max(max_val, val)
            color = colors[mapper]
    
            panel.add_patch(mplpatches.Rectangle([0, y_start+0.05], val, 0.9, linewidth=0, color=color, label=mapper))
        ylabels.append(mapper)
    
    #panel.legend(loc="lower right")
    
    
    panel.tick_params(axis='both', which='both',\
                        bottom='on',labelbottom='on',\
                        left='off', labelleft='on',\
                        right='off', labelright='off',\
                        top='off', labeltop='off')   
    panel.set_ylim(-0.05,  len(sorted_vals)+0.05)
    if stat_name == "speed":
        panel.set_xlim(0, 8000)
    elif stat_name == "memory":
        panel.set_xlim(0, 110)
    elif stat_name == "runtime":
        panel.set_xlim(0, 50)
    else:
        panel.set_xlim(0, max_val*1.02)

    panel.set_yticks([x+0.5 for x in range(len(sorted_vals))])
    panel.set_yticklabels(ylabels)
    
    unit_names = {"speed" : "reads/second/thread",
             "runtime" : "hours",
             "memory" : "gbytes" }
    graph_names = { "1000gp" : "1KGP/GRCh38",
                    "hgsvc" : "HGSVC/GRCh38",
                    "yeast" : "yeast" }

    #panel.set_yticks([(x*3)+1 for x in range(len(mapper_names))])
    #panel.set_yticklabels(mapper_names)
    panel.set_xlabel(stat_name + " (" + unit_names[stat_name] + ")")
    panel.set_title(graph_names[graph_name] + " " + stat_name)
    
    #plt.show()
    plt.savefig(graph_name + "_" + stat_name+"_"+read_name+"_plot.svg", dpi=400)

def main():
    parser = argparse.ArgumentParser(description="Plot stats from a tsv file (named speed_report_[reads].tsv) that is formatted:\
        graph   algorithm   reads   pairing speed   total_wall_clock_time   max_resident_set_size(kbytes)")
    parser.add_argument("graph", help="which graph [1000gp/hgsvc/yeast]")
    parser.add_argument("reads", help="which reads [novaseq6000/hiseqxten/hiseq2500]")
    parser.add_argument("stat", help="which stat to plot [speed/memory/runtime]")
    args = parser.parse_args()

    with open("speed_report_"+args.reads+".tsv", "r") as infile:
        labels = infile.readline().split()
    
        #{graph : {mapper name : (single, paired)}} 
        all_stats = {}
    
    
        for line in infile:
            l = line.split()
            if l[0] == args.graph :
                mapper=l[1]
                pairing = l[3]
                if args.stat == "speed":
                    stat= float(l[4])
                elif args.stat == "runtime":
                    time_str = l[5].split(":")
                    if len(time_str) == 2:
                        time_str = ["0"] + time_str
                    stat = 0 if len(time_str) == 1 else float(time_str[0]) + (float(time_str[1])/60.0) + (float(time_str[2])/3600.0) 
                elif args.stat == "memory":
                    stat = float(l[6]) / 1000000.0
                else:
                    print("Bad stat")
                    assert(False)
        
                all_stats[mapper+"_"+pairing] = stat
    

    plot_histogram(args.stat,  args.graph, all_stats, args.reads)

main()
