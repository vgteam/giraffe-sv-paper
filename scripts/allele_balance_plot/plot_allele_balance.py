import io
import json
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
from scipy import stats
import re
import functools

INDEL_SIZE_LIMIT = 40
colors_original={ "bowtie2":"#aaaa00", "bwa_mem":"#eedd88", "giraffe":"#44bb99", "giraffe_fast":"#99ddff", "giraffe_primary":"#77aadd", "graphaligner":"#bbcc33", "hisat2":"#ee8866", "vg_map":"#ffaabb", "minimap2":"#dddddd"}


def mean(array):
    if len(array) == 0:
        return 0
    return float(sum(array)) / float(len(array))

def error(array):
    return stats.sem(array)

def get_fractions(filename):
    #filename is a vcf file with allele depth
    #for each mapper, map indel len to [ref count, alt count, [fraction of alt for each site]]
    all_counts = {"giraffe":{}, "map":{}, "bwa":{}}
    indel_count = 0
    snp_count = 0
    with open(filename, "r") as vcf_file:
        for line in vcf_file:
            if line[0] != "#":
                tokens = line.split()
                ref_allele = tokens[3]
                alt_alleles = tokens[4].split(",")
                all_alleles = [ref_allele] + alt_alleles
                qual = float(tokens[5])
                info = tokens[7].split(";")
                info_dict = dict()
                for x in info:
                    if "=" in x:
                        info_dict[x.split("=")[0]] = x.split("=")[1]

                gt_labels = tokens[8].split(":")
                gt_vals_giraffe = tokens[9].split(":")
                gt_vals_map = tokens[10].split(":")
                gt_vals_bwa = tokens[11].split(":")
                gts = {"giraffe" : {gt_labels[i] : gt_vals_giraffe[i] for i in range(len(gt_labels))},
                       "map"     : {gt_labels[i] : gt_vals_map[i]     for i in range(len(gt_labels))},
                       "bwa"      : {gt_labels[i] : gt_vals_bwa[i]    for i in range(len(gt_labels))}}
                giraffe_gt = re.split("\||/",gts["giraffe"]["GT"])
                map_gt = re.split("\||/",gts["map"]["GT"])
                bwa_gt = re.split("\||/",gts["bwa"]["GT"])


                
                giraffe_is_het =  len(giraffe_gt) == 2 and giraffe_gt[0] != giraffe_gt[1]
                map_is_het     =  len(map_gt) == 2 and map_gt[0] != map_gt[1]
                bwa_is_het     =  len(bwa_gt) == 2 and bwa_gt[0] != bwa_gt[1]

                is_het = giraffe_is_het and map_is_het and bwa_is_het
                pass_filter = (int(info_dict.get("MQ", 40)) >=40 and int(info_dict.get("DP", 25)) >= 25)



                if is_het and pass_filter:
                    for mapper in ["giraffe", "map", "bwa"]:

                        read_counts = gts[mapper]["AD"].split(",")

                        alt_allele = alt_alleles[0]
                        indel_len =  len(alt_allele) - len(ref_allele)
                        if indel_len >= INDEL_SIZE_LIMIT:
                            indel_size = INDEL_SIZE_LIMIT
                        elif indel_len <= -INDEL_SIZE_LIMIT:
                            indel_size = -INDEL_SIZE_LIMIT
                        else:
                            indel_size=indel_len

                        ref_count = 0
                        alt_count = 0

                        for genotype_num in range(len(read_counts)):
                            count = 0 if read_counts[genotype_num] == "." else int(read_counts[genotype_num])

                            if genotype_num == 0:
                                #ref
                                ref_count += count
                            else:
                                #alt
                                alt_count += count
                                if count != 0:
                                    if indel_size not in all_counts[mapper]:
                                        all_counts[mapper][indel_size] = [0,0,[]]

                                    all_counts[mapper][indel_size][1] += count

                        
                        if (ref_count != 0 or alt_count != 0):
                            alt_fraction = float(alt_count) / float(ref_count+alt_count)
                            if indel_size not in all_counts[mapper]:
                                all_counts[mapper][indel_size] = [0,0,[]]
                            all_counts[mapper][indel_size][0] += ref_count
                            all_counts[mapper][indel_size][2].append( alt_fraction )


        #Each of these is a list of (indel length, average over total counts, average of fractions, error)
        all_fractions_giraffe = sorted([ (key, float(val[1]) / float(val[0]+val[1]), mean(val[2]), error(val[2])) for key, val in all_counts["giraffe"].items()])
        all_fractions_map     = sorted([ (key, float(val[1]) / float(val[0]+val[1]), mean(val[2]), error(val[2])) for key, val in all_counts["map"].items()])
        all_fractions_bwa     = sorted([ (key, float(val[1]) / float(val[0]+val[1]), mean(val[2]), error(val[2])) for key, val in all_counts["bwa"].items()])
        print(all_fractions_giraffe)
        print(all_fractions_map)
        print(all_fractions_bwa)
        with open ("indel_fractions.json", "w") as out_file:
            json.dump([all_fractions_giraffe, all_fractions_map, all_fractions_bwa], out_file)
        
        return (all_fractions_giraffe, all_fractions_map, all_fractions_bwa)


(all_fractions_giraffe, all_fractions_map, all_fractions_bwa) = get_fractions("all_calls.vcf")
with open("indel_fractions.json", "r") as in_file:
    all_fractions = json.load(in_file)


#[indel len, avg of read counts, avg of fractions, error, number of sites]
all_fractions_giraffe = all_fractions[0]
all_fractions_map = all_fractions[1]
all_fractions_bwa = all_fractions[2]

indel_lengths_giraffe   = [x[0] for x in all_fractions_giraffe]
indel_fractions_giraffe = [x[1] for x in all_fractions_giraffe]
indel_error_giraffe     = [x[3] for x in all_fractions_giraffe]
indel_lengths_map       = [x[0] for x in all_fractions_map]
indel_fractions_map     = [x[1] for x in all_fractions_map]
indel_error_map         = [x[3] for x in all_fractions_map]
indel_lengths_bwa       = [x[0] for x in all_fractions_bwa]
indel_fractions_bwa     = [x[1] for x in all_fractions_bwa]
indel_error_bwa         = [x[3] for x in all_fractions_bwa]


#Plot figure


#Make this look more like the R roc plots

matplotlib.rcParams["font.family"]="sans-serif"
matplotlib.rcParams['font.sans-serif']="Arial"
matplotlib.rcParams['font.size']=12.0
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['axes.titlesize']="x-large"
matplotlib.rcParams['axes.labelsize']="large"
matplotlib.rcParams['axes.labelcolor'] = "black"

fig_width = 15
fig_height = 5
plt.figure(figsize=(fig_width,fig_height))
panel_width=0.9
panel_height=0.8
panel=plt.axes([0.05, 0.1, panel_width, panel_height])



panel.plot([-90, 50], [0.5, 0.5], color='black',linestyle='dashed') #Line at 0.5

#Combine indels > 40bp

panel.scatter(indel_lengths_bwa, indel_fractions_bwa, label="bwa mem", color=colors_original["bwa_mem"])
panel.plot(indel_lengths_bwa, indel_fractions_bwa, color=colors_original["bwa_mem"])
panel.errorbar(indel_lengths_bwa, indel_fractions_bwa, yerr=indel_error_bwa, color=colors_original["bwa_mem"])

panel.scatter(indel_lengths_map, indel_fractions_map, label="vg map", color=colors_original["vg_map"])
panel.plot(indel_lengths_map, indel_fractions_map, color=colors_original["vg_map"])
panel.errorbar(indel_lengths_map, indel_fractions_map, yerr=indel_error_map, color=colors_original["vg_map"])

panel.scatter(indel_lengths_giraffe, indel_fractions_giraffe, label="vg giraffe", color=colors_original["giraffe"])
panel.plot(indel_lengths_giraffe, indel_fractions_giraffe, color=colors_original["giraffe"])
panel.errorbar(indel_lengths_giraffe, indel_fractions_giraffe, yerr=indel_error_giraffe, color=colors_original["giraffe"])



panel.legend(loc="lower left",frameon=False)


panel.set_xlabel("Insertion or deletion length")
panel.set_ylabel("Fraction of alternate allele")
#panel.set_ylim(0,1)
panel.set_xlim(-INDEL_SIZE_LIMIT-2, INDEL_SIZE_LIMIT+2)
panel.set_xticks([-40,-30,-20,-10,0,10,20,30,40])
panel.set_xticklabels(["<-40", "-30", "-20", "-10", "0","10","20","30",">40"])
    
plt.savefig("allele_balance.svg")
