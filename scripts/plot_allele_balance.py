import io
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
from scipy import stats
import re

FILTER_VARIANTS = False
INDEL_SIZE_LIMIT = 40

def mean(array):
    return float(sum(array)) / float(len(array))

def error(array):
    return stats.sem(array)

def get_fractions(filename):
    #filename is a vcf file with allele depth
    all_counts = {"giraffe":{}, "map":{}, "bwa":{}}
    indel_count = 0
    snp_count = 0
    with open(filename, "r") as vcf_file:
        for line in vcf_file:
            if line[0] != "#":
                tokens = line.split()
                ref_allele = tokens[3]
                alt_alleles = tokens[4].split(",")
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
                gt_vals_na19239 = tokens[12].split(":")
                gt_vals_na19239_only = tokens[13].split(":")
                gts = {"giraffe"      : {gt_labels[i] : gt_vals_giraffe[i]      for i in range(len(gt_labels))},
                       "map"          : {gt_labels[i] : gt_vals_map[i]          for i in range(len(gt_labels))},
                       "bwa"          : {gt_labels[i] : gt_vals_bwa[i]          for i in range(len(gt_labels))},
                       "na19239"      : {gt_labels[i] : gt_vals_na19239[i]      for i in range(len(gt_labels))},
                       "na19239_only" : {gt_labels[i] : gt_vals_na19239_only[i] for i in range(len(gt_labels))}}

                na19239_genotype = re.split("\||/",gts["na19239"]["GT"])
                giraffe_genotype = re.split("\||/",gts["giraffe"]["GT"])
                map_genotype = re.split("\||/",gts["map"]["GT"])
                bwa_genotype = re.split("\||/",gts["bwa"]["GT"])
                is_het = len(giraffe_genotype) == 2 and giraffe_genotype[0] != giraffe_genotype[1] and \
                         len(map_genotype) == 2 and map_genotype[0] != map_genotype[1] and \
                         len(bwa_genotype) == 2 and bwa_genotype[0] != bwa_genotype[1]

                #This variant is in the graph if NA19239 didn't have it
                #is_in_graph = gts["na19239_only"]["GT"] == "0/0" or gts["na19239_only"]["GT"] == "0|0" or \
                #              gts["na19239_only"]["GT"] == "./." or gts["na19239_only"]["GT"] == ".|."

                pass_filter = not FILTER_VARIANTS or (int(info_dict.get("MQ", 40)) >=40 and int(info_dict.get("DP", 25)) >= 25)
                if  is_het and pass_filter:

                    for mapper in ["giraffe", "map", "bwa"]:

                        genotype = re.split("\||/", gts[mapper]["GT"])
                        allele_num = int(genotype[1]) if genotype[1] != "0" else int(genotype[0])

                        counts = gts[mapper]["AD"].split(",")
                        ref_count = 0 if counts[0] == "." else int(counts[0])
                        alt_count = 0 if counts[allele_num] == "." else int(counts[allele_num])
                        alt_fraction = float(ref_count) / float(ref_count+alt_count)

                        indel_len = len(alt_alleles[allele_num-1]) - len(ref_allele)
                        if indel_len >= INDEL_SIZE_LIMIT:
                            indel_size = INDEL_SIZE_LIMIT
                        elif indel_len <= -INDEL_SIZE_LIMIT:
                            indel_size = -INDEL_SIZE_LIMIT
                        else:
                            indel_size=indel_len
                        
                        if ref_count!=0 or alt_count!=0:
                            if indel_size not in all_counts[mapper]:
                                all_counts[mapper][indel_size] = [0,0,[]]
                            all_counts[mapper][indel_size][0] += ref_count
                            all_counts[mapper][indel_size][1] += alt_count
                            all_counts[mapper][indel_size][2].append(alt_fraction)    

        all_fractions_giraffe = sorted([ (key, float(val[1]) / float(val[0]+val[1]), error(val[2])) for key, val in all_counts["giraffe"].items()])
        all_fractions_map     = sorted([ (key, float(val[1]) / float(val[0]+val[1]), error(val[2])) for key, val in all_counts["map"].items()])
        all_fractions_bwa     = sorted([ (key, float(val[1]) / float(val[0]+val[1]), error(val[2])) for key, val in all_counts["bwa"].items()])
        print(all_fractions_giraffe)
        print(all_fractions_map)
        print(all_fractions_bwa)
        with open ("indel_fractions.json", "w") as out_file:
            json.dump([all_fractions_giraffe, all_fractions_map, all_fractions_bwa], out_file)
        
        return (all_fractions_giraffe, all_fractions_map, all_fractions_bwa)


#(all_fractions_giraffe, all_fractions_map, all_fractions_bwa) = get_fractions("1kg_filtered_all.vcf")
with open("indel_fractions.json", "r") as in_file:
    all_fractions = json.load(in_file)

all_fractions_giraffe = all_fractions[0]
all_fractions_map = all_fractions[1]
all_fractions_bwa = all_fractions[2]

indel_lengths_giraffe   = [x[0] for x in all_fractions_giraffe]
indel_fractions_giraffe = [x[1] for x in all_fractions_giraffe]
indel_error_giraffe     = [x[2] for x in all_fractions_giraffe]
indel_lengths_map       = [x[0] for x in all_fractions_map]
indel_fractions_map     = [x[1] for x in all_fractions_map]
indel_error_map         = [x[2] for x in all_fractions_map]
indel_lengths_bwa       = [x[0] for x in all_fractions_bwa]
indel_fractions_bwa     = [x[1] for x in all_fractions_bwa]
indel_error_bwa         = [x[2] for x in all_fractions_bwa]


#Plot figure

fig_width = 15
fig_height = 5
plt.figure(figsize=(fig_width,fig_height))
panel_width=0.9
panel_height=0.8
panel=plt.axes([0.05, 0.1, panel_width, panel_height])


panel.plot([-90, 50], [0.5, 0.5], color='black',linestyle='dashed') #Line at 0.5

#Combine indels > 40bp

panel.scatter(indel_lengths_giraffe, indel_fractions_giraffe, label="vg giraffe", color="blue")
panel.plot(indel_lengths_giraffe, indel_fractions_giraffe, color="blue")
panel.errorbar(indel_lengths_giraffe, indel_fractions_giraffe, yerr=indel_error_giraffe, color="blue")

panel.scatter(indel_lengths_map, indel_fractions_map, label="vg map", color="red")
panel.plot(indel_lengths_map, indel_fractions_map, color="red")
panel.errorbar(indel_lengths_map, indel_fractions_map, yerr=indel_error_map, color="red")

panel.scatter(indel_lengths_bwa, indel_fractions_bwa, label="bwa mem", color="green")
panel.plot(indel_lengths_bwa, indel_fractions_bwa, color="green")
panel.errorbar(indel_lengths_bwa, indel_fractions_bwa, yerr=indel_error_bwa, color="green")


panel.legend(loc="lower left")


panel.set_xlabel("Insertion or deletion length")
panel.set_ylabel("Fraction of alternate allele")
#panel.set_ylim(0,1)
panel.set_xlim(-INDEL_SIZE_LIMIT-2, INDEL_SIZE_LIMIT+2)
    
plt.savefig("allele_balance_plot_1kg.png", dpi=500)
