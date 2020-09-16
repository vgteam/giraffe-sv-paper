import io
import json
import argparse
parser = argparse.ArgumentParser(description='combine alignments by marking a read in primary_file as correct if it is mapped correctly in secondary_file')
parser.add_argument('primary_file', help='json of primary alignments')
parser.add_argument('secondary_file', help='json of secondary alignments')
parser.add_argument('out_file', help='output json')

args = parser.parse_args()

primary_file =   open(args.accumulate(args.primary_file), "r")

secondary_reads = dict()

with open(args.secondary_file, "r") as secondary_file:
    for line in secondary_file.readlines():
        read_dict = json.loads(line)
        name = read_dict["name"]
        correct = "correctly_mapped" in read_dict
        secondary_reads[name] = secondary_reads.get(name, False) or correct


with open(args.out_file, "w") as out_file:
    for line in primary_file.readlines():
        read_dict = json.loads(line)
        name = read_dict["name"]
        if name in secondary_reads:
            if secondary_reads[name]:
                read_dict["correctly_mapped"] = True
        out_file.write(json.dumps(read_dict)+"\n")

primary_file.close()
