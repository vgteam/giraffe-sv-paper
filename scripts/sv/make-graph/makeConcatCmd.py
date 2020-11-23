import argparse


parser = argparse.ArgumentParser(description='Prepare the concat '
                                 'chunks command')
parser.add_argument('-b', help='the bed file from vg chunk', required=True)
parser.add_argument('-c', help='the chromosome name', required=True)
parser.add_argument('-o', help='the output vg file', required=True)
args = parser.parse_args()

# temporary concatenated graphs
concat_vg = args.o + '_temp.vg'
concat_vg2 = args.o + '_temp2.vg'

# list of chunks to concatenate
aug_chunks_paths = []

print('set -e')

# write intermediate concat every X file
# to avoid command being too long
cpt = 0
max_cpt = 500

# read input bed file
bed_in = open(args.b, 'r')
for line in bed_in:
    line = line.rstrip().split('\t')
    if line[0] == args.c:
        cpt += 1
        aug_chunks_paths.append(line[3])
        # remove alt paths
        print('vg paths -r -Q chr -v {} > {}'.format(line[3], concat_vg))
        print('mv {} {}'.format(concat_vg, line[3]))
        if cpt == max_cpt:
            # write an intermediate command and continue
            print('vg concat -p {} > {}'.format(' '.join(aug_chunks_paths),
                                                concat_vg))
            # prepare for potentially more commands
            print('mv {} {}'.format(concat_vg, concat_vg2))
            aug_chunks_paths = [concat_vg2]
            cpt = 0
bed_in.close()

# last concat command if there are enough chunks saved
if len(aug_chunks_paths) > 1:
    print('vg concat -p {} > {}'.format(' '.join(aug_chunks_paths),
                                        concat_vg))

# rename final concatenated file
print('mv {} {}'.format(concat_vg, args.o))

# remove temporary file
print ('rm -f ' + concat_vg2)
