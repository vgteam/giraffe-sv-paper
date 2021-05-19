import fileinput

# read sample names
samples = {}
with open('sample.list', 'r') as inf:
    for line in inf:
        line = line.rstrip()
        samples[line] = True

# read stream
samp_cols = []
for line in fileinput.input():
    line = line.rstrip()
    # skip headers and write appropriate tsv headers
    if line[0] == '#':
        if line[:6] == '#CHROM':
            line = line.split('\t')
            headers = []
            for ii in range(len(line)):
                if line[ii] in samples:
                    samp_cols.append(ii)
                    headers.append(line[ii])
            print('id\t' + '\t'.join(headers))
        continue
    # subset samples
    out_gts = []
    line = line.split('\t')
    nb_non_ref = 0
    for ii in samp_cols:
        gt = line[ii].split('|')
        gt_code = -1
        if gt.count('0') == 2:
            gt_code = 0
        elif '1' in gt and '0' in gt:
            gt_code = 1
        elif gt.count('1') == 2:
            gt_code = 2
        if gt_code > 0:
            nb_non_ref += 1
        out_gts.append(str(gt_code))
    if nb_non_ref > 3:
        print(line[2] + '\t' + '\t'.join(out_gts))

