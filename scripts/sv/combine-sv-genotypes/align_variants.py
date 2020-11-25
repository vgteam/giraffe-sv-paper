import vcf
from Bio import Align
from Bio.Seq import Seq
import argparse
from pyfaidx import Fasta


parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input gzipped VCF', required=True)
parser.add_argument('-f', help='fasta of the reference', required=True)
parser.add_argument('-o', default='out.vcf', help='Output VCF')
parser.add_argument('-l', default='none', help='log level')
parser.add_argument('-m', default=30, help='min SV size')
args = parser.parse_args()

ref_fa = Fasta(args.f, build_index=False)

aligner = Align.PairwiseAligner()
aligner.gap_score = -5
aligner.extend_gap_score = 0
aligner.match_score = 1
aligner.mismatch_score = -1
# aligner.end_gap_score = -10
# aligner.end_extend_gap_score = -1000
aligner.mode = 'global'
line_size = 200

vcf_file = args.i
vcf_out_file = args.o
min_sv_size = args.m
# log: none, all, tsv
log = args.l

# vcf_file = 'hg38-hsvlr_srdedup17_aug.NWD101273-cram.split.vcf.gz'
# vcf_out_file = 'hg38-hsvlr_srdedup17_aug.NWD101273-cram.split.al.vcf'
# min_sv_size = 30
# # log: none, all, tsv
# log = 'none'
safety_checks = False


# split long sequences
# e.g. using long unique kmers
def splitRefAlt(seq_ref, seq_alt, c=1):
    seq_ref = str(seq_ref)
    seq_alt = str(seq_alt)
    if len(seq_ref) < 500 or len(seq_alt) < 500:
        return([[seq_ref, seq_alt]])
    k = int(c * min(len(seq_ref), len(seq_alt))/200)
    step = 2 * k
    # sample kmers from ref
    kmers = []
    ref_pos = {}
    alt_pos = {}
    # sample kmers
    for ii in range(0, len(seq_ref)-k, step):
        k_cur = seq_ref[ii:(ii+k)]
        kmers.append(k_cur)
        ref_pos[k_cur] = ii
        alt_pos[k_cur] = -1
    # scan reference to flag non-unique kmers
    for ii in range(len(seq_ref)-k):
        k_cur = seq_ref[ii:(ii+k)]
        if k_cur in ref_pos and ii != ref_pos[k_cur]:
            # if found before someweher else, flag
            ref_pos[k_cur] = -1
    # scan alt for kmers
    for ii in range(len(seq_alt)-k):
        k_cur = seq_alt[ii:(ii+k)]
        if k_cur in kmers and ref_pos[k_cur] != -1:
            if alt_pos[k_cur] == -1:
                alt_pos[k_cur] = ii
            else:
                # if already found in alt, flag
                ref_pos[k_cur] = -1
    # unique kmers
    kmers_u = []
    for kk in kmers:
        if ref_pos[kk] != -1 and alt_pos[kk] != -1:
            kmers_u.append(kk)
    # check for consistent orders
    consistent_markers = True
    for ii in range(len(kmers_u)-1):
        if alt_pos[kmers_u[ii]] != -1 and \
           alt_pos[kmers_u[ii+1]] != -1 and \
           alt_pos[kmers_u[ii+1]] <= alt_pos[kmers_u[ii]]:
            # problem! next kmers is mapped before in alt
            if log == 'all':
                print("\tInconsistent markers, abort!")
            consistent_markers = False
    # split at markers
    if consistent_markers:
        ret = []
        ref_cur = 0
        alt_cur = 0
        for k_cur in kmers_u:
            if ref_cur == ref_pos[k_cur] and alt_cur == alt_pos[k_cur]:
                continue
            ret.append([seq_ref[ref_cur:ref_pos[k_cur]],
                        seq_alt[alt_cur:alt_pos[k_cur]]])
            ref_cur = ref_pos[k_cur]
            alt_cur = alt_pos[k_cur]
        if ref_cur < len(seq_ref) or alt_cur < len(seq_alt):
            ret.append([seq_ref[ref_cur:],
                        seq_alt[alt_cur:]])
        return(ret)
    else:
        return [[seq_ref, seq_alt]]


# gapless extension then alignment
def extendAlign(seq_ref, seq_alt):
    alf = ['', '', '']
    seqal_min = 0
    # try to extend from the left
    left_trim = 0
    while(left_trim < len(seq_ref) and left_trim < len(seq_alt)
          and seq_ref[left_trim] == seq_alt[left_trim]):
        left_trim += 1
    if left_trim > 0:
        if log == 'all':
            print('\t{}bp left extension'.format(left_trim))
        alf[0] = seq_ref[:left_trim]
        alf[1] = '|' * left_trim
        alf[2] = seq_alt[:left_trim]
    seq_ref = seq_ref[left_trim:]
    seq_alt = seq_alt[left_trim:]
    # try to extend from the right
    right_trim = 0
    while(right_trim < len(seq_ref) and right_trim < len(seq_alt)
          and seq_ref[-right_trim-1] == seq_alt[-right_trim-1]):
        right_trim += 1
    alf_right = ['', '', '']
    if right_trim > 0:
        if log == 'all':
            print('\t{}bp right extension'.format(right_trim))
        alf_right[0] = seq_ref[(-right_trim):]
        alf_right[1] = '|' * right_trim
        alf_right[2] = seq_alt[(-right_trim):]
        seq_ref = seq_ref[:(-right_trim)]
        seq_alt = seq_alt[:(-right_trim)]
    # sequence alignment (if necessary)
    if len(seq_ref) == 0:
        if log == 'all':
            print('\tFully extended, skip alignment')
        alf[0] += '-' * len(seq_alt)
        alf[1] += '-' * len(seq_alt)
        alf[2] += seq_alt
    elif len(seq_alt) == 0:
        if log == 'all':
            print('\tFully extended, skip alignment')
        alf[0] += seq_ref
        alf[1] += '-' * len(seq_ref)
        alf[2] += '-' * len(seq_ref)
    else:
        if log == 'all':
            print('\tSequence alignment {}x{}'.format(len(seq_ref),
                                                      len(seq_alt)))
        if log == 'tsv':
            seqal_min = min(len(seq_ref), len(seq_alt))
        als = aligner.align(seq_ref, seq_alt)
        alf_s = format(als[0]).split('\n')
        alf[0] += alf_s[0]
        alf[1] += alf_s[1]
        alf[2] += alf_s[2]
    if right_trim > 0:
        alf[0] += alf_right[0]
        alf[1] += alf_right[1]
        alf[2] += alf_right[2]
    return {'alf': alf, 'seqal_min': seqal_min}


# align REF and ALT
# eventually split large sequence
# to speed up alignment
def alignRefAlt(seq_ref, seq_alt):
    ra_split_1 = splitRefAlt(seq_ref, seq_alt)
    if log == 'all':
        print('\t' + str(len(ra_split_1)) + ' splits')
    ra_split_2 = []
    for ra in ra_split_1:
        ra_split_2 += splitRefAlt(ra[0], ra[1], c=5)
    if log == 'all':
        print('\t' + str(len(ra_split_2)) + ' splits')
    ra_split = []
    for ra in ra_split_2:
        ra_split += splitRefAlt(ra[0], ra[1], c=10)
    if log == 'all':
        print('\t' + str(len(ra_split)) + ' splits')
    alf = ['', '', '']
    nb_seqal = 0
    max_seqal_min = 0
    for ra in ra_split:
        al_o = extendAlign(ra[0], ra[1])
        if al_o['seqal_min'] > 0:
            nb_seqal += 1
            max_seqal_min = max(max_seqal_min, al_o['seqal_min'])
        alf_s = al_o['alf']
        alf[0] += alf_s[0]
        alf[1] += alf_s[1]
        alf[2] += alf_s[2]
    return {'alf': alf, 'log': {'n.splits': len(ra_split),
                                'nb.seqal': nb_seqal,
                                'max.seqal.min': max_seqal_min}}


# parse alignment into variants
def parseAlignment(alf):
    ref_pos = 0
    ii = 0
    n_vars = 0
    vars = []
    while ii < len(alf[0]):
        if alf[1][ii] == '|':
            # match: just move on base to the right
            ref_pos += 1
            ii += 1
        elif alf[0][ii] == '-' or alf[0][ii] == ' ':
            # gap on ref: insertions
            ins_seq = alf[2][ii]
            ii += 1
            while ii < len(alf[0]) and (alf[0][ii] == '-' or
                                        alf[0][ii] == ' '):
                ins_seq += alf[2][ii]
                ii += 1
            n_vars += 1
            vars.append(['INS', ref_pos, ins_seq])
            if log == 'all':
                print('\t{}bp INS at {}: {}'.format(len(ins_seq),
                                                    ref_pos,
                                                    ins_seq))
        elif alf[2][ii] == '-' or alf[2][ii] == ' ':
            # gap on alt: deletion
            del_seq = alf[0][ii]
            ii += 1
            while ii < len(alf[0]) and (alf[2][ii] == '-' or
                                        alf[2][ii] == ' '):
                del_seq += alf[0][ii]
                ii += 1
            n_vars += 1
            vars.append(['DEL', ref_pos, len(del_seq)])
            if log == 'all':
                print('\t{}bp DEL at {}: {}'.format(len(del_seq),
                                                    ref_pos,
                                                    del_seq))
            ref_pos += len(del_seq)
        else:
            # no math and no gap: SNP
            # print('\tSNP at {}: {}/{}'.format(ref_pos,
            #                                   alf[0][ii],
            #                                   alf[2][ii]))
            vars.append(['SNV', ref_pos, alf[2][ii]])
            n_vars += 1
            ii += 1
            ref_pos += 1
    return({'vars': vars, 'n_vars': n_vars})


tsv_out = '\t'.join(['{}'] * 6)
if log == 'tsv':
    print(tsv_out.format('ref.size', 'alt.size', 'n.splits',
                         'nb.seqal', 'max.seqal.min',
                         'nb.variants'))

vcf_reader = vcf.Reader(filename=vcf_file, compressed=True)
vcf_outf = open(vcf_out_file, 'w')
vcf_writer = vcf.Writer(vcf_outf, vcf_reader)

for record in vcf_reader:
    # REF and ALT alleles
    seq_ref = Seq(str(record.REF))
    seq_alt = Seq(str(record.ALT[0]))
    if len(seq_ref) < min_sv_size or len(seq_alt) < min_sv_size:
        vcf_writer.write_record(record)
        continue
    if log in ['var', 'all']:
        print('{}\tref: {}, alt: {}'.format(record.ID,
                                            len(seq_ref),
                                            len(seq_alt)))
    ara = alignRefAlt(seq_ref, seq_alt)
    alf = ara['alf']
    if log == 'all':
        for ii in range(1+int(len(alf[0])/line_size)):
            jj = min((ii+1)*line_size, len(alf[0]))
            ii = ii*line_size
            print(alf[0][ii:jj])
            print(alf[1][ii:jj])
            print(alf[2][ii:jj])
            print('+')
    # safety checks?
    if safety_checks:
        ref_al = alf[0].replace('-', '')
        alt_al = alf[2].replace('-', '')
        if ref_al != seq_ref:
            print(seq_ref)
            print(ref_al)
        if alt_al != seq_alt:
            print(seq_alt)
            print(alt_al)
    # Parse
    parse_o = parseAlignment(alf)
    # write records
    ref_pos = record.POS
    for var in parse_o['vars']:
        if var[1] == 0:
            base_pad = ref_fa[record.CHROM][ref_pos + var[1] - 2]
            base_pad = base_pad.seq.upper()
        else:
            base_pad = seq_ref[var[1] - 1]
        if var[0] == 'DEL':
            record.POS = ref_pos + var[1] - 1
            record.REF = base_pad + seq_ref[var[1]:(var[1]+var[2])]
            record.ALT[0] = base_pad
            # record.ALT[0] = seq_ref[var[1]]
            vcf_writer.write_record(record)
        if var[0] == 'INS':
            # print('{} {}'.format(len(seq_ref), var[1]))
            record.POS = ref_pos + var[1] - 1
            record.ALT[0] = base_pad + var[2]
            record.REF = base_pad
            vcf_writer.write_record(record)
    # TSV log
    if log == 'tsv':
        print(tsv_out.format(len(seq_ref),
                             len(seq_alt),
                             ara['log']['n.splits'],
                             ara['log']['nb.seqal'],
                             ara['log']['max.seqal.min'],
                             parse_o['n_vars']))

vcf_outf.close()
