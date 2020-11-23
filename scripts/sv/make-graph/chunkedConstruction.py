import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna
from pyfaidx import Fasta
import vcf
import subprocess
import os
import shutil


parser = argparse.ArgumentParser(description='Construcion without augmentation (e.g. for chunk that take too long augmenting).')
parser.add_argument('-r', help='the reference fasta', required=True)
parser.add_argument('-c', help='the VCF with variants to vg construct', required=True)
parser.add_argument('-a', help='the VCF with variants to should have been augmented', required=True)
parser.add_argument('-o', help='the output augmented .vg graph', required=True)
parser.add_argument('-d', help='debug mode: keep temporary files', action='store_true')
args = parser.parse_args()

vg = 'vg'

# temporary files
temp_vg = args.o + '_temp.vg'
temp_fa = args.o + '_temp.fa'
temp_vcf = args.o + '_temp.vcf'
# log file for all the vg commands
log_err = open(args.o + '.log', 'w')

# chunk region from the output .vg filename (e.g. mpmap_augment_chunks/augmented/chunk_chr16_3839328_4394815.vg)
fn = args.o.rstrip('.vg').split('/')[-1]
fn = fn.split('_')
chunk_seqn = fn[1]
chunk_s = int(fn[2])
chunk_e = int(fn[3])
print('Chunk {} from {} to {}.'.format(chunk_seqn, chunk_s, chunk_e))

## extract reference sequence
ref = Fasta(args.r)
ref_seq = ref[chunk_seqn][chunk_s:chunk_e]
ref_farec = SeqRecord(MutableSeq(ref_seq.seq.upper(), generic_dna), id=chunk_seqn,
                 description='')
SeqIO.write(ref_farec, temp_fa, "fasta")

# write the VCF for this chunk in a temp file
vcfo = open(temp_vcf, 'w')
nb_var_to_construct = 0

# extract variants to add using vg construct
vcfi = open(args.c, 'r')
vcf_reader = vcf.Reader(vcfi)
vcf_writer = vcf.Writer(vcfo, vcf_reader)
for record in vcf_reader:
    if record.CHROM == chunk_seqn and record.POS + 32 < chunk_e and record.POS - 32 > chunk_s:
        record.POS = record.POS - chunk_s
        vcf_writer.write_record(record)
        nb_var_to_construct += 1
vcfi.close()

# extract variants that should have been augmented
vcfi = open(args.a, 'r')
vcf_reader = vcf.Reader(vcfi)
for record in vcf_reader:
    if record.CHROM == chunk_seqn and record.POS + 32 < chunk_e and record.POS - 32 > chunk_s:
        record.POS = record.POS - chunk_s
        vcf_writer.write_record(record)
        nb_var_to_construct += 1
vcfi.close()

vcfo.close()

## vg construct
print('Make graph with vg construct adding ' + str(nb_var_to_construct) + ' variants...')
cons_cmd = [vg, 'construct', '-t', '1', '-r', temp_fa, '-a', '-f', '-S']
if(nb_var_to_construct>0):
    ## add variants from the VCF
    cons_cmd += ['-v', temp_vcf]
with open(temp_vg, 'w') as outf:
    subprocess.run(cons_cmd, stdout=outf, stderr=log_err, check=True)


shutil.copy(temp_vg, args.o)

# close log
log_err.close()
print('Check the log at ' + args.o + '.log')

# remove temporary files if not in "debug" mode
if not args.d:
    if(os.path.isfile(temp_vg)):
        os.remove(temp_vg)
    if(os.path.isfile(temp_fa)):
        os.remove(temp_fa)
    if(os.path.isfile(temp_fa + '.fai')):
        os.remove(temp_fa + '.fai')
    if(os.path.isfile(temp_vcf)):
        os.remove(temp_vcf)
