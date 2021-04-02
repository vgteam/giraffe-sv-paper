#!/usr/bin/env python3
"""
assess_mappability.py: determine the fraction of k-mers that could be correctly
mapped by a perfect mapper.

A perfect mapper maps k-mers correctly with probability equal to their
multiplicity in the reference.

We select a random sample of n k-mers from random strands of a FASTA, excluding
those containing Ns or other ambiguous bases. We then scan the FASTA and count
how many times each of those k-mers appears, across both strands.

We compute mappability as:

((unique k-mers) + (2-copy k-mers)/2 + (3-copy k-mers)/3 + ...)/(total k-mers)

Re-uses sample code and documentation from 
<http://users.soe.ucsc.edu/~karplus/bme205/f12/Scaffold.html>
"""

import argparse, sys, os, itertools, math, collections, random, re, time, multiprocessing, multiprocessing.pool
from Bio import SeqIO
import enlighten
from bloom_filter import BloomFilter

def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Now add all the options to it
    parser.add_argument("fasta", type=str,
        help="name of the FASTA to read; must be indexable")
    parser.add_argument("-n", type=int, default=10000,
        help="the number of k-mers to count")
    parser.add_argument("-k", type=int, default=150,
        help="the length of each k-mer")
    parser.add_argument("--thread_count", type=int, default=multiprocessing.cpu_count(),
        help="number of k-mer counting threads to use")
    parser.add_argument("--batch_size", type=int, default=10000,
        help="number of forward-strand k-mer candidates to count in each batch")
    parser.add_argument("--bloom_error", type=int, default=1E-4,
        help="error rate on the bloom filter")
   
    return parser.parse_args(args)
    
# What characters are we allowed to have in k-mers?
ACGT = {'A', 'C', 'G', 'T'}
    
def all_ACGT(sequence):
    """
    Return true if the given sequence is all A, C, G, and T, and false otherwise.
    """
    
    for c in sequence:
        if c not in ACGT:
            return False
    return True
    
    
def count_kmers(k, acceptable, sequence):
    """
    Count all the kmers that are in the given bloom filter and found in the
    given sequence. Returns a Counter of the kmers that match the filter, and
    the total forward-strand kmers looked at.
    """
   
    counts = collections.Counter()
    total = 0
    
    if len(sequence) >= k:
        # There is anything to count.
        
        # Pre upper case everything
        sequence = sequence.upper()
        # And pre reverse complement
        rc_sequence = sequence.reverse_complement()
       
        for i in range(len(sequence) - (k - 1)):
            # Pull out every candidate kmer, in upper case
            candidate = sequence[i:i + k]
            
            total += 1
            
            if candidate in acceptable:
                if all_ACGT(candidate):
                    # Count it
                    counts[candidate] += 1
                
        for i in range(k, len(sequence) + 1):
            rc_candidate = rc_sequence[i - k:i]
            if rc_candidate in acceptable:
                if all_ACGT(rc_candidate):
                    # Count it
                    counts[rc_candidate] += 1
    
    return counts, total

def main(args):
    """
    Parses command line arguments, and does the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    if options.k == 0:
        sys.stderr.write("Cannot use empty k-mers\n")
        sys.exit(1)
        
    if options.n == 0:
        sys.stderr.write("Cannot use no k-mers\n")
        sys.exit(1)
        
    sys.stderr.write('Load FASTA...\n')
    
    # We access the same strings so many times that if we index here we use many times the genome's size in memory.
    # There may be non-cached string views involved.
    # Just load the whole FASTA.
    index = SeqIO.to_dict(SeqIO.parse(options.fasta, "fasta"))
    
    sys.stderr.write('Analyze sequence sizes...\n')
    
    # Compute lengths of all sequences, for sampling
    sequence_names = list(index.keys())
    sequence_lengths = [len(index[name]) for name in sequence_names]
    
    # Weight sequences by how many full-length kmers fit in them
    sequence_weights = [max(l - (options.k - 1), 0) for l in sequence_lengths]
    # And get ready to look up how many start positions we have available in each sequence
    available_starts_by_name = dict(zip(sequence_names, sequence_weights))
    
    if sum(sequence_weights) == 0:
        # We can't sample anything actually
        sys.stderr.write("No long enough sequences in file\n")
        sys.exit(1)
        
    sys.stderr.write('Sample {}-mers...\n'.format(options.k))
        
    with enlighten.Counter(total=options.n, desc='Sample', unit='{}-mers'.format(options.k)) as bar:
    
        # This will be all the k-mers we want to count.
        # We use a counter because if we sample the same k-mer multiple times at the
        # sampling stage we want to count it multiple times for mappability
        # assessment.
        kmers = collections.Counter()
        kmers_sampled = 0
        while kmers_sampled < options.n:
            # Sample a k-mer
            sequence_name = random.choices(sequence_names, sequence_weights)[0]
            kmer_start = random.randint(0, available_starts_by_name[sequence_name])
            kmer = index[sequence_name].seq[kmer_start:kmer_start + options.k]
            
            # Convert to upper case
            kmer = kmer.upper()
            
            if not all_ACGT(kmer):
                # Reject this one for having unacceptable letters
                continue
            
            strand = random.choice([False, True])
            if strand:
                # Flip half the kmers to the reverse strand
                kmer = kmer.reverse_complement()
            
            # Note that we are looking for this k-mer
            kmers[kmer] += 1
            # And that we sampled one.
            kmers_sampled += 1
            bar.update()
            
    # Make the bloom filter
    acceptable = BloomFilter(max_elements=options.n, error_rate=options.bloom_error)
    for kmer in kmers.keys():
        acceptable.add(kmer)
            
    sys.stderr.write('Count {}-mers...\n'.format(options.k))
            
    # Now traverse the whole FASTA and count
    counts = collections.Counter()
    
    with enlighten.Counter(total=sum(sequence_weights), desc='Count', unit='{}-mers'.format(options.k)) as bar:
    
        # We will do the counting in processes
        processes = multiprocessing.pool.Pool(options.thread_count)
        
        # We put the AsyncResults in this queue and handle them as they become ready.
        # A straggler could make it get a bit big.
        result_queue = collections.deque()
        
        def handle_result(result):
            """
            Process the return value from a counting job. Runs in main thread.
            """
            
            # See if something is done
            reply_counts, reply_kmers_processed = result
            # If so, handle it
            for kmer, count in reply_counts.items():
                if kmer in kmers:
                    # Not a bloom filter false positive
                    counts[kmer] += count
            bar.update(reply_kmers_processed)
        
        for sequence in index.values():
            # Where is the past-end for the k-mer start positions?
            start_past_end = max(len(sequence) - (options.k - 1), 0)
        
            # Where will the next batch start in this sequence
            cursor = 0
            
            while cursor < start_past_end:
                # Work out how big a batch to make
                this_batch_size = min(options.batch_size, start_past_end - cursor)
                
                # Find the piece of sequence to process.
                # Make sure to provide the tail end where k-mer starts aren't.
                part = sequence.seq[cursor:cursor + this_batch_size + (options.k - 1)]
                
                async_result = processes.apply_async(count_kmers, (options.k, acceptable, part))
                
                result_queue.append(async_result)
                
                cursor += this_batch_size
                
                while len(result_queue) > 0 and result_queue[0].ready():
                    # Pop off tasks that are done
                    handle_result(result_queue[0].get())
                    result_queue.popleft()
                    
                while len(result_queue) > options.thread_count * 4:
                    # Too many things waiting. Wait and pick some up before continuing.
                    handle_result(result_queue[0].get())
                    result_queue.popleft()
            
        while len(result_queue) > 0:
            # Collect all the jobs left at the end
            handle_result(result_queue[0].get())
            result_queue.popleft()
                
    # Bucket k-mers by multiplicity in the genome.
    # We know they all appear at least once.
    kmers_by_multiplicity = collections.defaultdict(list)
    for kmer, multiplicity in counts.items():
        kmers_by_multiplicity[multiplicity].append(kmer)
        
    # Count up the total number of kmers with each multiplicity, properly
    # weighting multiple sampling
    kmers_with_multiplicity = {m: sum((kmers[k] for k in ks)) for m, ks in kmers_by_multiplicity.items()}
    
    # Compute mappability
    effective_mapped = sum((count / m for m, count in kmers_with_multiplicity.items()))
    possible_mapped = sum(kmers.values())
    
    print("Expect to map {:.2f} {}-mers of out of {} total, for mappability of {:.2f}%".format(
        effective_mapped, options.k, options.n, effective_mapped / possible_mapped * 100))
        
    for m in sorted(kmers_with_multiplicity.keys()):
        print("{} copies: \t{} sampled {}-mers".format(m, kmers_with_multiplicity[m], options.k))
    
    return 0

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
