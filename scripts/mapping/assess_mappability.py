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

import argparse, sys, os, itertools, math, collections, random, re, threading, queue, time
from Bio import SeqIO
import enlighten

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
    parser.add_argument("--thread_count", type=int, default=10,
        help="number of k-mer counting threads to use")
   
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
    
    
def count_kmers(options, kmers, counts, task_queue, reply_queue):
    """
    Count all the kmers that are in kmers and found in the (sequence, first
    start, past-last start) batches available from task_queue. Puts counts into
    counter in counts. Sends back numbers of completed items via reply_queue.
    """
    
    while True:
        batch = task_queue.get()
        if batch is None:
            # Done!
            return
            
        (sequence, start, past_end) = batch
        
        for i in range(start, past_end):
    
            # Pull out every candidate kmer, in upper case
            candidate = sequence.seq[i:i + options.k].upper()
            
            if not all_ACGT(candidate):
                # Reject this one for having unacceptable letters
                continue
                
            if candidate in kmers:
                # Count it
                counts[candidate] += 1
                
            rc_candidate = candidate.reverse_complement()
            if rc_candidate in kmers:
                # And its reverse strand
                counts[rc_candidate] += 1
                
        reply_queue.put(past_end - start)
        task_queue.task_done()

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
            
    sys.stderr.write('Count {}-mers...\n'.format(options.k))
            
    # Now traverse the whole FASTA and count
    counts = collections.Counter()
    
    # We put regions to do into this queue, as tuples of (sequence object, first start, last start)
    task_queue = queue.Queue(1000)
    # We have an unlimited reply queue of processed position counts, so we
    # don't deadlock and one thread can handle filling one and emptying the
    # other.
    reply_queue = queue.Queue()
    BATCH_SIZE=10000
        
    with enlighten.Counter(total=sum(sequence_weights), desc='Count', unit='{}-mers'.format(options.k)) as bar:
        
        # Start threads to run the tasks
        threads = []
        for i in range(options.thread_count):
            threads.append(threading.Thread(target=count_kmers, args=(options, kmers, counts, task_queue, reply_queue)))
            threads[-1].start()
            
        for sequence in index.values():
            for i in range(0, max(len(sequence) - (options.k - 1), 0), BATCH_SIZE):
                # Compose the batch
                batch = (sequence, i, min(i + BATCH_SIZE, len(sequence) - (options.k - 1)))
                
                while True:
                    try:
                        # Send it off
                        task_queue.put_nowait(batch)
                        # Go get a new one
                        break
                    except queue.Full:
                        # Check for replies and then try to queue again
                        try:
                            while True:
                                # Update the progress bar with anything completed
                                number_completed = reply_queue.get_nowait()
                                bar.update(number_completed)
                                reply_queue.task_done()
                        except queue.Empty:
                            # Until we run out of completions here now.
                            
                            # Then there's probably nothing to do.
                            time.sleep(0)
                    
                
        for t in threads:
            # Add Nones to tell the threads to stop
            task_queue.put(None)
          
        any_alive = True
        while any_alive:
            # Until all the threads stop
            
            # Yield
            time.sleep(0)
            
            try:
                while True:
                    # Update the progress bar with anything completed.
                    # It should all be completed now that the threads are done.
                    number_completed = reply_queue.get_nowait()
                    bar.update(number_completed)
            except queue.Empty:
                # Until we run out of completions here now.
                pass
                
            # Poll the other threads
            any_alive = False
            for t in threads:
                if t.is_alive():
                    any_alive = True
                    break
            
            
        for t in threads:
            # Join all the threads
            t.join()
            
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
