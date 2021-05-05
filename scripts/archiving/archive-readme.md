# Giraffe Paper Data

These data files were used in preparing the manuscript for *Genotyping common, large structural variations in 5,202 genomes using pangenomes, the Giraffe mapper, and the vg toolkit*, and will be useful for anyone looking to repeat or build on the experiments.

<!-- The archiving pipeline will turn links under the following URL into relative links: https://cgl.gi.ucsc.edu/data/giraffe/ -->

## Code
All files should be readable by vg v1.32.0, with commit SHA1 hash 095c529f8c70521ca60cf0435bac1b0b4ffd1f6d. The tool is available on Github at https://github.com/vgteam/vg/releases/tag/v1.32.0. A [fully static binary for GNU/Linux 3.2.0](https://cgl.gi.ucsc.edu/data/giraffe/code/vg/v1.32.0/vg) and a [gzip-compressed source tar file](https://cgl.gi.ucsc.edu/data/giraffe/code/vg/v1.32.0/vg-v1.32.0.tar.gz) are available.

The code to reproduce the analysis in the paper can be found [on the web](https://github.com/vgteam/giraffe-sv-paper) or [in the archive](software/code/giraffe-sv-paper/).

## Read Mapping
For the mapping experiments, we used a series of graphs.

For each graph, the graph itself is available in ".vg" and ".xg" files. The ".xg" files are all in XG format, while the ".vg" files are in either VPKG-encapsulated PackedGraph or VG Protobuf format, depending on the graph.

Graphs have other associated index files.

The ".snarls" and ".trivial.snarls" files contain information on the snarl decomposition of the graph, excluding or including "trivial" snarls with no internal nodes, respectively. These are not needed for mapping, but can be exported to JSON with "vg view -Rj file.snarls".

The ".gcsa" and ".gcsa.lcp" files contain the full-text indexes needed for mapping with VG-MAP.

The ".dist" files, required for mapping with Giraffe, contain the Giraffe distance indexes.

Graphs can also have one or more sets of ".gbwt", ".gg", and ".min" files, which contain the haplotype database, the graph node sequence information, and the minimizer lookup information, respectively. A set of matched files is required for mapping with Giraffe. The ".gg" file can be regenerated from the GBWT and the graph, and obviates the need to pass the full graph to Giraffe.

These matched sets come in several versions. The "raw" set, with no additional extension, contains the original haplotypes articulated in the input VCF (for VCF-based graphs) or the named paths of the input alignment (for alignment-based graphs). The "full" set starts with the raw graph and covers connected components without haplotype data with synthesized haplotypes as described in the manuscript. The "cover.<number>" and "sampled.<number>" sets use haplotypes produced to cover graph elements, or by sampling from the raw haplotypes with additional haplotypes for connected components with no data, as described in the manuscript. Each is built to have the included number’s worth of haplotype coverage overall. The "cover" sets, with no associated number, use haplotypes produced to cover graph elements up to an overall coverage of 16. Of these sets, the "sampled.64" set is recommended for general use.

The available graphs are as follows:
* [The human 1000 Genomes Project Liftover No Seg Dupes (1000GPlons) graph](https://cgl.gi.ucsc.edu/data/giraffe/mapping/graphs/for-NA19239/1000gplons/hs38d1/)
    * With simulation samples removed, as used for mapping, in the files named "`1000GPlons_hs38d1_filter.*`".
    * Without those samples removed, as used for simulation, in the files named "`1kg_hs37d5.*`".
* [The human HGSVC graph](https://cgl.gi.ucsc.edu/data/giraffe/mapping/graphs/for-NA19240/hgsvc/hs38d1/), in the files named "`HGSVC_hs38d1.*`". The simulation samples were not removed from this graph.
* [The linear human reference version 38 graph (a negative control for the HGSVC graph)](https://cgl.gi.ucsc.edu/data/giraffe/mapping/graphs/generic/primary/hs38d1/), in the files named "`primaryhs38d1.*`".
* [The yeast graph used as a mapping target](https://cgl.gi.ucsc.edu/data/giraffe/mapping/graphs/generic/cactus/yeast_subset/), in the files named "`yeast_subset.*`".
* [The graph used to simulate yeast reads](https://cgl.gi.ucsc.edu/data/giraffe/mapping/graphs/generic/cactus/yeast_all/), in the files named "`yeast_all.*`".
* [The linear yeast *Saccharomyces cerevisiae* S288c reference graph (a negative control for the yeast graph)](https://cgl.gi.ucsc.edu/data/giraffe/mapping/graphs/generic/primary/S288C/), in the files named "`primaryS288C.*`".

Additionally, the simulated reads used in the mapping experiments are available. These reads are available in vg’s Graph Alignment/Map (GAM) format, annotated with their true positions along named paths, and can be converted to interleaved FASTQ with the command `vg view -aiX file.gam`, or to JSON with the command `vg view -aj file.gam`. All files contain interleaved, paired reads.
* For the 1000GPlons graph, reads simulated to resemble sample NA19239 as sequenced by different Illumina sequencing technologies are available:
   * [For HiSeq2500](https://cgl.gi.ucsc.edu/data/giraffe/reads/sim/for-NA19239/1000gplons/hs38d1/hiseq2500/out_sim_gbwt/sim.gam)
   * [For HiSeq X Ten](https://cgl.gi.ucsc.edu/data/giraffe/reads/sim/1000gplons/hs38d1/hs37d5/hiseqxten/out_sim_gbwt/sim.gam)
   * [For NovaSeq 6000](https://cgl.gi.ucsc.edu/data/giraffe/reads/sim/1000gplons/hs38d1/hs37d5/novaseq6000/out_sim_gbwt/sim.gam).
* For the HGSVC graph, reads simulated to resemble sample NA19240 as sequenced by different Illumina sequencing technologies are available:
   * [For HiSeq2500](https://cgl.gi.ucsc.edu/data/giraffe/reads/sim/for-NA19240/hgsvc/grch38/hiseq2500/out_sim_gbwt/sim.gam)
   * [For HiSeq X Ten](https://cgl.gi.ucsc.edu/data/giraffe/reads/sim/for-NA19240/hgsvc/grch38/hiseqxten/out_sim_gbwt/sim.gam)
   * [For NovaSeq 6000](https://cgl.gi.ucsc.edu/data/giraffe/reads/sim/for-NA19240/hgsvc/grch38/novaseq6000/out_sim_gbwt/sim.gam)
* For the yeast graph, reads simulated to resemble different yeast strains as sequenced by a single Illumina sequencing technology are available:
   * [For strain DBVPG6044](https://cgl.gi.ucsc.edu/data/giraffe/reads/sim/yeast/sim-DBVPG6044.gam)
   * [For strain DBVPG6765]https://cgl.gi.ucsc.edu/data/giraffe/(reads/sim/yeast/sim-DBVPG6765.gam)
   * [For strain N44](https://cgl.gi.ucsc.edu/data/giraffe/reads/sim/yeast/sim-N44.gam)
   * [For strain UWOPS034614](https://cgl.gi.ucsc.edu/data/giraffe/reads/sim/yeast/sim-UWOPS034614.gam)
   * [For strain UWOPS919171](https://cgl.gi.ucsc.edu/data/giraffe/reads/sim/yeast/sim-UWOPS919171.gam)
   * [For strain Y12](https://cgl.gi.ucsc.edu/data/giraffe/reads/sim/yeast/sim-Y12.gam)
   * [For strain YPS138](https://cgl.gi.ucsc.edu/data/giraffe/reads/sim/yeast/sim-YPS138.gam)

## Structural variant genotyping

During the evaluation of the SV genotyping accuracy we reproduced analysis from Hickey et al. on the HGSVC graph and the GIAB graph (although now using GIAB v0.6). The indexes, including new Giraffe indexes, for these pangenomes are available:
   * [For the HGSVC graph](https://cgl.gi.ucsc.edu/data/giraffe/calling/hgsvc/)
   * [For the GIAB graph](https://cgl.gi.ucsc.edu/data/giraffe/calling/giab/)

There are also [indexes for the pangenome containing the combined SV catalogs](https://cgl.gi.ucsc.edu/data/giraffe/calling/combined-sv-graph/).

There is [a summary of the structural variants genotyped across 2,000 samples from the MESA cohort and their frequencies](https://cgl.gi.ucsc.edu/data/giraffe/calling/vggiraffe-sv-mesa-svsites.csv.gz)

The structural variants found in the high-coverage 1000 Genomes Project dataset are:
   * [SV-eQTLs using the subset of samples with RNA-seq information from the GEUVADIS dataset](calling/vggiraffe-sv-eqtl-geuvadis.FDR01.csv)
   * [SVs with strong inter-super-population frequency patterns (more than 10% shift in allele frequency compared to the median frequency across all super-populations)](https://cgl.gi.ucsc.edu/data/giraffe/calling/vggiraffe-sv-superpop-af-diff-med10.csv.gz)
   * [SVs overlapping coding, untranslated, promoter or introns of protein-coding genes](https://cgl.gi.ucsc.edu/data/giraffe/calling/vggiraffe-sv-2504kgp-pcgenes.csv.gz)
   * [Information about each SV site, including frequencies of the most expressed and second most expressed alleles, and frequency in the 5 super-populations](https://cgl.gi.ucsc.edu/data/giraffe/calling/vggiraffe-sv-2504kgp-svsites.csv.gz)
   * [The allele counts for each of the 2,504 unrelated samples across SV sites](https://cgl.gi.ucsc.edu/data/giraffe/calling/vggiraffe-sv-2504kgp-svsite-ac.tsv.gz)

