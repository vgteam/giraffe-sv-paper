# Giraffe Paper Products

This directory contains reusable work products for the manuscript *Pangenomics enables genotyping common structural variants in 5,202 diverse genomes*. It is organized as follows:

* `HGSVC_*`: Minimal graph indexes needed for running Giraffe against the HGSVC structural variant graph, with 6 embedded haplotypes and a path cover of graph regions without haplotypes. Giraffe will be able to relatively quickly regenerate the `.min` and `.gg` index files at runtime.
    * `HGSVC_hs38d1.full.gbwt`: Haplotype index
    * `HGSVC_hs38d1.dist`: Distance index
    * `HGSVC_hs38d1.xg`: Graph in xg format



