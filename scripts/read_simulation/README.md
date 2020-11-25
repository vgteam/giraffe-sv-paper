## Real reads used for training the error model

Down-sampled to 1M, enough for training the error model

- novaseq6000-ERR3239454-shuffled-1m.fq.gz
- hiseq2500-ERR309934-shuffled-1m.fq.gz
- hiseqxten-SRR6691663-shuffled-1m.fq.gz
    
## Read simulation

For each sequencing techonology (see above), we make three files:

- `sim.fq.gz` FASTQ file with simulated reads
- `sim.gam` truth: simulated reads annotated with path in the graph
- `true.pos` truth: reformatted information about the true positions of simulate reads


```sh
for SEQ in novaseq6000 hiseq2500 hiseqxten
do
    ## simulate 1M reads from 1KGP graph for sample NA19239
    NREADS=1000000 SEQ=$SEQ sh simulate-1kgp.sh
    ## simulate 1M reads from HGSVC graph for sample NA19240
    NREADS=1000000 SEQ=$SEQ sh simulate-hgsvc.sh
done
```
