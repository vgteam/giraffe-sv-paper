# Plotting Scripts

These scripts create roc plots, qq plots, and a bar plot comparing speed and runtime. 
They are meant to be run on stats files produced by the speed and accuracy scripts in the   mapping directory.

## Running the scripts


### ROC and QQ Plots

Each of the `plot-roc*` and `plot-qq*` R scripts are run on a tsv file that is formatted:

```
correct mq  score   aligner
1   60  0   Giraffe
...
```

Each line in the stats file represents one read. The line in the example represents a read that was mapped by Giraffe correctly with mapping quality 60 and score 0.

The roc and qq plotting scripts are run on the stats file as:

```
./plot-roc.R stats.tsv plot.svg
```

The `plot-rocs.sh` script will run `plot-roc.R` using the stats filed generated from the `_accuracy` scripts in the `mapping` directory. It will produce separate roc plots for each sequencing technology for single and paired end reads.

### Speed, Runtime, and Memory Plots
Bar plots of speed, runtime, and memory can be run using `plot_real_speed.py`. This script requires tsv files formatted: 
```
graph   algorithm   reads   pairing speed   total_wall_clock_time   max_resident_set_size(kbytes)
```
The script is run:
```
python plot_real_speed.py [graph] [reads] [stat]
```
And requires a stat file named `speed_report_[reads].tsv`
