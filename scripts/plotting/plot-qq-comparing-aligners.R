#!/usr/bin/env Rscript

# plot-qq.R <stats TSV> <destination image file> [<comma-separated "aligner" names to include> [title]]

list.of.packages <- c("tidyverse", "ggrepel", "svglite")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require("tidyverse")
require("ggrepel")

# Read in the combined toil-vg stats.tsv, listing:
# correct, mapq, aligner (really graph name), read name, count
dat <- read.table(commandArgs(TRUE)[1], header=T)

if (! ("count" %in% names(dat))) {
    # If the count column is not present, add i
    dat$count <- rep(1, nrow(dat))
}

if (length(commandArgs(TRUE)) > 2) {
    # A set of aligners to plot is specified. Parse it.
    aligner.set <- unlist(strsplit(commandArgs(TRUE)[3], ","))
    # Subset the data to those aligners
    dat <- dat[dat$aligner %in% aligner.set,]
    # And restrict the aligner factor levels to just the ones in the set
    dat$aligner <- factor(dat$aligner, levels=aligner.set)
}

# Determine title
title <- ''
if (length(commandArgs(TRUE)) > 3) {
    title <- commandArgs(TRUE)[4]
}

# Determine the order of aligners, based on sorting in a dash-separated tag aware manner
aligner.names <- levels(factor(dat$aligner))
name.lists <- aligner.names %>% (function(name) map(name,  (function(x) as.list(unlist(strsplit(x, "-"))))))
# Transpose name fragments into a list of vectors for each position, with NAs when tag lists end early
max.parts <- max(sapply(name.lists, length))
name.cols <- list()
for (i in 1:max.parts) {
    name.cols[[i]] <- sapply(name.lists, function(x) if (length(x) >= i) { x[[i]] } else { NA })
}
name.order <- do.call(order,name.cols)
aligner.names <- aligner.names[name.order]
dat$aligner <- factor(dat$aligner, levels=aligner.names)
name.lists <- name.lists[name.order]

# Determine colors for aligners
light.colors <- c("#aaaa00", "#eedd88", "#44bb99", "#99ddff", "#77aadd", "#bbcc33", "#ee8866", "#ffaabb", "#dddddd", "#1f78b4","#e31a1c","#33a02c","#6600cc","#ff8000","#5c415d"    ,"#458b74","#698b22","#008b8b")
bold.colors <- light.colors # Use the same colors for everything
# We have to go through both lists together when assigning colors, because pe and non-pe versions of a condition need corresponding colors.
cursor <- 1
# This will map from non-pe condition name string to color index.
colors <- c()
# Given a string aligner condition name without "-pe", assign colors
assign_color <- function(name) {
    if (! name %in% names(colors)) {
        # No colors assigned for this pair of conditions, so assign them.
        
        if (cursor > length(bold.colors)) {
            write(colors, stderr())
            write(aligner.names, stderr())
            stop('Ran out of colors! Too many conditions!')
        }
        
        # We always assign pe and non-pe colors in lockstep, whichever we see first.
        # We need two entries for -se and no tag which are the same.
        new.colors <- c(bold.colors[cursor], light.colors[cursor], light.colors[cursor])
        names(new.colors) <- c(paste(name, 'pe', sep='-'), paste(name, 'se', sep='-'), name)
        
        # Commit changes to enclosing scope
        colors <<- c(colors, new.colors)
        cursor <<- cursor + 1
    }
}

# We want to pre-assign colors to aligners for this experiment, if those aligners show up
# Must use humanized names here.
pre.assign <- c("Bowtie2", "BWA", "Giraffe", "GiraffeFast", "GiraffePrimary", "GraphAligner", "Hisat2", "Map", "Minimap2")
for (i in 1:length(pre.assign)) {
    assign_color(pre.assign[[i]])
}


# Now assign colors for all the conditions we actually have
for (i in 1:length(name.lists)) {
    # For each name
    name.parts <- unlist(name.lists[[i]])
    if (name.parts[length(name.parts)] == "pe") {
        # Drop the pe tag if present
        name.parts <- name.parts[-c(length(name.parts))]
    }
    if (name.parts[length(name.parts)] == "se") {
        # Drop the se tag if present
        name.parts <- name.parts[-c(length(name.parts))]
    }
    
    # Join up to a string again
    name <- paste(name.parts, collapse='-')
    
    assign_color(name)
}

# Make colors a vector in the same order as the actually-used aligner names
colors <- colors[aligner.names]

dat$bin <- cut(dat$mq, c(-Inf,seq(0,60,1),Inf))

x <- as.data.frame(summarize(group_by(dat, bin, aligner), N=n(), mapq=mean(mq), mapprob=mean(1-10^(-mapq/10)), observed=weighted.mean(correct, count)))

dat.plot <- ggplot(x, aes(1-mapprob+1e-9, 1-observed+1e-9, color=aligner, size=N, weight=N, label=round(mapq,2))) +
    scale_color_manual(values=colors, guide=guide_legend(title=NULL, ncol=2)) +
    scale_y_log10("measured error", limits=c(5e-7,2), breaks=c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0)) +
    scale_x_log10("error estimate", limits=c(5e-7,2), breaks=c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0)) +
    scale_size_continuous("number", guide=guide_legend(title=NULL, ncol=4)) +
    geom_point() +
    geom_smooth() +
    geom_abline(intercept=0, slope=1, linetype=2) +
    theme_bw()
    
if (title != '') {
    # And a title
    dat.plot + ggtitle(title)
}

filename <- commandArgs(TRUE)[2]
ggsave(filename, height=4, width=7)
