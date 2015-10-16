#!Rscript
# paired.R
# Run msisensor on a bunch of files
# Randall Johnson
# BSP CCR Genetics Core at Frederick Natonal Library
# Leidos Biomedical Research, Inc
# Created May 19, 2014
# Last Modified September 5, 2014


### args ###
# 1st argument: file pairs...three columns - normal bam file, tumor bam file, output file name
# 2nd argument: final output file name

args <- commandArgs(TRUE)

# read in file pairs
pairs <- read.table(args[1], sep = '\t', header = TRUE)

# run analyses
commands <- paste('msisensor msi -d microsats.sub.list -n', pairs$normal,
                  '-t', pairs$tumor, '-o', pairs$out, '-b 12')

# fetch critical information from the main file (file names are stored in pairs$out)
out <- data.frame(sample = as.character(pairs$out),
                  sites = NA,
                  somatic = NA,
                  pct = NA,
                  stringsAsFactors = FALSE)

for(i in 1:length(out$sample))
{
    # run MSIsensor
    system(commands[i])

    # collect data
    tmp <- unlist(read.table(out$sample[i], header = TRUE))
    out$sites[i] <- tmp[1]
    out$somatic[i] <- tmp[2]
    out$pct[i] <- tmp[3]

    # remove extraneous files
    system(paste('rm ', out$sample[i], '*', sep = ''))
}

# write output
write.table(out, file = args[2], row.names = FALSE, quote = FALSE)
