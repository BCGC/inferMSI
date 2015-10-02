# decaplex.R
# Infer MSI status for a tumor sample with decaplex model
# Randy Johnson
# CCR Collaborative Bioinformatics Resource at Frederick National Laboratory
# Leidos Biomedical Research, Inc


######### Setup #########

library(inferMSI)
library(magrittr)

# script arguments
args <- commandArgs(TRUE)

# for testing
if(FALSE)
{
    args <- "i=tmp/TI025842000_S8_L001.realigned.bam.repeatseq o=MSIresults.out" # do one file
    args <- "wd=tmp" # do all repeatseq files in this directory
}

# parse arguments (if any were given)
if(length(args) > 0)
{
    tmp <- strsplit(args, ' ')[[1]] %>%
           strsplit(split = '=')


    args <- lapply(tmp, `[`, 2)
    names(args) <- sapply(tmp, `[`, 1)
}else{
    args <- list()
}

# default arguements
if(is.null(args$i))
    args$i <- system(paste('ls', args$wd, '| grep repeatseq$'), intern = TRUE)

if(length(args$i) == 0)
    stop("Missing input files")

if(is.null(args$o))
    args$o <- 'MSIresults.out'

# load model
data(euroDecaplex)


######### Infer status #########

out <- inferMSI(args$i, euroDecaplex, 0.2)

write.table(out, file = args$o, sep = '\t', quote = FALSE, )
