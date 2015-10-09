# read.repeatseq.R
# read repeatseq data from file <f>
# Randy Johnson
# CCR Collaborative Bioinformatics Resrouce at Frederick National Laboratory
# Leidos Biomedical Research, Inc

read.repeatseq <- function(f, markers = NULL)
{
    if(!file.exists(f))
        stop(paste("File", f, "seems not to exist"))

    # read in main lines of info
    tmp <- readLines(f)
    tmp <- gsub('~', '', tmp[which(substr(tmp, 1, 1) == '~')])

    # break out different parts:
    ### site info
    tmp <- strsplit(tmp, ' REF:')
    siteInfo <- sapply(tmp, `[`, 1)

    ### reference allele
    tmp <- strsplit(sapply(tmp, `[`, 2), ' A:')
    ref <- sapply(tmp, `[`, 1)

    ### observed alleles
    tmp <- strsplit(sapply(tmp, `[`, 2), ' C:')
    obs <- sapply(tmp, `[`, 1)
    alleles <- lapply(strsplit(sapply(tmp, `[`, 1), '] ', fixed = TRUE), function(x)
                      {
                          # split observations apart
                          x <- strsplit(gsub(']', '', x, fixed = TRUE), '[', fixed = TRUE)
                          retval <- suppressWarnings(as.numeric(sapply(x, `[`, 2)))
                          names(retval) <- sapply(x, `[`, 1)
                          # order by allele size
                          retval <- retval[order(names(retval))]
                          # convert to proportion
                          retval <- retval / sum(retval)
                      })

    ### concordance
    tmp <- strsplit(sapply(tmp, `[`, 2), ' D:')
    concordance <- suppressWarnings(as.numeric(sapply(tmp, `[`, 1)))

    ### total reads
    tmp <- strsplit(sapply(tmp, `[`, 2), ' R:')
    totalReads <- suppressWarnings(as.numeric(sapply(tmp, `[`, 1)))

    ### total reads after filtering
    tmp <- strsplit(sapply(tmp, `[`, 2), ' S:')
    nReads <- suppressWarnings(as.numeric(sapply(tmp, `[`, 1)))

    ### total reads with no CIGAR sequence present
    tmp <- strsplit(sapply(tmp, `[`, 2), ' M:')
    noCIGAR <- suppressWarnings(as.numeric(sapply(tmp, `[`, 1)))

    ### average mapping quality
    tmp <- strsplit(sapply(tmp, `[`, 2), ' GT:')
    qc <- suppressWarnings(as.numeric(sapply(tmp, `[`, 1)))

    ### genotype
    tmp <- strsplit(sapply(tmp, `[`, 2), ' L:')
    geno <- sapply(tmp, `[`, 1)

    ### likelihodd of genotype
    lik <- suppressWarnings(as.numeric(sapply(tmp, `[`, 2)))


    ##### put everything into a list #####
    retval <- list()
    for(j in 1:length(siteInfo))
    {
        retval[[j]] <- list(siteInfo = siteInfo[j],
                            ref = ref[j],
                            alleles = alleles[[j]],
                            concordance = concordance[j],
                            totalReads = totalReads[j],
                            nReads = nReads[j],
                            noCIGAR = noCIGAR[j],
                            qc = qc[j],
                            geno = geno[j],
                            lik = lik[j])
    }

    if(!is.null(markers))
    {
        retvalNames <- unlist(sapply(retval, `[`, 'siteInfo'))

        if(length(which(retvalNames %in% names(markers))) == 0)
        {
            if(length(which(gsub('^chr', '', retvalNames) %in% names(markers))) == 0)
            {
                if(length(which(retvalNames %in% gsub('^chr', '', names(markers)))) == 0)
                {
                    warning('Marker names and repeatseq names do not match')
                    names(retval) <- retvalNames
                }else{
                    names(retval) <- markers[paste('chr', retvalNames, sep = '')]
                }
            }else{
                names(retval) <- markers[gsub('^chr', '', retvalNames)]
            }
        }else{
            names(retval) <- markers[retvalNames]
        }
    }else{
        names(retval) <- unlist(sapply(retval, `[`, 'siteInfo'))
    }

    return(retval)
}
