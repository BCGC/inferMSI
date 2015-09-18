# repeatseq.snakemake
# Randy Johnson
# This Snakefile will generate repeatseq files given a bam files
# invoke with snakemake -s repeatseq.snakemake
# CCR Collaborative Bioinformatics Resrouce at Frederick National Laboratory
# Leidos Biomedical Research, Inc

from glob import glob

# set these variables to the appropriate paths
reference = "~/hg19/human_g1k_v37.fasta"
regions = "../../Data/repeatSeq/b37.2014.sub_10.regions"

# run on all bam files in the current working directory
bams = glob("*.bam")
rptseq = expand("{b}.repeatseq", b=bams)

rule all:
    input: rptseq


rule repeatseq:
    input: "{indivs}.bam"
    output: "{indivs}.bam.repeatseq",
            temp("{indivs}.bam.vcf")
    shell: "repeatseq -repeatseq {input} {reference} {regions}"
