# Methods summary and general notes for the inferMSI R package

## Installation
To install this R package type the commands listed below into a terminal window, noting that a few other software packages may/will be needed (i.e. [R](https://cran.r-project.org), [MSIsensor](https://github.com/ding-lab/msisensor), [repeatseq](https://github.com/adaptivegenome/repeatseq) and [Snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home)). Installation of some additional R packages is also required prior to installation of this package, including: diptest, mclust, boot, missForest and magrittr.

```
git clone https://github.com/johnsonra/inferMSI/
cd inferMSI
R CMD INSTALL inferMSI
```

## Input files
All input files for inferMSI functions are assumed to be output from [repeatseq](https://github.com/adaptivegenome/repeatseq), using the "-repeatseq" argument. This format lists a summary of the length all mapped reads, giving inferMSI a measure of the distribution of microsatellite variability in the sample. For MSS samples this should be very similar to what would be returned in the default output (i.e. a genotype with one or two germline variants), but MSI samples, by nature of their instability, can be expected to have a more heterogeneous mixture of both somatic and germline variants.

## Developing a prediction model using training data
When developing a prediction model using trainMSI(), a list of germline input files and a list of somatic input files with known MSI status is analyzed:
- The distribution of germline variants in the population under study is characterized using the mclust package.
- Absolute deviations from the population distribution of germline variants are calculated for all markers in each somatic sample (missing scores for markers with insufficient read depth are imputed using the missForest() function).
- A binary classification prediction model is then generated using logistic regression to model known MSI status as a function of individual marker scores from each sample. The final model is an average of a series of prediction models generated using a leave-one-out validation scheme. The ability of each submodel to predict MSI status of the excluded sample is used to validate the sensitivity and specificity of the final model.

Additional details for trainMSI(), including specific function arguments and returned values can be accessed in R (after loading the inferMSI package) with the following command:

```
?trainMSI
```

## Inference of MSI status
The MSI status of additional samples using inferMSI() is inferred using the model obtained from trainMSI().
- Scores are obtained in the same manner, comparing the observed population distribution with the distribution of variants observed in each somatic sample, which may not collapse nicely to a standard genotype for markers that have become unstable.
- The scores are then combined with the beta values in the prediction model to obtain a predicted likelihood of microsatellite instability. Each value returned from inferMSI() can be interpreted as the probability (or a measure of this probability at the least) that the sample has significant instabilities. Models with noisy training data will suffer from a lack of sensitivity, resulting in lower over all MSI scores.

Additional details for inferMSI(), including specific function arguments and returned values can be accessed in R (after loading the inferMSI package) with the following command:

```
?inferMSI
```
