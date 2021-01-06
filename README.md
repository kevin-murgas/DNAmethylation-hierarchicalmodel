# DNAmethylation
## Code repository for Murgas et al. 2021, "A Bayesian Hierarchical Model for DNA Methylation in Colorectal Cancer" (in preparation)

### Overview
This study takes DNA methylation data, in the form of Infinium EPIC (850K) methylation array processed beta-values, from a cohort of 21 colorectal tumor patients. The multiple-sampling procedure provides two tumor samples for each patient, and a matching normal sample in some patients. These data then serve as inputs to a hierarchical random-effects model, which is fit to the logit-transformed beta-values (M-values) at each CpG site independently, using a Bayesian MCMC algorithm (Stan), including TCGA-based weakly informative prior probabilities. The resulting model fits are analyzed to determine sites with decreased variance within the tumor level compared to the normal tissue variance (i.e. preferential conservation within tumor). Subsequent analyses include significance bootstrapping and pathway analysis on significantly conserved genes.

### Code Files
Listed based on relevant order of processing

1. **Pre-processing:** 
    * LoadDataAndQC.R
        * This file is used to pre-process the raw EPIC array IDAT files into beta values for each sample
        * Processing methods are based on minfi R package 'noob' pre-process method and various QC approaches
        * Processed beta values are stored in a file myFA.Rdata, containing a data structure FullAnnotation, with each row corresponding to a CpG site including relevant data annotation for that site and subsequent beta values for each sample included in the cohort.
2. **TCGA-based weakly informative prior estimation:**
    * estimateTCGApriors.R
        * This file loads in a set of relevant TCGA data from the COAD cohort, specifically only patients with 450K methylation array data for at least 2 tumor samples or 1 tumor and 1 normal, so as to support empirical estimation of hierarchical variance parameters
        * Beta-values (pre-processed) are analyzed such that individual hierarchical mean and variance parameters can be estimated via statistical variance calculation, with an estimate for each component of the model at each CpG in the 450K array
        * The resulting estimates are fit with various statistical distributions as informative priors, and prior relaxation is performed by algebraically manipulating distribution parameters such that variance is increased while mode (peak) is maintained
    * prior_sensitivity.Rmd
        * This R-markdown file performs analysis on different sets of priors, after fitting the Stan model on 850K data (below, step 3)
        * Various comparisons of the prior distributions to the posterior median distributions are made to assess sensitivity of posterior results to the priors used
3. **Stan Hierarchical Model Fitting:**
    * model_TCGApriors.stan
        * This file contains Stan code for the model, defining inputs, priors, and the hierarchical model
    * StanCParallel.R
        * This file is used to parallelize Stan fits on blocks of 10000 CpG sites (in 87 chunks to cover all 866836 ~ 850K sites)
        * An "empty" Stan fit is produce by running with 0 iterations, which creates a model basis for Stan to efficiently fit the hierarchical model at all sites
        * For each CpG site in a given block, data is extracted as the set of beta values for all samples at that site along with patient labels and tumor indicator for use in the Stan algorithm
        * After each site is fit with hierarchical model, relevant statistics (mean, SEM, median, n_eff, Rhat) for each parameter and log-posterior of the fit are extracted and saved
    * jobTCGA.sh
        * This is an example job bash script for running the StanCParallel.R code as a batch processing job, based on SLURM scheduling syntax
    * GetPosteriorsSingleSites.R
        * This file essentially does the same as StanCParallel.R, but fits on a specified subset of sites (manually entered in code), so as to examine individual site fit results
4. **Concatenate Fit Results**
    * loadSCRuns.R
        * This script loads each of the 87 blocks of StanCParallel.R fit results and concatenates each parameter into a single data structure for all 850K sites, then saves these together in a file FullResultsTCGA.Rdata
5. **Analysis of Fit Results**
    * StanAnalysis.Rmd
        * This R-markdown file loads the full parameter results from the hierarchical model fits and presents a basic summary of the results
        * Importantly (for the next script to use), this script also prepares the data by processing CpG sites on the gene level and performing bootstrapping, which are saved as three files: prepped_data.csv, gened_data.csv, permed_data.csv
    * Paper_Figures.Rmd
        * This R-markdown file essentially performs all of the data summarization and statistical analysis presented in our paper.
        * Depends on the prepared data from StanAnalysis.Rmd

### Contact for questions
Any questions please direct to Kevin Murgas, the owner of this repository, on Github or via email at kevin.murgas@stonybrookmedicine.edu


