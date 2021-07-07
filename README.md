# DNAmethylation
## Code repository for Murgas et al. 2021, "A Bayesian Hierarchical Model for DNA Methylation in Colorectal Cancer" (in preparation)

### Overview
This study takes DNA methylation data, in the form of Infinium EPIC (850K) methylation array processed beta-values, from a cohort of 21 colorectal tumor patients. The multiple-sampling procedure provides two tumor samples for each patient, and a matching normal sample in some patients. These data then serve as inputs to a hierarchical random-effects model, which is fit to the logit-transformed beta-values (M-values) at each CpG site independently, using a Bayesian MCMC algorithm (Stan), including TCGA-based weakly informative prior probabilities. The resulting model fits are analyzed to determine sites with decreased variance within the tumor level compared to the normal tissue variance (i.e. preferential conservation within tumor). Subsequent analyses include significance bootstrapping and pathway analysis on significantly conserved genes.


### Tutorial for Analysis
1. **Download and prepare data**
    * The data used in this study is available from NCBI GEO GSE166212 [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166212)
    * Utilizing the Normalized_Beta_values xlsx file (containing processed methylation beta values for the CpG sites in the EPIC array across all samples in our data) is recommended
        * Alternatively, raw .IDAT data can be processed using minfi methods, see LoadDataAndQC.R script below. Note the data is saved as "myFA_bulkonly.Rdata" because the raw .IDAT files include some non-bulk tumor samples, and we only work with the bulk samples in our study.
    * Load the data into Rstudio using the script "PrepareDataFromGEO.R"
        * This script formats the columns of the Normalized_Beta_values file into a data.frame and saves as "myFA_bulkonly.Rdata" for model analysis
2. **(optional) Define priors in Stan model file**
    * We provide a Stan model file containing TCGA-based weakly informative priors: model_TCGApriors.stan
        * These priors were defined using an independent cohort from TCGA COAD project to empirically estimate distributions for each fixed effect parameter and random effect variance hyperparameter in the model. These estimates were then relaxed by maintaining the mode but increase the variance of each distribution, arriving at weakly-informative priors
    * Priors can be manually modified by editing the Stan file under the model{} section
3. **Run Hierarchical Model Fit with Stan**
    * **Note:** we recommend use of a high-performance compute cluster to analyze in batch jobs and using parallelization on multiple CPUs in order to process the data expediently. The scripts provided utilize this approach to split the 866,091 sites of the EPIC array into 87 chunks of 10,000 sites, each of which is processed in parallel. Using 24 cores per chunk, the entire dataset can be modeled in approximately 12-24 hours real time.
        * The shell script jobTCGA.sh provides an example of bash code to submit batch jobs using SLURM job scheduling. Parameters such as the partition and number of cores can be adjusted, and alternative scheduling software such as SGE will require slightly modified syntax.
    * Running the model fits on multiple sites is accomplished by the script StanCParallel.R, which takes an integer input to specify the chunk for which to offset and take 10,000 sites from the full data set to model.
        * **Note:** the patient ID and normal/tumor tissue status of each sample is determined based on the column names, be sure to modify this if you are working with different data, in order to specify the patient and tissue labels to the Stan algorithm
        * In addition to estimating the fixed effect parameters and variance hyperparameters of the model at each CpG site, a relative conservation score is computed using the posterior distributions of sigmaP and sigmaT
        * After fitting each model, the code will save a file (StanCParResults_X.Rdata, with X being the chunk index) to a specified output folder. The results files contain model parameter estimates as well as convergence diagnostics.
    * For the purposes of running the model on a small number of sites, this can be done on a local computer using the script GetPosteriorsSingleSites.R
        * This allows for examining individual CpG site model fits, for example to assess fits compared to the data or convergence
4. **Combine Model Fit Results**
    * After running the model fits in batches, the results for each chunk can be combined into a single results file using the loadSCRuns.R script
        * This script simply takes each results file and concatenates into a single data frame and saves as "FullResultsTCGA.Rdata"
        * Be sure to set the working directory in the first few lines of the script to the correct destination in your computer
5. **Analysis of Model Results**
    * Given the results files, one can explore different aspects of the data. The scripts StanAnalysis.Rmd and Paper_Figures.Rmd are two examples of the analysis we apply in our study
        * These analyses include:
            * Examining the overall distribution of parameter and hyperparameter estimates, as well as convergence diagnostics (we implement a cutoff filter based on the Rhat parameter for each CpG site, where Rhat>1.1 would indicate a poor model convergence, which typically occurs in <1% of CpG sites in our data)
            * Comparison of posterior estimates to the priors provided to the model
            * Examining the distribution of the relative conservation score across all CpG sites, and grouping CpG sites by functional regions and CpG-island relationship
            * Examining the spatial distribution of relative conservation score with respect to gene start site (defined by Ensembl, and correcting for gene direction)
            * Calculation of gene conservation score (averaging all CpG sites in each gene) and bootstrapping to establishing signficance of conservation
            * Examination of gene average conservation compared to known colorectal cancer gene lists
    * Additional analysis was performed using the Reactome Pathway Analysis online tool [link](https://reactome.org/PathwayBrowser/#/TOOL=AT), inputting a list of significantly conserved genes
        * Reactome utilizes a hypergeometric test to assess pathway over-representation in the list of genes
    * We welcome other researchers to consider, suggest, and perform further analyses on the model fit data, which include fixed effect parameters that could be used for differential methylation analysis for example. Additionally we welcome using our analysis pipeline for new data sets and potentially modifying our hierarchical model to accomodate different data structures with different hierarchical levels.


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


