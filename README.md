sccomp - Outlier-aware and count-based compositional analysis of
single-cell data
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidyseurat/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/tidyseurat/actions/)
<!-- badges: end -->

# <img src="inst/logo-01.png" height="139px" width="120px" />

Cell omics such as single-cell genomics, proteomics and microbiomics
allow the characterisation of tissue and microbial community
composition, which can be compared between conditions to identify
biological drivers. This strategy has been critical to unveiling markers
of disease progression such as cancer and pathogen infection. For cell
omic data, no method for differential variability analysis exists, and
methods for differential composition analysis only take a few
fundamental data properties into account. Here we introduce sccomp, a
generalised method for differential composition and variability analyses
able to jointly model data count distribution, compositionality,
group-specific variability and proportion mean-variability association,
with awareness against outliers. Sccomp is an extensive analysis
framework that allows realistic data simulation and cross-study
knowledge transfer. Here, we demonstrate that mean-variability
association is ubiquitous across technologies showing the inadequacy of
the very popular Dirichlet-multinomial modelling and provide mandatory
principles for differential variability analysis. We show that sccomp
accurately fits experimental data, with a 50% incremental improvement
over state-of-the-art algorithms. Using sccomp, we identified novel
differential constraints and composition in the microenvironment of
primary breast cancer.

# Installation

## (simple) Suggested for single-cell and CyTOF analyses

**Bioconductor**

``` r
if (!requireNamespace("BiocManager")) {
   install.packages("BiocManager")
 }
 BiocManager::install("sccomp")
```

**Github**

``` r
devtools::install_github("stemangiola/sccomp")
```

## (more complex and efficient, until further optimisation of the default installation) Suggested for microbiomics

**Github**

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
check_cmdstan_toolchain()
install_cmdstan(cores = 2)
# Then, check the correct cmdstanr installation here
# https://mc-stan.org/cmdstanr/articles/cmdstanr.html

# Then install sccomp with the cmdstanr branch
devtools::install_github("stemangiola/sccomp@cmdstanr")
```

# Analysis

`sccomp` can model changes in composition and variability. Normally the
furmula for variability is either `~1`, which assumes that the
cell-group variability is independent on any covariate, or
`~ factor_of_interest`, which assumes that the model is dependent on the
factor of interest only. However, more complex models for variability
are possible, is the sample size is large. In any case the model for
variability must be a subset of the model for composition.

## From Seurat Object

``` r
res =
  seurat_obj |>
  sccomp_glm( 
   formula_composition = ~ type, 
    formula_variability = ~ 1, 
    sample, 
    cell_group 
  )
```

## From SingleCellExperiment Object

``` r
res =
  sce_obj |>
  sccomp_glm( 
    formula_composition = ~ type, 
    formula_variability = ~ 1, 
    sample, 
    cell_group 
  )
```

## From data.frame

``` r
res =
  seurat_obj[[]] |>
  sccomp_glm(
    formula_composition = ~ type, 
    formula_variability = ~ 1, 
    sample, 
    cell_group 
  )
```

## From counts

``` r
res =
  counts_obj |>
  sccomp_glm( 
    formula_composition = ~ type, 
    formula_variability = ~ 1, 
    .sample = sample,
    .cell_group = cell_group,
    .count = count
  )
```

    ## sccomp says: outlier identification first pass - step 1/3 [ETA: ~20s]

    ## Finished in  2.4 seconds.
    ## Running standalone generated quantities after 1 MCMC chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## sccomp says: outlier identification second pass - step 2/3 [ETA: ~60s]

    ## Finished in  9.0 seconds.
    ## Running standalone generated quantities after 1 MCMC chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## sccomp says: outlier-free model fitting - step 3/3 [ETA: ~20s]

<<<<<<< HEAD
    ## Running MCMC with 6 parallel chains...
    ## 
    ## Chain 1 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 1 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[299] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 1

    ## Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 1 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[299] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 1

    ## Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 1 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 1

    ## Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 1 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[6] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 1

    ## Chain 2 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 2 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[4] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 2

    ## Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 2 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[4] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 2

    ## Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 2 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[5] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 2

    ## Chain 3 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 3 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 3

    ## Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 3 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[10] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 3

    ## Chain 4 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 4 Exception: Exception: beta_binomial_lpmf: Second prior sample size parameter[328] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 4

    ## Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 4 Exception: Exception: beta_binomial_lpmf: Second prior sample size parameter[328] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 4

    ## Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 4 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 4

    ## Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 4 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[10] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 4

    ## Chain 5 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 5 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 5 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[26] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 5 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 5 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 5

    ## Chain 5 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 5 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[26] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 5 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 5 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 5

    ## Chain 5 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 5 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 5 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 5 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 5

    ## Chain 5 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 5 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[20] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 5 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 5 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 5

    ## Chain 6 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 6 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 6 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[296] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 6 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 6 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 6

    ## Chain 6 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 6 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[296] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 6 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 6 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 6

    ## Chain 6 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 6 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 6 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 6 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 6

    ## Chain 6 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 6 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[22] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 6 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 6 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 6

    ## Chain 1 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 2 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 6 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 3 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 4 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 1 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 5 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 2 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 6 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 1 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 1 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 5 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 3 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 4 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 2 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 2 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 6 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 6 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 1 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 5 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 5 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 3 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 3 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 4 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 4 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 6 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 1 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 2 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 3 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 5 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 4 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 6 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 3 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 5 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 1 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 6 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 2 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 3 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 4 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 1 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 5 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 6 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 3 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 1 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 4 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 5 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 6 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 2 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 3 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 1 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 4 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 5 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 6 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 1 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 1 finished in 9.1 seconds.
    ## Chain 3 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 5 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 6 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 6 finished in 9.0 seconds.
    ## Chain 3 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 4 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 3 finished in 9.5 seconds.
    ## Chain 2 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 5 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 5 finished in 9.5 seconds.
    ## Chain 4 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 2 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 4 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 4 finished in 10.6 seconds.
    ## Chain 2 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 2 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 2 finished in 12.9 seconds.
    ## 
    ## All 6 chains finished successfully.
    ## Mean chain execution time: 10.1 seconds.
    ## Total execution time: 13.3 seconds.
=======
    ## sccomp says: the composition design matrix has columns: (Intercept), typecancer

    ## sccomp says: the variability design matrix has columns: (Intercept)

``` r
res
```

    ## # A tibble: 72 × 14
    ##    cell_group parameter   covariate c_lower c_effect c_upper   c_pH0    c_FDR
    ##    <chr>      <chr>       <chr>       <dbl>    <dbl>   <dbl>   <dbl>    <dbl>
    ##  1 B1         (Intercept) <NA>       0.982     1.15   1.32   0       0       
    ##  2 B1         typecancer  type      -1.18     -0.867 -0.551  0       0       
    ##  3 B2         (Intercept) <NA>       0.390     0.693  0.982  0.00200 0.000119
    ##  4 B2         typecancer  type      -1.09     -0.625 -0.160  0.0340  0.00413 
    ##  5 B3         (Intercept) <NA>      -0.573    -0.291 -0.0134 0.253   0.0192  
    ##  6 B3         typecancer  type      -0.726    -0.304  0.0915 0.292   0.0758  
    ##  7 BM         (Intercept) <NA>      -1.22     -0.935 -0.652  0       0       
    ##  8 BM         typecancer  type      -0.722    -0.308  0.0898 0.291   0.0650  
    ##  9 CD4 1      (Intercept) <NA>       0.225     0.408  0.585  0.0140  0.000860
    ## 10 CD4 1      typecancer  type      -0.0757    0.176  0.428  0.580   0.164   
    ## # … with 62 more rows, and 6 more variables: v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, count_data <list>

## With contrasts

``` r
res =
  counts_obj |>
  sccomp_glm( 
    formula_composition = ~ 0 + type, 
    formula_variability = ~ 1, 
    contrasts =  c(typecancer - typebenign, typebenign - typecancer),
    .sample = sample,
    .cell_group = cell_group,
    .count = count
  )
```

    ## sccomp says: outlier identification first pass - step 1/3 [ETA: ~20s]

    ## sccomp says: outlier identification second pass - step 2/3 [ETA: ~60s]

    ## sccomp says: outlier-free model fitting - step 3/3 [ETA: ~20s]

    ## sccomp says: the composition design matrix has columns: typebenign, typecancer

    ## sccomp says: the variability design matrix has columns: (Intercept)
>>>>>>> master

``` r
res
```

<<<<<<< HEAD
    ## # A tibble: 72 × 8
    ##    cell_group parameter   c_lower c_effect c_upper  c_pH0   c_FDR count_data
    ##    <chr>      <chr>         <dbl>    <dbl>   <dbl>  <dbl>   <dbl> <list>    
    ##  1 B1         (Intercept)  0.550     0.707  0.873  0      0       <tibble>  
    ##  2 B1         typecancer  -1.19     -0.885 -0.583  0      0       <tibble>  
    ##  3 B2         (Intercept)  0.124     0.368  0.623  0.0866 0.00447 <tibble>  
    ##  4 B2         typecancer  -1.12     -0.621 -0.150  0.0408 0.00491 <tibble>  
    ##  5 B3         (Intercept) -0.603    -0.411 -0.206  0.0223 0.00154 <tibble>  
    ##  6 B3         typecancer  -0.583    -0.207  0.145  0.483  0.133   <tibble>  
    ##  7 BM         (Intercept) -1.31     -1.10  -0.881  0      0       <tibble>  
    ##  8 BM         typecancer  -0.751    -0.334  0.0684 0.253  0.0644  <tibble>  
    ##  9 CD4 1      (Intercept)  0.361     0.489  0.618  0      0       <tibble>  
    ## 10 CD4 1      typecancer  -0.0715    0.164  0.410  0.616  0.178   <tibble>  
    ## # … with 62 more rows
=======
    ## # A tibble: 72 × 9
    ##    cell_group parameter        covariate c_lower c_effect c_upper  c_pH0   c_FDR
    ##    <chr>      <chr>            <chr>       <dbl>    <dbl>   <dbl>  <dbl>   <dbl>
    ##  1 B1         typebenign - ty… type       0.615     0.945  1.31   0      0      
    ##  2 B1         typecancer - ty… type      -1.31     -0.945 -0.615  0      0      
    ##  3 B2         typebenign - ty… type       0.170     0.642  1.09   0.0340 0.00428
    ##  4 B2         typecancer - ty… type      -1.09     -0.642 -0.170  0.0340 0.00428
    ##  5 B3         typebenign - ty… type      -0.145     0.224  0.582  0.448  0.136  
    ##  6 B3         typecancer - ty… type      -0.582    -0.224  0.145  0.448  0.136  
    ##  7 BM         typebenign - ty… type      -0.0439    0.348  0.752  0.217  0.0540 
    ##  8 BM         typecancer - ty… type      -0.752    -0.348  0.0439 0.217  0.0540 
    ##  9 CD4 1      typebenign - ty… type      -0.396    -0.150  0.102  0.666  0.226  
    ## 10 CD4 1      typecancer - ty… type      -0.102     0.150  0.396  0.666  0.226  
    ## # … with 62 more rows, and 1 more variable: count_data <list>

Of the output table, the estimate columns startwith the prefix `c_`
indicate `composition`.

## Suggested settings for single-cell RNA sequencing

We reccommend to set `bimodal_mean_variability_association  = TRUE`. The
bimodality of the mean-variability association can be confirmed from the
plots$credible_intervals_2D (see below).

## Suggested settings for CyTOF and microbiome data

We reccommend to set `bimodal_mean_variability_association  = FALSE`
(Default).
>>>>>>> master

## Visualise data + inference

``` r
plots = plot_summary(res) 
```

<<<<<<< HEAD
    ## Running standalone generated quantities after 6 MCMC chains, 1 chain at a time ...
    ## 
    ## Chain 1 finished in 0.0 seconds.
    ## Chain 2 finished in 0.0 seconds.
    ## Chain 3 finished in 0.0 seconds.
    ## Chain 4 finished in 0.0 seconds.
    ## Chain 5 finished in 0.0 seconds.
    ## Chain 6 finished in 0.0 seconds.
    ## 
    ## All 6 chains finished successfully.
    ## Mean chain execution time: 0.0 seconds.
    ## Total execution time: 9.2 seconds.

=======
>>>>>>> master
    ## Joining, by = c("sample", "cell_group")
    ## Joining, by = c("cell_group", "type")

Plot of group proportion, faceted by groups. The blue boxplots represent
the posterior predictive check. If the model is likely be descriptively
adequate to the data, the blue boxplot should roughly overlay with the
black boxplot, which represent the observed data. The outliers are
coloured in red. A boxplot will be returned for every (discrete)
covariates present in `formula_composition`. The color coding represent
the significant associations for composition and/or variability.

``` r
plots$boxplot
```

    ## [[1]]

![](inst/figures/unnamed-chunk-12-1.png)<!-- -->

Plot of estimates of differential composition (c\_) on the x axis, and
differential variability (v\_) on the y axis. The error bars represent
95% credible intervals. The dashed lines represent the minimal effect
that the hypothesis test is based on. An effect is labelled as
significant if bigger than the minimal effect according to the 95%
credible interval. Facets represent the covariates in the model.

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-13-1.png)<!-- -->

## Visualisation of the MCMC chains from the posterior distribution

It is possible to directly evaluate the posterior distribution. In this
example we plot the Monte Carlo chain for the slope parameter of the
first cell type. We can see that has converged and is negative with
probability 1.

``` r
attr(res, "fit")$draws(variables = c("beta[2,1]")) %>% 
mcmc_trace()
```

![](inst/figures/unnamed-chunk-14-1.png)<!-- -->

## Differential variability

We can model the cell-group variability also dependent on type, and so
test differences in variability

``` r
res = 
  counts_obj |>
  sccomp_glm( 
    formula_composition = ~ type, 
    formula_variability = ~ type, 
    .sample = sample,
    .cell_group = cell_group,
    .count = count
  )
```

    ## sccomp says: outlier identification first pass - step 1/3 [ETA: ~20s]

    ## Finished in  2.5 seconds.
    ## Running standalone generated quantities after 1 MCMC chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## sccomp says: outlier identification second pass - step 2/3 [ETA: ~60s]

    ## Finished in  10.4 seconds.
    ## Running standalone generated quantities after 1 MCMC chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## sccomp says: outlier-free model fitting - step 3/3 [ETA: ~20s]

<<<<<<< HEAD
    ## Running MCMC with 6 parallel chains...
    ## 
    ## Chain 1 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 1 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 1

    ## Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 1 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[32] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 1

    ## Chain 2 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 2 Exception: Exception: beta_binomial_lpmf: Second prior sample size parameter[281] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 2

    ## Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 2 Exception: Exception: beta_binomial_lpmf: Second prior sample size parameter[281] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 2

    ## Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 2 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 2

    ## Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 2 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[3] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 2

    ## Chain 3 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 3 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[24] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 3

    ## Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 3 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[24] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 3

    ## Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 3 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 3

    ## Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 3 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[11] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 3

    ## Chain 4 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 4 Exception: Exception: beta_binomial_lpmf: Second prior sample size parameter[10] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 4

    ## Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 4 Exception: Exception: beta_binomial_lpmf: Second prior sample size parameter[10] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 4

    ## Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 4 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 4

    ## Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 4 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 4

    ## Chain 5 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 5 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 5 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 5 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 5 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 5

    ## Chain 5 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 5 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 5 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 5 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 5

    ## Chain 6 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 6 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 6 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[300] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 6 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 6 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 6

    ## Chain 6 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 6 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[300] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 6 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 6 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 6

    ## Chain 6 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 6 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 6 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 6 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 6

    ## Chain 6 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:

    ## Chain 6 Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[2] is 0, but must be positive finite! (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 46, column 0 to line 51, column 5) (in '/stornext/HPCScratch/mangiola.s_HPC_scratch/.Rtemp/Rtmp0tPgFn/model-631856d421f2.stan', line 169, column 4 to line 176, column 6)

    ## Chain 6 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,

    ## Chain 6 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

    ## Chain 6

    ## Chain 1 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 4 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 2 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 6 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 1 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 5 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 4 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 6 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 2 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 3 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 1 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 1 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 5 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 4 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 4 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 6 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 6 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 2 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 2 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 3 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 5 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 5 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 4 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 1 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 6 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 5 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 3 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 2 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 3 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 4 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 6 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 1 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 5 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 3 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 4 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 6 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 3 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 5 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 2 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 1 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 4 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 6 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 3 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 5 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 4 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 2 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 3 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 6 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 1 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 4 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 5 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 3 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 4 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 6 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 4 finished in 12.4 seconds.
    ## Chain 2 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 1 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 5 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 6 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 6 finished in 12.7 seconds.
    ## Chain 3 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 5 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 5 finished in 13.4 seconds.
    ## Chain 3 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 3 finished in 14.0 seconds.
    ## Chain 2 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 1 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 1 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 1 finished in 15.4 seconds.
    ## Chain 2 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 2 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 2 finished in 16.3 seconds.
    ## 
    ## All 6 chains finished successfully.
    ## Mean chain execution time: 14.0 seconds.
    ## Total execution time: 16.5 seconds.
=======
    ## sccomp says: the composition design matrix has columns: (Intercept), typecancer

    ## sccomp says: the variability design matrix has columns: (Intercept), typecancer
>>>>>>> master

``` r
res
```

<<<<<<< HEAD
    ## # A tibble: 72 × 13
    ##    cell_group parameter   c_lower c_effect c_upper    c_pH0    c_FDR v_lower
    ##    <chr>      <chr>         <dbl>    <dbl>   <dbl>    <dbl>    <dbl>   <dbl>
    ##  1 B1         (Intercept)  0.587     0.752  0.916  0        0          -5.11
    ##  2 B1         typecancer  -1.17     -0.850 -0.501  0.000751 0.000250    1.54
    ##  3 B2         (Intercept)  0.173     0.417  0.662  0.0388   0.00171    -4.45
    ##  4 B2         typecancer  -1.04     -0.559 -0.0567 0.0746   0.0205      1.48
    ##  5 B3         (Intercept) -0.561    -0.363 -0.152  0.0593   0.00548    -5.62
    ##  6 B3         typecancer  -0.545    -0.128  0.301  0.636    0.184       1.76
    ##  7 BM         (Intercept) -1.28     -1.05  -0.833  0        0          -6.10
    ##  8 BM         typecancer  -0.734    -0.301  0.145  0.324    0.0821      1.24
    ##  9 CD4 1      (Intercept)  0.399     0.534  0.673  0        0          -5.50
    ## 10 CD4 1      typecancer  -0.0465    0.230  0.515  0.417    0.141       1.76
    ## # … with 62 more rows, and 5 more variables: v_effect <dbl>, v_upper <dbl>,
    ## #   v_pH0 <dbl>, v_FDR <dbl>, count_data <list>
=======
    ## # A tibble: 72 × 14
    ##    cell_group parameter   covariate c_lower c_effect c_upper    c_pH0     c_FDR
    ##    <chr>      <chr>       <chr>       <dbl>    <dbl>   <dbl>    <dbl>     <dbl>
    ##  1 B1         (Intercept) <NA>       1.03      1.24   1.45   0        0        
    ##  2 B1         typecancer  type      -1.26     -0.915 -0.564  0        0        
    ##  3 B2         (Intercept) <NA>       0.517     0.766  1.00   0        0        
    ##  4 B2         typecancer  type      -1.06     -0.635 -0.194  0.0265   0.00428  
    ##  5 B3         (Intercept) <NA>      -0.554    -0.317 -0.0577 0.167    0.0104   
    ##  6 B3         typecancer  type      -0.647    -0.196  0.242  0.512    0.129    
    ##  7 BM         (Intercept) <NA>      -1.23     -0.918 -0.591  0.000501 0.0000501
    ##  8 BM         typecancer  type      -0.723    -0.282  0.176  0.350    0.0761   
    ##  9 CD4 1      (Intercept) <NA>       0.226     0.407  0.592  0.0158   0.00109  
    ## 10 CD4 1      typecancer  type      -0.0857    0.187  0.473  0.538    0.158    
    ## # … with 62 more rows, and 6 more variables: v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, count_data <list>
>>>>>>> master

Plot 1D significance plot

``` r
plots = plot_summary(res)
```

<<<<<<< HEAD
    ## Running standalone generated quantities after 6 MCMC chains, 1 chain at a time ...
    ## 
    ## Chain 1 finished in 0.0 seconds.
    ## Chain 2 finished in 0.0 seconds.
    ## Chain 3 finished in 0.0 seconds.
    ## Chain 4 finished in 0.0 seconds.
    ## Chain 5 finished in 0.0 seconds.
    ## Chain 6 finished in 0.0 seconds.
    ## 
    ## All 6 chains finished successfully.
    ## Mean chain execution time: 0.0 seconds.
    ## Total execution time: 9.7 seconds.

=======
>>>>>>> master
    ## Joining, by = c("sample", "cell_group")
    ## Joining, by = c("cell_group", "type")

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-16-1.png)<!-- -->

Plot 2D significance plot. Data points are cell groups. Error bars are
the 95% credible interval. The dashed lines represent the default
threshold fold change for which the probabilities (c_pH0, v_pH0) are
calculated. pH0 of 0 represent the rejection of the null hypothesis,
that no effect is observed.

This plot is provided only if differential variability has been tested.
The differential variability estimates are reliable only if the linear
association between mean and variability for `(intercept)` (left-hand
side facet) is satisfied. A scatterplot (beside the Intercept) is
provided for each of the categories of interest. The for each category
of interest, the composition and variability effects should be generally
uncorrelated.

``` r
plots$credible_intervals_2D
```

![](inst/figures/unnamed-chunk-17-1.png)<!-- -->
