sccomp - Outlier-aware and count-based compositional analysis of
single-cell data
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidyseurat/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/tidyseurat/actions/)

<!-- badges: end -->

<a href="https://www.youtube.com/watch?v=R_lt58We9nA&ab_channel=RConsortium" target="_blank">
<img src="https://img.youtube.com/vi/R_lt58We9nA/mqdefault.jpg" alt="Watch the video" width="280" height="180" border="10" />
</a>

# <img src="inst/logo-01.png" height="139px" width="120px"/>

`sccomp` tests differences in cell type proportions from single-cell
data. It is robust against outliers, it models continuous and discrete
factors, and capable of random-effect/intercept modelling.

Please cite [PNAS - sccomp: Robust differential composition and
variability analysis for single-cell
data](https://www.pnas.org/doi/full/10.1073/pnas.2203828120)

## Characteristics

- Complex linear models with continuous and categorical covariates
- Multilevel modelling, with population fixed and random
  effects/intercept
- Modelling data from counts
- Testing differences in cell-type proportionality
- Testing differences in cell-type specific variability
- Cell-type information share for variability adaptive shrinkage
- Testing differential variability
- Probabilistic outlier identification
- Cross-dataset learning (hyperpriors).

# Installation

**Bioconductor**

``` r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("sccomp")
```

**Github**

``` r
devtools::install_github("stemangiola/sccomp")
```

| Function                           | Description                                                                                                                 |
|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| `sccomp_estimate`                  | Fit the model onto the data, and estimate the coefficients                                                                  |
| `sccomp_remove_outliers`           | Identify outliers probabilistically based on the model fit, and exclude them from the estimation                            |
| `sccomp_test`                      | Calculate the probability that the coefficients are outside the H0 interval (i.e. test_composition_above_logit_fold_change) |
| `sccomp_replicate`                 | Simulate data from the model, or part of the model                                                                          |
| `sccomp_predict`                   | Predicts proportions, based on the mode, or part of the model                                                               |
| `sccomp_remove_unwanted_variation` | Removes the variability for unwanted factors                                                                                |
| `plot`                             | Plors summary plots to asses significance                                                                                   |

# Analysis

`sccomp` can model changes in composition and variability. By default,
the formula for variability is either `~1`, which assumes that the
cell-group variability is independent of any covariate or
`~ factor_of_interest`, which assumes that the model is dependent on the
factor of interest only. The variability model must be a subset of the
model for composition.

## Binary factor

### From proportions

If counts are available we strongly discourage the use of proportions,
as an important source of uncertainty (i.e. for rare groups/cell-types)
is not modeled.

Proportions should be bigger than 0. Assuming that 0s derive from a
precision threshold (e.g. deconvolution), 0s are converted to the
smaller non 0 value.

``` r
res = 
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    .sample = sample,
    .cell_group = cell_group,
    .count = proportion, 
    bimodal_mean_variability_association = TRUE,
    cores = 1,
    verbose = FALSE
  ) |> 
  sccomp_remove_outliers(verbose = FALSE) |> 
    sccomp_test(test_composition_above_logit_fold_change = 0.2)
```

    ## # A tibble: 72 × 18
    ##    cell_group parameter   factor c_lower c_effect c_upper  c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>       <chr>    <dbl>    <dbl>   <dbl>  <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercept) <NA>    0.811     1.07   1.32   0      0           NaN
    ##  2 B1         typecancer  type   -0.810    -0.466 -0.117  0.0650 0.0105      NaN
    ##  3 B2         (Intercept) <NA>    0.228     0.535  0.828  0.0170 0.00127     NaN
    ##  4 B2         typecancer  type   -0.907    -0.475 -0.0285 0.114  0.0293      NaN
    ##  5 B3         (Intercept) <NA>   -0.685    -0.419 -0.153  0.0500 0.00339     NaN
    ##  6 B3         typecancer  type   -0.676    -0.314  0.0565 0.29   0.0887      NaN
    ##  7 BM         (Intercept) <NA>   -1.39     -1.13  -0.864  0      0           NaN
    ##  8 BM         typecancer  type   -0.548    -0.203  0.160  0.49   0.138       NaN
    ##  9 CD4 1      (Intercept) <NA>    0.126     0.310  0.509  0.128  0.0206      NaN
    ## 10 CD4 1      typecancer  type   -0.0734    0.176  0.453  0.571  0.169       NaN

    ## # ℹ 62 more rows
    ## # ℹ 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>, v_R_k_hat <dbl>,
    ## #   count_data <list>

# Visualisation

## Summary plots

``` r
plots = plot(res) 
```

    ## Joining with `by = join_by(cell_group, sample)`
    ## Joining with `by = join_by(cell_group, type)`

### Visualise the descriptive adequacy of the model

A plot of group proportion, faceted by groups. The blue boxplots
represent the posterior predictive check. If the model is likely to be
descriptively adequate to the data, the blue box plot should roughly
overlay with the black box plot, which represents the observed data. The
outliers are coloured in red. A box plot will be returned for every
(discrete) covariate present in `formula_composition`. The colour coding
represents the significant associations for composition and/or
variability.

``` r
plots$boxplot
```

    ## [[1]]

![](inst/figures/unnamed-chunk-10-1.png)<!-- -->

### Mean/variability association

Plot 2D significance plot. Data points are cell groups. Error bars are
the 95% credible interval. The dashed lines represent the default
threshold fold change for which the probabilities (c_pH0, v_pH0) are
calculated. pH0 of 0 represent the rejection of the null hypothesis that
no effect is observed.

This plot is provided only if differential variability has been tested.
The differential variability estimates are reliable only if the linear
association between mean and variability for `(intercept)` (left-hand
side facet) is satisfied. A scatterplot (besides the Intercept) is
provided for each category of interest. The for each category of
interest, the composition and variability effects should be generally
uncorrelated.

``` r
plots$credible_intervals_2D
```

![](inst/figures/unnamed-chunk-11-1.png)<!-- -->

### Effects significance

A plot of estimates of differential composition (c\_) on the x-axis and
differential variability (v\_) on the y-axis. The error bars represent
95% credible intervals. The dashed lines represent the minimal effect
that the hypothesis test is based on. An effect is labelled as
significant if bigger than the minimal effect according to the 95%
credible interval. Facets represent the covariates in the model.

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-12-1.png)<!-- -->

## Visualisation of the MCMC chains from the posterior distribution

It is possible to directly evaluate the posterior distribution. In this
example, we plot the Monte Carlo chain for the slope parameter of the
first cell type. We can see that it has converged and is negative with
probability 1.

``` r
res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")
```

![](inst/figures/unnamed-chunk-13-1.png)<!-- -->

Of the output table, the estimate columns start with the prefix `c_`
indicate `composition`, or with `v_` indicate `variability` (when
formula_variability is set).

## Contrasts

``` r
counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ 0 + type, 
    .sample = sample,
    .cell_group = cell_group, 
    .count = proportion,
    bimodal_mean_variability_association = TRUE,
    cores = 1 ,
    verbose = FALSE
  ) |> 
  sccomp_remove_outliers(verbose = FALSE) |> 
    sccomp_test(
      contrasts =  c("typecancer - typehealthy", "typehealthy - typecancer"),
      test_composition_above_logit_fold_change = 0.2
    )
```

    ## # A tibble: 60 × 18
    ##    cell_group  parameter factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>       <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature  typecanc… <NA>    -2.07    -1.66   -1.24  0       0            NA
    ##  2 B immature  typeheal… <NA>     1.24     1.66    2.07  0       0            NA
    ##  3 B mem       typecanc… <NA>    -2.39    -1.78   -1.06  0       0            NA
    ##  4 B mem       typeheal… <NA>     1.06     1.78    2.39  0       0            NA
    ##  5 CD4 cm S10… typecanc… <NA>    -1.25    -0.948  -0.642 0       0            NA
    ##  6 CD4 cm S10… typeheal… <NA>     0.642    0.948   1.25  0       0            NA
    ##  7 CD4 cm hig… typecanc… <NA>     0.514    1.62    2.64  0.00400 1.27e-3      NA
    ##  8 CD4 cm hig… typeheal… <NA>    -2.64    -1.62   -0.514 0.00400 1.27e-3      NA
    ##  9 CD4 cm rib… typecanc… <NA>     0.553    1.21    1.85  0.00200 6.25e-4      NA
    ## 10 CD4 cm rib… typeheal… <NA>    -1.85    -1.21   -0.553 0.00200 6.25e-4      NA
    ## # ℹ 50 more rows
    ## # ℹ 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>, v_R_k_hat <dbl>,
    ## #   count_data <list>


## Categorical factor (e.g. Bayesian ANOVA)

This is achieved through model comparison with `loo`. In the following
example, the model with association with factors better fits the data
compared to the baseline model with no factor association. For
comparisons `check_outliers` must be set to FALSE as the leave-one-out
must work with the same amount of data, while outlier elimination does
not guarantee it.

If `elpd_diff` is away from zero of \> 5 `se_diff` difference of 5, we
are confident that a model is better than the other
[reference](https://discourse.mc-stan.org/t/interpreting-elpd-diff-loo-package/1628/2?u=stemangiola).
In this case, -79.9 / 11.5 = -6.9, therefore we can conclude that model
one, the one with factor association, is better than model two.

``` r
library(loo)

# Fit first model
model_with_factor_association = 
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    .sample =  sample, 
    .cell_group = cell_group, 
    .count = proportion,
    bimodal_mean_variability_association = TRUE,
    cores = 1, 
    enable_loo = TRUE, # Needed for model comparison and ANOVA
    verbose = FALSE
  )


# Fit second model
model_without_association = 
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ 1, 
    .sample =  sample, 
    .cell_group = cell_group, 
    .count = proportion,
    bimodal_mean_variability_association = TRUE,
    cores = 1 , 
    enable_loo = TRUE, # Needed for model comparison and ANOVA
    verbose = FALSE
  )


# Compare models
loo_compare(
  model_with_factor_association |> attr("fit") |> loo(),
  model_without_association |> attr("fit") |> loo()
)
```

    ##        elpd_diff se_diff
    ## model1   0.0       0.0  
    ## model2 -83.2      15.2

## Differential variability, binary factor

We can model the cell-group variability also dependent on the type, and
so test differences in variability

``` r
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    formula_variability = ~ type,
    .sample = sample,
    .cell_group = cell_group,
    .count = proportion,
    bimodal_mean_variability_association = TRUE,
    cores = 1 ,
    verbose = FALSE
  ) |> 
  sccomp_remove_outliers(verbose = FALSE) |> 
    sccomp_test(
      test_composition_above_logit_fold_change = 0.2
    )
```

``` r
res
```

    ## # A tibble: 60 × 18
    ##    cell_group  parameter factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>       <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature  (Interce… <NA>     0.605    0.942  1.30   0       0           NaN
    ##  2 B immature  typeheal… type     1.09     1.58   2.06   0       0           NaN
    ##  3 B mem       (Interce… <NA>    -1.36    -1.03  -0.697  0       0           NaN
    ##  4 B mem       typeheal… type     1.68     2.15   2.68   0       0           NaN
    ##  5 CD4 cm S10… (Interce… <NA>     1.80     2.06   2.32   0       0           NaN
    ##  6 CD4 cm S10… typeheal… type     0.501    0.823  1.14   0.00100 1.00e-4     NaN
    ##  7 CD4 cm hig… (Interce… <NA>    -1.05    -0.509  0.0229 0.123   1.39e-2     NaN
    ##  8 CD4 cm hig… typeheal… type    -2.37    -1.61  -0.849  0       0           NaN
    ##  9 CD4 cm rib… (Interce… <NA>     0.198    0.566  0.950  0.0260  3.05e-3     NaN
    ## 10 CD4 cm rib… typeheal… type    -2.35    -1.87  -1.36   0       0           NaN
    ## # ℹ 50 more rows

    ## # ℹ 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>, v_R_k_hat <dbl>,
    ## #   count_data <list>

# Suggested settings

## For single-cell RNA sequencing

We recommend setting `bimodal_mean_variability_association  = TRUE`. The
bimodality of the mean-variability association can be confirmed from the
plots\$credible_intervals_2D (see below).

## For CyTOF and microbiome data

We recommend setting `bimodal_mean_variability_association  = FALSE`
(Default).


# Visualisation

## Summary plots

``` r
plots = plot(res) 
```

    ## Joining with `by = join_by(cell_group, sample)`
    ## Joining with `by = join_by(cell_group, type)`

A plot of group proportion, faceted by groups. The blue boxplots
represent the posterior predictive check. If the model is likely to be
descriptively adequate to the data, the blue box plot should roughly
overlay with the black box plot, which represents the observed data. The
outliers are coloured in red. A box plot will be returned for every
(discrete) covariate present in `formula_composition`. The colour coding
represents the significant associations for composition and/or
variability.

``` r
plots$boxplot
```

    ## [[1]]

![](inst/figures/unnamed-chunk-15-1.png)<!-- -->

A plot of estimates of differential composition (c\_) on the x-axis and
differential variability (v\_) on the y-axis. The error bars represent
95% credible intervals. The dashed lines represent the minimal effect
that the hypothesis test is based on. An effect is labelled as
significant if bigger than the minimal effect according to the 95%
credible interval. Facets represent the covariates in the model.

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-16-1.png)<!-- -->

## Visualisation of the MCMC chains from the posterior distribution

It is possible to directly evaluate the posterior distribution. In this
example, we plot the Monte Carlo chain for the slope parameter of the
first cell type. We can see that it has converged and is negative with
probability 1.

``` r
res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")
```

![](inst/figures/unnamed-chunk-17-1.png)<!-- -->

Plot 1D significance plot

``` r
plots = plot(res)
```

    ## Joining with `by = join_by(cell_group, sample)`
    ## Joining with `by = join_by(cell_group, type)`

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-18-1.png)<!-- -->

Plot 2D significance plot. Data points are cell groups. Error bars are
the 95% credible interval. The dashed lines represent the default
threshold fold change for which the probabilities (c_pH0, v_pH0) are
calculated. pH0 of 0 represent the rejection of the null hypothesis that
no effect is observed.

This plot is provided only if differential variability has been tested.
The differential variability estimates are reliable only if the linear
association between mean and variability for `(intercept)` (left-hand
side facet) is satisfied. A scatterplot (besides the Intercept) is
provided for each category of interest. The for each category of
interest, the composition and variability effects should be generally
uncorrelated.

``` r
plots$credible_intervals_2D
```

![](inst/figures/unnamed-chunk-19-1.png)<!-- -->

# Multilevel modelling

`sccomp` is cabable of estimating population (i.e. fixed) and group
(i.e. random) effects. The formulation is analogous to the `lme4`
package and `brms`.

!! For now, only one grouping is allowed (e.g. group2\_\_).

``` r
res = 
  seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ type + continuous_covariate + (type | group2__), 
    formula_variability = ~ type,
    .sample = sample,
    .cell_group = cell_group,
    bimodal_mean_variability_association = TRUE,
    cores = 1,
    verbose = FALSE,
    variational_inference = FALSE # For this more complex model use full HMC inference
  ) |> 
  sccomp_remove_outliers(variational_inference = FALSE, verbose = FALSE) |> 
    sccomp_test(
      test_composition_above_logit_fold_change = 0.2
    )

res
```

    ## # A tibble: 210 × 18
    ##    cell_group parameter  factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>      <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature (Intercep… <NA>     0.461  1.00      1.58  0.00425 1.08e-3   1336.
    ##  2 B immature typehealt… type     0.614  1.31      1.94  0.00100 5.01e-4   1527.
    ##  3 B immature continuou… conti…  -0.187  0.0599    0.340 0.847   6.12e-1   4872.
    ##  4 B immature (Intercep… <NA>    -0.606 -0.0682    0.399 0.723   6.31e-1     NA 
    ##  5 B immature typehealt… <NA>    -0.552 -0.00762   0.570 0.771   6.59e-1     NA 
    ##  6 B immature (Intercep… <NA>    -0.491  0.00427   0.466 0.803   7.05e-1     NA 
    ##  7 B immature typehealt… <NA>    -0.358  0.146     0.742 0.590   5.32e-1     NA 
    ##  8 B mem      (Intercep… <NA>    -1.07  -0.427     0.395 0.270   5.61e-2   2104.
    ##  9 B mem      typehealt… type     0.656  1.52      2.30  0.00325 1.42e-3   2344.
    ## 10 B mem      continuou… conti…  -0.225  0.0734    0.388 0.792   5.93e-1   4059.
    ## # ℹ 200 more rows
    ## # ℹ 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>, v_R_k_hat <dbl>,
    ## #   count_data <list>


# Removal of unwanted variation

After you model your dataset, you can remove the unwanted variation from
your input data, **for visualisation purposes**

We decide to just keep the type population (i.e. fixed) effect for
abundance, and do not keep it for variability.

``` r
res |> sccomp_remove_unwanted_variation(~type)
```

    ## sccomp says: calculating residuals

    ## sccomp says: regressing out unwanted factors

    ## # A tibble: 600 × 5
    ##    sample       cell_group adjusted_proportion adjusted_counts logit_residuals
    ##    <chr>        <chr>                    <dbl>           <dbl>           <dbl>
    ##  1 10x_6K       B immature              0.0674           316.          -0.545 
    ##  2 10x_8K       B immature              0.170           1278.           0.548 
    ##  3 GSE115189    B immature              0.134            315.           0.229 
    ##  4 SCP345_580   B immature              0.0827           476.          -0.299 
    ##  5 SCP345_860   B immature              0.141            905.           0.288 
    ##  6 SCP424_pbmc1 B immature              0.102            273.          -0.0679
    ##  7 SCP424_pbmc2 B immature              0.182            544.           0.635 
    ##  8 SCP591       B immature              0.0311            17.7         -1.35  
    ##  9 SI-GA-E5     B immature              0.0278           116.          -0.620 
    ## 10 SI-GA-E7     B immature              0.0989           726.           0.729 
    ## # ℹ 590 more rows

## The old framework

The new tidy framework was introduced in 2024, two, understand the
differences and improvements. Compared to the old framework, please read
this [blog
