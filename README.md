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

# Characteristics

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

Of the output table, the estimate columns start with the prefix `c_`
indicate `composition`, or with `v_` indicate `variability` (when
formula_variability is set).

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
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_remove_outliers(cores = 1, verbose = FALSE) |> # Optional
  sccomp_test()

res
```

    ## # A tibble: 72 × 18
    ##    cell_group parameter  factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>      <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercep… <NA>    0.905    1.25    1.63   0       0           NaN
    ##  2 B1         typecancer type   -1.27    -0.783  -0.332  0.00600 3.  e-3     NaN
    ##  3 B2         (Intercep… <NA>    0.397    0.815   1.19   0.00100 5.26e-5     NaN
    ##  4 B2         typecancer type   -1.40    -0.818  -0.222  0.0200  7.63e-3     NaN
    ##  5 B3         (Intercep… <NA>   -0.661   -0.311   0.0608 0.288   2.91e-2     NaN
    ##  6 B3         typecancer type   -0.888   -0.370   0.130  0.257   7.27e-2     NaN
    ##  7 BM         (Intercep… <NA>   -1.33    -0.942  -0.569  0       0           NaN
    ##  8 BM         typecancer type   -0.854   -0.330   0.205  0.31    1.01e-1     NaN
    ##  9 CD4 1      (Intercep… <NA>    0.0815   0.414   0.766  0.091   8.74e-3     NaN
    ## 10 CD4 1      typecancer type   -0.347    0.0738  0.504  0.713   3.31e-1     NaN
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

We can plot the relationship between abundance and variability. As we
can see below, they are positively correlated, you also appreciate that
this relationship is by model for single cell RNA sequencing data.

`sccomp` models, these relationship to obtain a shrinkage effect on the
estimates of both the abundance and the variability. This shrinkage is
adaptive as it is modelled jointly, thanks for Bayesian inference.

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

## Contrasts

``` r
counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ 0 + type, 
    .sample = sample,
    .cell_group = cell_group, 
    .count = proportion,
    bimodal_mean_variability_association = TRUE,
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_test( contrasts =  c("typecancer - typehealthy", "typehealthy - typecancer"))
```

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
    enable_loo = TRUE
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
    enable_loo = TRUE
  )

# Compare models
loo_compare(
  model_with_factor_association |> attr("fit") |> loo(),
  model_without_association |> attr("fit") |> loo()
)
```

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
    cores = 1, verbose = FALSE
  )
```

    ## # A tibble: 72 × 14
    ##    cell_group parameter   factor c_lower c_effect c_upper c_n_eff c_R_k_hat
    ##    <chr>      <chr>       <chr>    <dbl>    <dbl>   <dbl>   <dbl>     <dbl>
    ##  1 B1         (Intercept) <NA>    1.01      1.27   1.52       NaN      3.83
    ##  2 B1         typecancer  type   -0.927    -0.600 -0.267      NaN      3.86
    ##  3 B2         (Intercept) <NA>    0.359     0.715  1.05       NaN      3.86
    ##  4 B2         typecancer  type   -1.12     -0.668 -0.168      NaN      3.84
    ##  5 B3         (Intercept) <NA>   -0.673    -0.370 -0.0775     NaN      3.85
    ##  6 B3         typecancer  type   -0.581    -0.161  0.271      NaN      3.87
    ##  7 BM         (Intercept) <NA>   -1.21     -0.912 -0.599      NaN      3.86
    ##  8 BM         typecancer  type   -0.933    -0.502 -0.0694     NaN      3.85
    ##  9 CD4 1      (Intercept) <NA>    0.193     0.382  0.549      NaN      3.88
    ## 10 CD4 1      typecancer  type    0.0359    0.260  0.515      NaN      3.86
    ## # ℹ 62 more rows
    ## # ℹ 6 more variables: v_lower <dbl>, v_effect <dbl>, v_upper <dbl>,
    ## #   v_n_eff <dbl>, v_R_k_hat <dbl>, count_data <list>

``` r
res
```

    ## # A tibble: 72 × 18
    ##    cell_group parameter  factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>      <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercep… <NA>    0.905    1.25    1.63   0       0           NaN
    ##  2 B1         typecancer type   -1.27    -0.783  -0.332  0.00600 3.  e-3     NaN
    ##  3 B2         (Intercep… <NA>    0.397    0.815   1.19   0.00100 5.26e-5     NaN
    ##  4 B2         typecancer type   -1.40    -0.818  -0.222  0.0200  7.63e-3     NaN
    ##  5 B3         (Intercep… <NA>   -0.661   -0.311   0.0608 0.288   2.91e-2     NaN
    ##  6 B3         typecancer type   -0.888   -0.370   0.130  0.257   7.27e-2     NaN
    ##  7 BM         (Intercep… <NA>   -1.33    -0.942  -0.569  0       0           NaN
    ##  8 BM         typecancer type   -0.854   -0.330   0.205  0.31    1.01e-1     NaN
    ##  9 CD4 1      (Intercep… <NA>    0.0815   0.414   0.766  0.091   8.74e-3     NaN
    ## 10 CD4 1      typecancer type   -0.347    0.0738  0.504  0.713   3.31e-1     NaN
    ## # ℹ 62 more rows
    ## # ℹ 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>, v_R_k_hat <dbl>,
    ## #   count_data <list>

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

# Suggested settings

## For single-cell RNA sequencing

We recommend setting `bimodal_mean_variability_association  = TRUE`. The
bimodality of the mean-variability association can be confirmed from the
plots\$credible_intervals_2D (see below).

## For CyTOF and microbiome data

We recommend setting `bimodal_mean_variability_association  = FALSE`
(Default).

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

# Removal of unwanted variation

After you model your dataset, you can remove the unwanted variation from
your input data, **for visualisation purposes**

We decide to just keep the type population (i.e. fixed) effect for
abundance, and do not keep it for variability.

``` r
res |> sccomp_remove_unwanted_variation(~type)
```

## Visualisation of the MCMC chains from the posterior distribution

It is possible to directly evaluate the posterior distribution. In this
example, we plot the Monte Carlo chain for the slope parameter of the
first cell type. We can see that it has converged and is negative with
probability 1.

``` r
res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")
```

![](inst/figures/unnamed-chunk-19-1.png)<!-- -->

Plot 1D significance plot

``` r
plots = res |> sccomp_test() |> plot()
```

    ## Joining with `by = join_by(cell_group, sample)`
    ## Joining with `by = join_by(cell_group, type)`

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-20-1.png)<!-- -->

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

![](inst/figures/unnamed-chunk-21-1.png)<!-- -->

## The old framework

The new tidy framework was introduced in 2024, two, understand the
differences and improvements. Compared to the old framework, please read
this [blog
post](https://tidyomics.github.io/tidyomicsBlog/post/2023-12-07-tidy-sccomp/).
