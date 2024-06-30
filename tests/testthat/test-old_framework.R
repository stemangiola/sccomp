library(dplyr)
library(sccomp)
data("seurat_obj")
data("sce_obj")
data("counts_obj")

result = 
  seurat_obj |>
  sccomp_glm(
    formula_composition = ~ 0 + type + (type | group__),
    formula_variability = ~ 0 + type,
    sample, cell_group,
    check_outliers = TRUE,
    approximate_posterior_inference = FALSE,
    contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
    cores = 1,
    mcmc_seed = 42,       
    max_sampling_iterations = 1000
  ) |> 
  suppressWarnings()

test_that("Generate data",{
  
  
  result |> 
    
    
    sccomp_replicate() |>
    nrow() |>
    expect_equal(600) 
  
  
})


test_that("remove unwanted variation",{
  
  result |>
    sccomp_remove_unwanted_variation(~ type) 
  
})

# test_that("plot test composition",{
#   
#   plot_summary(result) 
#   
#   
# })

test_that("plot test composition",{
  
  result |> 
    test_contrasts() |> 
    lifecycle::expect_deprecated()  
  
  
})
