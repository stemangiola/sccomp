library(readr)

glm_multi_beta_binomial = read_file("inst/stan/glm_multi_beta_binomial.stan") 
save(glm_multi_beta_binomial, file="data/glm_multi_beta_binomial.rda", compress = "xz")

glm_multi_beta_binomial_generate = read_file("inst/stan/glm_multi_beta_binomial_generate_date.stan") 
save(glm_multi_beta_binomial_generate, file="data/glm_multi_beta_binomial_generate.rda", compress = "xz")

glm_multi_beta_binomial_simulate = read_file("inst/stan/glm_multi_beta_binomial_simulate_data.stan") 
save(glm_multi_beta_binomial_simulate, file="data/glm_multi_beta_binomial_simulate.rda", compress = "xz")

devtools::document()
devtools::install()
