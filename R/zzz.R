.onAttach <- function(lib, pkg) {
	packageStartupMessage(
		"NOTE: hdbayes depends on the R package cmdstanr, which is not available on CRAN.
      The cmdstanr package may be installed at https://mc-stan.org/cmdstanr/."
	)
}