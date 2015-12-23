wasserman_normal_sim <- function(m, pi0, xi_min, xi_max, seed=NULL){
  if (!is.null(seed)) set.seed(seed)

  X   <- runif(m, min=xi_min, max=xi_max)
  H   <- rbinom(m,1,1-pi0)
  Z   <- rnorm(m, H*X)
  pvalue <- 1-pnorm(Z)
  simDf <- data.frame(pvalue=pvalue, filterstat=X,H=H, Z=Z)
}


sim <- wasserman_normal_sim(10000,0.85, 0, 3, seed=1)
ihw_res1 <- ihw(sim$pvalue, sim$filterstat, .1, nbins=10)

expect_true(all(apply(weights(ihw_res1, levels_only=T),2, IHW:::total_variation) <= ihw_res1@regularization_term + 10^(-12)))

# try same simulation with other value of alpha
ihw_res1_lower_alpha <- ihw(sim$pvalue, sim$filterstat, .01, nbins=10)
testthat::expect_less_than( rejections(ihw_res1_lower_alpha), rejections(ihw_res1))


sim$group <- as.factor(IHW:::groups_by_filter(sim$filterstat, 10))
ihw_res2 <- ihw(sim$pvalue, sim$group, .1)

expect_equal(rejections(ihw_res1), rejections(ihw_res2))

# quick test for show method
expect_equal(capture.output(ihw_res1), capture.output(ihw_res2))

plot_ihw(ihw_res2)


ihw_res3 <- ihw(sim$pvalue, sim$group, .1, covariate_type="nominal")

expect_true(all(apply(weights(ihw_res3, levels_only=T),2, IHW:::uniform_deviation) <= ihw_res3@regularization_term + 10^(-12)))

# now test small inputs.
sim <- wasserman_normal_sim(200,0.85, 0, 3, seed=1)
ihw_res_small <- ihw(sim$pvalue, sim$filterstat, .1)
plot_ihw(ihw_res_small, scale="nominal")

ihw_res_small2 <- ihw(sim$pvalue, sim$filterstat, .05, nbins=2)