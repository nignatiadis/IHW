wasserman_normal_sim <- function(m, pi0, xi_min, xi_max, seed=NULL){
  if (!is.null(seed)) set.seed(seed)

  X   <- runif(m, min=xi_min, max=xi_max)
  H   <- rbinom(m,1,1-pi0)
  Z   <- rnorm(m, H*X)
  pvalue <- 1-pnorm(Z)
  simDf <- data.frame(pvalue=pvalue, filterstat=X,H=H, Z=Z)
}


sim <- wasserman_normal_sim(10000,0.85, 0, 3, seed=1)
sim$group <- as.factor(IHW:::groups_by_filter(sim$filterstat, 10))

ihw_res1 <- ihw(sim$pvalue, sim$filterstat, .1, nbins=10)

expect_true(all(apply(weights(ihw_res1, levels_only=T),2, IHW:::total_variation) <= ihw_res1@regularization_term + 10^(-12)))


# now test the formula interface
ihw_res1_formula1 <- ihw(pvalue~filterstat, data=sim, alpha=.1, nbins=10)
ihw_res1_formula2 <- ihw(sim$pvalue~sim$filterstat, alpha=.1, nbins=10)

expect_equal(rejections(ihw_res1), rejections(ihw_res1_formula1))
expect_equal(rejections(ihw_res1), rejections(ihw_res1_formula2))


# try same simulation with other value of alpha
ihw_res1_lower_alpha <- ihw(sim$pvalue, sim$filterstat, .01, nbins=10)
testthat::expect_less_than( rejections(ihw_res1_lower_alpha), rejections(ihw_res1))

# try with only 1 fold
expect_message(ihw_res1_single_fold <- ihw(sim$pvalue, sim$filterstat, .1, nbins=10, nfolds=1))

ihw_res2 <- ihw(sim$pvalue, sim$group, .1)

expect_equal(rejections(ihw_res1), rejections(ihw_res2))


plot(ihw_res2)


ihw_res3 <- ihw(sim$pvalue, sim$group, .1, covariate_type="nominal")

expect_true(all(apply(weights(ihw_res3, levels_only=T),2, IHW:::uniform_deviation) <= ihw_res3@regularization_term + 10^(-12)))

# now test small inputs.
sim_small <- wasserman_normal_sim(200,0.85, 0, 3, seed=1)
ihw_res_small <- ihw(sim_small$pvalue, sim_small$filterstat, .1)
plot(ihw_res_small, scale="nominal")

ihw_res_small2 <- ihw(sim_small$pvalue, sim_small$filterstat, .05, nbins=2)


# now test ihwResult class getters or methods

# nbins
expect_equal(nbins(ihw_res1), 10L)
expect_equal(nbins(ihw_res1), nbins(ihw_res1_single_fold))

# nfolds
expect_equal(nfolds(ihw_res1), 5L) # the default choice
expect_equal(nfolds(ihw_res1_single_fold), 1L) # the default choice

n <- nrow(ihw_res1)
expect_equal(n, 10000)
expect_equal(nrow(ihw_res1_single_fold), n)

# methods which return vector equal to number of hypotheses

mymethods <- c(adj_pvalues, weights, thresholds, pvalues, weighted_pvalues, covariates,
        groups_factor ,rejected_hypotheses)

lengths <- sapply(mymethods, function(f) length(f(ihw_res1)))
expect_true(all(lengths == n))

# check if weight budget fullfilled
expect_equal( sum(weights(ihw_res1)) , nrow(ihw_res1))
expect_equal( sum(weights(ihw_res1_single_fold)) , nrow(ihw_res1_single_fold))

# check if levels_only works
ws_sorted1 <- sort(unique(as.numeric(weights(ihw_res1, levels_only=TRUE))))
ws_sorted2 <- sort(unique(weights(ihw_res1, levels_only=FALSE)))
expect_equal(ws_sorted1, ws_sorted2)

expect_equal(dim(weights(ihw_res1, levels_only=TRUE)), c(nbins(ihw_res1), nfolds(ihw_res1)))
expect_equal(dim(weights(ihw_res1_single_fold, levels_only=TRUE)),
         c(nbins(ihw_res1_single_fold), nfolds(ihw_res1_single_fold)))

# same for thresholds
ts_sorted1 <- sort(unique(as.numeric(thresholds(ihw_res1, levels_only=TRUE))))
ts_sorted2 <- sort(unique(thresholds(ihw_res1, levels_only=FALSE)))
expect_equal(ts_sorted1, ts_sorted2)

expect_equal(dim(thresholds(ihw_res1, levels_only=TRUE)), c(nbins(ihw_res1), nfolds(ihw_res1)))
expect_equal(dim(thresholds(ihw_res1_single_fold, levels_only=TRUE)),
         c(nbins(ihw_res1_single_fold), nfolds(ihw_res1_single_fold)))

expect_equal(covariate_type(ihw_res1), "ordinal")

expect_equal(alpha(ihw_res1), 0.1)
expect_equal(alpha(ihw_res1_lower_alpha), 0.01)

expect_true(is.data.frame(as.data.frame(ihw_res1)))

# quick test for show method
expect_equal(capture.output(ihw_res1), capture.output(ihw_res2))

# now let's also test if ECDF method runs
sim3 <- wasserman_normal_sim(2000,0.85, 0, 3, seed=1)
ihw_naive <- ihw(sim3$pvalue, sim3$filterstat, .1, nfolds=1L, nbins=3L, lambdas=Inf, distrib_estimator="ECDF")
# should have increased rejections compared to BH
# also opportunity to test get_bh_threshold
expect_less_than( sum(sim3$pvalue <= get_bh_threshold(sim3$pvalue, .1)), rejections(ihw_naive))




#--------- Filtered method ------------------------------------------#
#-- Test if IHW with input of only p-values <= threshold , still works
#--------------------------------------------------------------------#

mgroups <- table(sim$group)
sim_filt <- subset(sim, sim$pvalue <= 0.5)
ihw_res1_filtered_single_fold <- ihw(sim_filt$pvalue, sim_filt$group,.1, 
                                     nfolds=1, m_groups=mgroups)

expect_equal(rejections(ihw_res1_single_fold), rejections(ihw_res1_filtered_single_fold))

t1 <- thresholds(ihw_res1_filtered_single_fold, levels_only=T)
t2 <- thresholds(ihw_res1_single_fold, levels_only=T)
expect_equal(t1,t2)

#-------------------------------------------------------------------
#--------- Check if manual definition of folds works ---------------
#-------------------------------------------------------------------
ihw_res_lambda <- ihw(sim$pvalue, sim$filterstat, .1, nbins=10, lambda=10)
folds <- ihw_res_lambda@df$fold
m_groups_tbl <- table(groups_factor(ihw_res_lambda), folds)
ihw_res_lambda_folds <- ihw(sim$pvalue, groups_factor(ihw_res_tmp), .1, m_groups= m_groups_tbl,
                     folds=folds, lambda=10)
ihw_res_lambda_folds2 <- ihw(sim$pvalue, groups_factor(ihw_res_tmp), .1,
                     folds=folds, lambda=10)
expect_equal(rejections(ihw_res_lambda), rejections(ihw_res_lambda_folds))
expect_equal(weights(ihw_res_lambda_folds), weights(ihw_res_lambda_folds2))

