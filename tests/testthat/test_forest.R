#--------- Helper functions ------------------------------------------#
#-- test that direct extension of helpers functions are indeed equivalent
#--------------------------------------------------------------------#

#groups_by_filter_multivariate extends groups_by_filter
covariates <- runif(300)
groups_vector <- groups_by_filter(covariates, 10, seed = 1)
groups_multivariate <- groups_by_filter_multivariate(matrix(covariates), 10, seed = 1)
expect_true(all(groups_vector == groups_multivariate))

#fill_nas_reorder_dataframe extends fill_nas_reorder
covariates_reorder_vector <- fill_nas_reorder(covariates, rep(T, 300), 1:length(covariates))
covariates_reorder_multivariate <- unlist(fill_nas_reorder_dataframe(as.matrix(covariates), rep(T, 300), 1:length(covariates)))
expect_true(all(covariates_reorder_vector == covariates_reorder_multivariate))


# regular simulation
wasserman_normal_sim <- function(m, pi0, xi_min, xi_max, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X <- runif(m, min = xi_min, max = xi_max)
  H <- rbinom(m, 1, 1 - pi0)
  Z <- rnorm(m, H * X)
  pvalue <- 1 - pnorm(Z)
  simDf <- data.frame(pvalue = pvalue, filterstat = X, H = H, Z = Z, noise = 0)
}
m <- 10000
sim <- wasserman_normal_sim(m, 0.85, 0, 3, seed = 1)

#create covariates of all above formats 
cov_num_vec <- sim$filterstat
cov_fac_vec <- IHW:::groups_by_filter(cov_num_vec, 4)
cov_num_mat <- as.matrix(cov_num_vec)
cov_num_mat_noise <- as.matrix(sim[, c("filterstat", "noise")])

#--------- IHW general ------------------------------------------#
#--- should reduce to BH---------------
#--------------------------------------------------------------------#
ihw_test1 <- ihw(sim$pvalue, cov_num_vec, .1, nbins = 1)
ihw_test2 <- ihw(sim$pvalue, cov_num_vec, .1, lambda = 0)
rej_bh <- sum(p.adjust(sim$pvalue, n = m, method = "BH") <= .1)

expect_equal(rejections(ihw_test1), rejections(ihw_test2))
expect_equal(rejections(ihw_test2), rej_bh)

#--------- IHW quantile ------------------------------------------#
# we want IHW to work for covariates input of the class:
# vector-factor, vector-numeric (old)
# matrix-numeric, 1 column and multiple columns (new)
#--------------------------------------------------------------------#

## ---run ihw with quantile stratification with different covariate class, check that yield equivalent results----
set.seed(1)
ihw_num_vec <- ihw(sim$pvalue, cov_num_vec, .1, nbins = 4, stratification_method = "quantiles")
set.seed(1)
ihw_fac_vec <- ihw(sim$pvalue, cov_fac_vec, .1)
set.seed(1)
ihw_num_mat <- ihw(sim$pvalue, cov_num_mat, .1, nbins = 4, stratification_method = "quantiles")

expect_true(all(groups_factor(ihw_num_vec) == groups_factor(ihw_num_mat)))
expect_true(all(ihw_num_vec@df$fold == ihw_num_mat@df$fold))

expect_equal(rejections(ihw_num_vec), rejections(ihw_fac_vec))
expect_equal(rejections(ihw_num_vec), rejections(ihw_num_mat))

#common sensem that it should not reject all hypotheses
expect_lt(rejections(ihw_num_vec), m)
## ---test the formula interface----

ihw_num_mat_noise <- ihw(sim$pvalue, cov_num_mat_noise, .1, stratification_method = "quantiles")
ihw_num_mat_noise_formula1 <- ihw(pvalue ~ filterstat + noise, data = sim, alpha = .1)
ihw_num_mat_noise_formula2 <- ihw(sim$pvalue ~ sim$filterstat + sim$noise, alpha = .1)

expect_equal(rejections(ihw_num_mat_noise), rejections(ihw_num_mat_noise_formula1))
expect_equal(rejections(ihw_num_mat_noise), rejections(ihw_num_mat_noise_formula1))
expect_lt(rejections(ihw_num_mat_noise), m)

#--------- IHW forest ------------------------------------------#
#--- run ihw with forest stratificatio---------------
#--------------------------------------------------------------------#

#set forest parameter
ntrees <- 2
n_censor_thres <- 2

#run ihw forest with different covariate format
set.seed(1)
ihw_forest_num_vec <- ihw(sim$pvalue, cov_num_vec, .1, stratification_method = "forest", ntrees = ntrees, n_censor_thres = n_censor_thres)
set.seed(1)
ihw_forest_num_mat <- ihw(sim$pvalue, cov_num_mat, .1, stratification_method = "forest", ntrees = ntrees, n_censor_thres = n_censor_thres)
set.seed(1)
ihw_forest_num_mat_noise <- ihw(sim$pvalue, cov_num_mat_noise, .1, stratification_method = "forest", ntrees = ntrees, n_censor_thres = n_censor_thres)

#should be equivalent
expect_equal(rejections(ihw_forest_num_vec), rejections(ihw_forest_num_mat)) 
expect_lt(rejections(ihw_forest_num_vec), m)

# check if weight budget fullfilled
expect_equal(sum(weights(ihw_forest_num_vec)), m)

#check if nrow works
expect_equal(m, nrow(ihw_forest_num_vec))

# methods which return vector equal to number of hypotheses
mymethods <- c(adj_pvalues, weights, thresholds, pvalues, weighted_pvalues, rejected_hypotheses) 
lengths <- sapply(mymethods, function(f) length(f(ihw_forest_num_vec)))
expect_true(all(lengths == m))
expect_equal(length(covariates(ihw_forest_num_vec)), m)

lengths <- sapply(mymethods, function(f) length(f(ihw_forest_num_mat_noise)))
expect_true(all(lengths == m))
expect_equal(nrow(covariates(ihw_forest_num_mat_noise)), m) 


#--------- IHW forest getters------------------------------------------#
#--- check IHW class getters---------------
#--------------------------------------------------------------------#

#IHW::groups_factor
expect_true(is.data.frame(groups_factor(ihw_forest_num_vec)))
expect_equal(dim(groups_factor(ihw_forest_num_vec)), c(m, nfolds(ihw_forest_num_vec) *ntrees *n_censor_thres))

#IHW::m_groups
expect_warning(ihw_forest_num_vec_m_groups <- m_groups(ihw_forest_num_vec))
expect_true(is.list(ihw_forest_num_vec_m_groups))
expect_equal(nfolds(ihw_forest_num_vec) *ntrees *n_censor_thres, length(ihw_forest_num_vec_m_groups))
expect_true(is.matrix(ihw_forest_num_vec_m_groups[[1]]))
expect_equal(ncol(ihw_forest_num_vec_m_groups[[1]]), nfolds(ihw_forest_num_vec))

# IHW::weights
ihw_forest_num_vec_weights <- weights(ihw_forest_num_vec, levels_only = T)
expect_true(is.list(ihw_forest_num_vec_weights))
expect_equal(length(ihw_forest_num_vec_weights), nfolds(ihw_forest_num_vec))

# IHW::stratification_breaks
expect_true(is.numeric(stratification_breaks(ihw_forest_num_vec)))

#------ some getters are disabled for forest and should throw errors----
expect_error(thresholds(ihw_forest_num_vec, levels_only = TRUE)) 
expect_error(plugin_fdr(ihw_forest_num_vec))
expect_error(plot(ihw_forest_num_mat_noise))
expect_error(plot(ihw_forest_num_vec))
