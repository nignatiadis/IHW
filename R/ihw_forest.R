# for internal use only
# similar to ihw.default, but for using stratification method "forest"
# because groups in ihw.default_forest has the structure of nested list instead of
# a simple vector as in ihw.default, a separate function was necessary
ihw_forest <- function(pvalues, covariates, alpha,
                       covariate_type,
                       nbins,
                       m_groups,
                       folds,
                       quiet,
                       nfolds,
                       nfolds_internal,
                       nsplits_internal,
                       lambdas,
                       seed,
                       distrib_estimator,
                       lp_solver,
                       adjustment_type,
                       null_proportion,
                       null_proportion_level,
                       return_internal,
                       ntrees,
                       tau,
                       nodedepth,
                       nodesize,
                       mtry,
                       ...) {

  # This function essentially wraps the lower level function ihw_internal
  # e.g. takes care of NAs, sorts pvalues and then
  # returns a nice ihw object
  if (!(covariate_type == "ordinal" & is.numeric(covariates))) {
    stop("the stratificaion_method forest only works with covariate_type = ordinal")
  }

  if (!is.null(m_groups)) {
    stop("m_groups should only be specified when the covariates are a factor")
  }

  if (nbins != "auto") {
    stop("with the stratificaion_method forest nbins can not be directly controlled. Please use the parameters nodedepth and nodesize instead")
  }

  if ((length(lambdas) == 1) & (lambdas[1] == "auto")) {
    nbins <- max(1, min(40, floor(length(pvalues) / 1500))) # rule of thumb..
    # default for forest without regularisation
    lambdas <- Inf
  }

  if (nfolds == 1) {
    message("Using only 1 fold! Only use this if you want to learn the weights, but NEVER for testing!")
  }

  # gracefully handle NAs
  nna <- !is.na(pvalues)

  if (all(nna)) {
    nna <- TRUE
  }

  # filter our p-values
  pvalues <- pvalues[nna]
  weights <- rep(NA, length(pvalues))
  # browser()
  covariates <- as.matrix(covariates) # cannot convert factor to matrix
  covariates <- covariates[nna, , drop = FALSE]

  # rename columns
  nvar <- ncol(covariates)
  if (nvar == 1) {
    colnames(covariates) <- "covariate"
  } else {
    colnames(covariates) <- seq_len(nvar)
  }

  if (any(is.na(covariates))) {
    stop("Covariates corresponding to non-NA p-values should never be NA. Aborting.")
  }

  if (is.null(folds)) {
    nfolds <- as.integer(nfolds)
    # for the forest sratification, the folds need to be pre-specified, because the construction depends on the p-values
    folds <- sample(seq_len(nfolds), length(pvalues), replace = TRUE)
  } else {
    nfolds <- length(unique(folds))
  }

  # start by making sure seed is correctly specified
  if (!is.null(seed)) {
    # http://stackoverflow.com/questions/14324096/setting-seed-locally-not-globally-in-r?rq=1
    tmp <- stats::runif(1)
    old <- .Random.seed
    on.exit({
      .Random.seed <<- old
    })
    set.seed(as.integer(seed)) # seed
  }

  # the stratification method is the main point of deviation of ihw.default_forest from ihw.default
  # groups in ihw.default_forest has the structure of nested list compared to a simple vector in ihw.default
  groups <- group_by_forest(pvalues, covariates, folds, ntrees, tau, nodedepth, nodesize, mtry, seed)

  penalty <- "uniform deviation"

  group_levels <- lapply(groups, function(groups_i) {
    lapply(groups_i, function(groups_i_t) levels(droplevels(groups_i_t)))
  })
  # sort pvalues globally
  order_pvalues <- order(pvalues) # TODO BREAK ties by filter statistic rank
  reorder_pvalues <- order(order_pvalues)

  sorted_groups <- lapply(groups, function(groups_i) {
    groups_i[order_pvalues, , drop = FALSE]
  })
  sorted_pvalues <- pvalues[order_pvalues]

  if (!is.null(folds)) {
    sorted_folds <- folds[order_pvalues]
  } else {
    sorted_folds <- NULL
  }

  m_groups <- lapply(sorted_groups, function(sorted_groups_i) {
    lapply(sorted_groups_i, function(sorted_groups_i_t) table(sorted_groups_i_t, sorted_folds))
  })

  # run diagnostic and give feedback
  group_levels_length <- lapply(group_levels, function(group_levels_i) {
    lapply(group_levels_i, function(group_levels_i_t) length(group_levels_i_t))
  })
  group_levels_length <- unlist(group_levels_length)
  nbins_ideal <- max(1, min(40, floor(length(pvalues) / 1500))) # rule of thumb..
  if (all(group_levels_length) < nbins_ideal) {
    message(paste("For", length(pvalues), "non-na pvalues, we recommend", nbins_ideal, " stratification bins for granular hypothesis weighting. Consider to increase nodedepth or decrease nodesize."))
  }

  # once we have groups, check whether they include enough p-values
  m_groups_unlist <- unlist(m_groups)
  if (all(m_groups_unlist < 2)) {
    stop("Bins of size < 2 are currently not allowed. Please decrease nodedepth or increase nodesize.")
  } else if (any(m_groups_unlist < 1000)) {
    message("We recommend that you supply (many) more than 1000 p-values for meaningful data-driven hypothesis weighting results. Consider to decrease nodedepth or increase nodesize.")
  }
  ntrees_all_comb <- ntrees # length(res[[1]])

  # loop over all folds and trees and run ihw_internal
  res <- lapply(seq_len(nfolds), function(i) {
    lapply(seq_len(ntrees_all_comb), function(t) {
      sorted_groups_i_t <- sorted_groups[[i]][, t, drop = TRUE]
      m_groups_i_t <- m_groups[[i]][[t]]
      ihw_internal(sorted_groups_i_t, sorted_pvalues, alpha, lambdas,
        m_groups_i_t,
        penalty = penalty,
        quiet = quiet,
        sorted_folds = sorted_folds,
        nfolds = nfolds,
        nfolds_internal = nfolds_internal,
        nsplits_internal = nsplits_internal,
        seed = NULL,
        distrib_estimator = distrib_estimator,
        lp_solver = lp_solver,
        adjustment_type = adjustment_type,
        null_proportion = null_proportion,
        null_proportion_level = null_proportion_level,
        folds_subset = i, # only run on the relevant fold
        ...
      )
    })
  })

  if (return_internal) {
    return(res)
  }

  # extract lambdas, weights from nested list
  ntrees_all_comb <- ntrees # length(res[[1]])

  # extract lambdas
  fold_lambdas <- matrix(NA, nrow = nfolds, ncol = ntrees_all_comb)
  for (i in seq_len(nfolds)) {
    for (t in seq_len(ntrees_all_comb)) {
      fold_lambdas[i, t] <- res[[i]][[t]]$lambda
    }
  }

  groups <- do.call(cbind, groups)

  # extract sorted weights
  sorted_weights <- matrix(NA, nrow = length(pvalues), ncol = nfolds * ntrees_all_comb)
  for (i in seq_len(nfolds)) {
    for (t in seq_len(ntrees_all_comb)) {
      sorted_weights[, (i - 1) * ntrees_all_comb + t] <- res[[i]][[t]]$sorted_weights
    }
  }

  weight_matrices_forest <- lapply(
    seq_len(nfolds),
    function(i) {
      lapply(
        seq_len(ntrees_all_comb),
        function(t) res[[i]][[t]]$weight_matrix[, i, drop = TRUE]
      )
    }
  )

  # average out individual trees (random forest approach)
  sorted_weights <- rowMeans(sorted_weights, na.rm = TRUE)
  sorted_weighted_pvalues <- mydiv(sorted_pvalues, sorted_weights)

  sorted_adj_p <- p.adjust(sorted_weighted_pvalues, method = adjustment_type, n = length(sorted_weighted_pvalues))

  # resort back to original form
  weights <- fill_nas_reorder(sorted_weights, nna, reorder_pvalues)
  pvalues <- fill_nas_reorder(sorted_pvalues, nna, reorder_pvalues)
  weighted_pvalues <- fill_nas_reorder(sorted_weighted_pvalues, nna, reorder_pvalues)
  adj_pvalues <- fill_nas_reorder(sorted_adj_p, nna, reorder_pvalues)
  folds <- factor(fill_nas_reorder(sorted_folds, nna, reorder_pvalues), levels = seq_len(nfolds))

  groups <- fill_nas_reorder_dataframe(groups, nna, seq_len(nrow(groups)))
  if (is.factor(covariates)) {
    covariates <- fill_nas_reorder(covariates, nna, seq_along(covariates))
  } else if (is.matrix(covariates)) {
    covariates <- fill_nas_reorder_dataframe(covariates, nna, seq_len(nrow(covariates)))
  }

  df <- data.frame(
    pvalue = pvalues,
    adj_pvalue = adj_pvalues,
    weight = weights,
    weighted_pvalue = weighted_pvalues,
    covariate = covariates,
    fold = as.factor(folds)
  )
  df <- cbind(df, groups)

  ihw_obj <- new("ihwResult",
    df = df,
    weights = matrix(NA), # there is no single corresponding weight matrix for with forest stratification
    alpha = alpha,
    nbins = integer(0), # there is no corresponding nbins for with forest stratification
    nfolds = nfolds,
    regularization_term = fold_lambdas,
    m_groups = integer(0), # there is no corresponding m_groups for with forest stratification
    penalty = penalty,
    covariate_type = covariate_type,
    adjustment_type = adjustment_type,
    reg_path_information = data.frame(),
    solver_information = list(),
    stratification_method = "forest",
    weight_matrices_forest = weight_matrices_forest,
    tau = as.double(tau)
  )

  ihw_obj
}

#' Stratify hypotheses based on random forest construction
#'
#' Hypotheses are stratified into bins based on random forest construction, alternative to \code{groups_by_filter_multivariate}
#'   groups are homogenous wrt to Storeys null proportion estimator
#'
#'  see https://doi.org/10.7717/peerj.6035 for details on BocaLeek construction
#' @param pvalues Numeric vector of unadjusted p-values.
#' @param covariates Matrix which contains the covariates (independent under the H0 of the p-value) for each test.
#' @param folds Integer vector, Pre-specify assignment of hypotheses into folds.
#' @param ntrees Integer, see same parameter in \code{\link[randomForestSRC]{rfsrc}}
#' @param tau Double, censoring threshold tau of the pvalues in the stratification method "forest". See more in group_by_forest
#' @param nodedepth Integer, see same parameter in \code{\link[randomForestSRC]{rfsrc}}
#' @param nodesize Integer, see same parameter in \code{\link[randomForestSRC]{rfsrc}}
#' @param mtry Integer, see same parameter in \code{\link[randomForestSRC]{rfsrc}}
#'               Use "auto" for automatic selection.
#' @param seed Integer, specifies random seed to be used
#' @importFrom randomForestSRC rfsrc
#' @importFrom stats predict
#' @return list of data.frames representing group membership for each test. list length equals the number of folds.
#' @examples
#'
#' save.seed <- .Random.seed
#' set.seed(1)
#' X <- runif(20000, min = 0, max = 2.5) # covariate
#' H <- rbinom(20000, 1, 0.1) # hypothesis true or false
#' Z <- rnorm(20000, H * X) # Z-score
#' folds <- sample(1:5, 20000, replace = TRUE)
#' .Random.seed <- save.seed
#' pvalue <- 1 - pnorm(Z) # pvalue
#' groups <- group_by_forest(pvalue, as.matrix(X), folds)
#' m_groups <- lapply(groups, function(group) {
#'   lapply(group, function(group_i) table(group_i, folds))
#' })
#' @export
group_by_forest <- function(pvalues, covariates, folds, ntrees = 10, tau = 0.5, nodedepth = 4, nodesize = 1000, mtry = "auto", seed = NULL) {
  m <- length(pvalues)
  nfolds <- length(unique(folds))

  if (mtry == "auto") mtry <- ceiling(0.9 * ncol(covariates)) # a lot of noise data => high mtry

  nodesize <- as.integer(nodesize)
  mtry <- as.integer(mtry)
  nodedepth <- as.integer(nodedepth)


  groups <- lapply(seq_len(nfolds), function(i) {
    # binary indicator from Boca and leek/storey
    data <- data.frame(
      indic = (pvalues >= tau),
      covariates = covariates
    )
    data_other_folds <- data[folds != i, ]

    # grow forest based on other folds
    forest_other_fold <- randomForestSRC::rfsrc(
      indic ~ . - indic,
      data = data_other_folds,
      ntree = ntrees,
      mtry = mtry,
      nodesize = nodesize,
      nodedepth = nodedepth,
      splitrule = "mse",
      block.size = FALSE,
      forest.wt = FALSE,
      seed = seed
    )

    # predict terminal nodes for all covariates based on the forest structure
    predict_groups <- stats::predict(forest_other_fold, data, membership = TRUE)

    groups <- predict_groups$membership
    groups <- as.data.frame(groups)
    groups[] <- lapply(groups, as.factor)

    names(groups) <- paste0("fold", i, "_tree_", seq_along(groups))
    groups
  })

  return(groups)
}
