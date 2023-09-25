#' @export
ihw <- function(...)
{
  UseMethod("ihw")
}




#' ihw: Main function for Independent Hypothesis Weighting
#'
#' Given a vector of p-values, a vector of covariates which are independent of the p-values under the null hypothesis and
#' a nominal significance level alpha, IHW learns multiple testing weights and then applies the weighted Benjamini Hochberg
#' (or Bonferroni) procedure.
#'
#'
#' @param pvalues  Numeric vector of unadjusted p-values.
#' @param covariates  Vector which contains the one-dimensional covariates (independent under the H0 of the p-value)
#'                for each test. Can be numeric or a factor. (If numeric it will be converted into factor by binning.)
#' @param alpha   Numeric, sets the nominal level for FDR control.
#' @param covariate_type  "ordinal" or "nominal" (i.e. whether covariates can be sorted in increasing order or not)
#' @param nbins  Integer, number of groups into which p-values will be split based on covariate. Use "auto" for
#'             automatic selection of the number of bins. Only applicable when covariates is not a factor.
#' @param m_groups Integer vector of length equal to the number of levels of the covariates (only to be specified
#'               when the latter is a factor/categorical). Each entry corresponds to the number of hypotheses to be tested in
#'               each group (stratum). This argument needs to be given when the complete vector of p-values is
#'               not available, but only p-values below a given threshold, for example because of memory reasons.
#'               See the vignette for additional details and an example of how this principle can be applied with
#'               numerical covariates.
#' @param folds  Integer vector or NULL. Pre-specify assignment of hypotheses into folds.
#' @param quiet  Boolean, if False a lot of messages are printed during the fitting stages.
#' @param nfolds Number of folds into which the p-values will be split for the pre-validation procedure
#' @param nfolds_internal  Within each fold, a second  (nested) layer of cross-validation can be conducted to choose a good
#'              regularization parameter. This parameter controls the number of nested folds.
#' @param lambdas  Numeric vector which defines the grid of possible regularization parameters.
#'				Use "auto" for automatic selection.
#' @param seed Integer or NULL. Split of hypotheses into folds is done randomly. To have the output of the function be reproducible,
#'	the seed of the random number generator is set to this value at the start of the function. Use NULL if you don't want to set the seed.
#' @param adjustment_type Character ("BH" or "bonferroni") depending on whether you want to control FDR or FWER.
#' @param null_proportion Boolean, if True (default is False), a modified version of Storey's estimator
#'        is used within each bin to estimate the proportion of null hypotheses.
#' @param null_proportion_level Numeric, threshold for Storey's pi0 estimation procedure, defaults to 0.5
#' @param return_internal Returns a lower level representation of the output (only useful for debugging purposes).
#' @param ... Arguments passed to internal functions.
#'
#' @return A ihwResult object.
#' @seealso ihwResult, plot,ihwResult-method, ihw.DESeqResults
#'
#' @examples
#'
#' save.seed <- .Random.seed; set.seed(1)
#' X   <- runif(20000, min=0, max=2.5)   # covariate
#' H   <- rbinom(20000,1,0.1)            # hypothesis true or false
#' Z   <- rnorm(20000, H*X)              # Z-score
#' .Random.seed <- save.seed
#' pvalue <- 1-pnorm(Z)                  # pvalue
#'
#' ihw_fdr <- ihw(pvalue, X, .1)        # Standard IHW for FDR control
#' ihw_fwer <- ihw(pvalue, X, .1, adjustment_type = "bonferroni")    # FWER control
#
#' table(H[adj_pvalues(ihw_fdr) <= 0.1] == 0) #how many false rejections?
#' table(H[adj_pvalues(ihw_fwer) <= 0.1] == 0)
#'
#'
#' @export
#' @aliases ihw
#'

ihw.default <- function(pvalues,
                        covariates,
                        alpha,
                        covariate_type = "ordinal",
                        nbins = "auto",
                        m_groups = NULL,
                        folds = NULL,
                        quiet = TRUE ,
                        nfolds = 5L,
                        nfolds_internal = nfolds - 1L,
                        lambdas = "auto",
                        seed = 1L,
                        adjustment_type = "BH",
                        null_proportion = FALSE,
                        null_proportion_level = 0.5,
                        return_internal = FALSE,
                        ...) {
  # This function essentially wraps the lower level function ihw_internal
  # e.g. takes care of NAs, sorts pvalues and then
  # returns a nice ihw object



  if (!adjustment_type %in% c("BH", "bonferroni")) {
    stop("IHW currently only works with BH or bonferroni types of multiple testing corrections")
  }


  if (is.null(folds)) {
    nfolds <- as.integer(nfolds)
  } else {
    nfolds <- length(unique(folds))
  }

  nfolds_internal <- as.integer(nfolds_internal)


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
  covariates <- covariates[nna]
  weights <- rep(NA, length(pvalues))

  if (any(is.na(covariates))) {
    stop("Covariates corresponding to non-NA p-values should never be NA. Aborting.")
  }


  if (covariate_type == "ordinal" & is.numeric(covariates)) {
    if (!is.null(m_groups)) {
      stop("m_groups should only be specified when the covariates are a factor")
    }

    if (nbins == "auto") {
      nbins <-
        max(1, min(200, floor(length(pvalues) / 700))) # rule of thumb..
    }
    groups <-
      as.factor(groups_by_filter(covariates, nbins, seed = seed))
    penalty <- "total variation"

  } else if (is.factor(covariates)) {
    groups <- droplevels(covariates)
    if (nbins != "auto" & nbins != nlevels(groups)) {
      warning(
        "Overwriting manually specified nbins, since it has to equal the number
				of levels for categorical covariates"
      )
    }

    nbins <- nlevels(groups)

    if (covariate_type == "nominal") {
      penalty <- "uniform deviation"
    } else if (covariate_type == "ordinal") {
      penalty <- "total variation"
    }
  } else {
    stop("Covariates are not of the appropriate type")
  }

  if ((length(lambdas) == 1) & (lambdas[1] == "auto")) {
    lambdas <- c(Inf, 1.0 / 2 ^ (1:13), 0)
  }

  if (nbins < 1) {
    stop("Cannot have less than one bin.")
  }


  group_levels <- levels(droplevels(groups))
  # sort pvalues globally
  order_pvalues <-
    order(pvalues) #TODO BREAK ties by filter statistic rank
  reorder_pvalues <- order(order_pvalues)

  sorted_groups <- groups[order_pvalues]
  sorted_pvalues <- pvalues[order_pvalues]

  if (!is.null(folds)) {
    sorted_folds <- as.factor(folds[order_pvalues])
  } else {
    sorted_folds <- NULL
  }

  if (is.null(m_groups)) {
    if (is.null(folds)) {
      m_groups <- table(sorted_groups)
    } else {
      m_groups <- table(sorted_groups, sorted_folds)
    }
  } else {
    # TODO: also check if m_groups entered by user is plausible based on observed m_groups and max p_value
  }

  #--- check if m_groups is correctly specified----------------------
  if (is.null(folds)) {
    if (length(group_levels) != length(m_groups)) {
      stop("Number of unique levels should be equal to length of m_groups.")
    }
  } else {
    if (any(dim(m_groups) != c(length(group_levels), nfolds))) {
      stop("m_groups should be of dimension nbins*nfolds")
    }
  }
  # once we have groups, check whether they include enough p-values
  if (nbins > 1 & any(m_groups < 200)) {
    message(
      "We recommend that you supply (many) more than 200 p-values for meaningful data-driven hypothesis weighting results."
    )
  }

  # start by making sure seed is correctly specified
  if (!is.null(seed)) {
    #http://stackoverflow.com/questions/14324096/setting-seed-locally-not-globally-in-r?rq=1
    tmp <- runif(1)
    old <- .Random.seed
    on.exit({
      .Random.seed <<- old
    })
    set.seed(as.integer(seed)) #seed
  }

  if (nbins == 1) {
    nfolds <- 1L
    nfolds_internal <- 1L
    message("Only 1 bin; IHW reduces to Benjamini Hochberg (uniform weights)")
    sorted_adj_p <-
      p.adjust(sorted_pvalues,
               method = adjustment_type,
               n = sum(m_groups))
    rjs <- sum(sorted_adj_p <= alpha)
    res <-
      list(
        lambda = 0,
        fold_lambdas = 0,
        rjs = rjs,
        sorted_pvalues = sorted_pvalues,
        sorted_weighted_pvalues = sorted_pvalues,
        sorted_adj_p = sorted_adj_p,
        sorted_weights = rep(1, length(sorted_pvalues)),
        sorted_groups = as.factor(rep(1, length(sorted_pvalues))),
        sorted_folds = as.factor(rep(1, length(sorted_pvalues))),
        weight_matrix = matrix(1)
      )
  } else if (nbins > 1) {
    res <- ihw_internal(
      sorted_groups,
      sorted_pvalues,
      alpha,
      lambdas,
      m_groups,
      penalty = penalty,
      quiet = quiet,
      sorted_folds = sorted_folds,
      nfolds = nfolds,
      nfolds_internal = nfolds_internal,
      seed = NULL,
      adjustment_type = adjustment_type,
      null_proportion = null_proportion,
      null_proportion_level = null_proportion_level,
      ...
    )
  }

  if (return_internal) {
    return(res)
  }

  # resort back to original form
  weights <-
    fill_nas_reorder(res$sorted_weights, nna, reorder_pvalues)
  pvalues <-
    fill_nas_reorder(res$sorted_pvalues, nna, reorder_pvalues)
  weighted_pvalues <-
    fill_nas_reorder(res$sorted_weighted_pvalues, nna, reorder_pvalues)
  adj_pvalues <-
    fill_nas_reorder(res$sorted_adj_p, nna, reorder_pvalues)
  groups     <-
    factor(fill_nas_reorder(res$sorted_groups, nna, reorder_pvalues),
           levels = group_levels)
  fold_levels <- levels(res$sorted_folds)
  folds      <-
    factor(fill_nas_reorder(res$sorted_folds, nna, reorder_pvalues),
           levels = fold_levels)
  covariates <-
    fill_nas_reorder(covariates, nna, 1:length(covariates))


  df <- data.frame(
    pvalue = pvalues,
    adj_pvalue = adj_pvalues,
    weight = weights,
    weighted_pvalue = weighted_pvalues,
    group = groups,
    covariate = covariates,
    fold = as.factor(folds)
  )

  ihw_obj <- new(
    "ihwResult",
    df = df,
    weights = res$weight_matrix,
    alpha = alpha,
    nbins = as.integer(nbins),
    nfolds = nfolds,
    regularization_term = res$fold_lambdas,
    m_groups = as.integer(m_groups),
    penalty = penalty,
    covariate_type = covariate_type,
    adjustment_type = adjustment_type,
    reg_path_information = data.frame(),
    solver_information = list()
  )

  ihw_obj

}

# operate on pvalues which have already been sorted and without any NAs
ihw_internal <-
  function(sorted_groups,
           sorted_pvalues,
           alpha,
           lambdas,
           m_groups,
           penalty = "total variation",
           return_num_rejection_list = FALSE,
           seed = NULL,
           quiet = TRUE,
           sorted_folds = NULL,
           nfolds = 10L,
           nfolds_internal = nfolds,
           adjustment_type = "BH",
           null_proportion = FALSE,
           null_proportion_level = 0.5,
           debug_flag = FALSE,
           ...) {
    folds_prespecified <- !is.null(sorted_folds)

    n_lambdas <- length(lambdas)
    m <- length(sorted_pvalues)
    split_sorted_pvalues <- split(sorted_pvalues, sorted_groups)
    m_groups_available <- sapply(split_sorted_pvalues, length)


    #  do the k-fold strategy
    if (!is.null(seed)) {
      #http://stackoverflow.com/questions/14324096/setting-seed-locally-not-globally-in-r?rq=1
      old <- .Random.seed
      on.exit({
        .Random.seed <<- old
      })
      set.seed(as.integer(seed)) #seed
    }

    if (!folds_prespecified) {
      sorted_folds <-
        factor(sample(1:nfolds, m, replace = TRUE), levels = 1:nfolds)
    }

    sorted_weights <- rep(NA, m)
    fold_lambdas <- rep(NA, nfolds)
    weight_matrix <- matrix(NA, nlevels(sorted_groups), nfolds)


    ws_mat_list <- list()
    rjs_mat_list <- list()
    for (i in seq_along(levels(sorted_folds))) {
      if (!quiet)
        message(paste("Estimating weights for fold", i))
      if (nfolds == 1) {
        filtered_sorted_groups <- sorted_groups
        filtered_sorted_pvalues <- sorted_pvalues
        filtered_split_sorted_pvalues <-
          split(filtered_sorted_pvalues, filtered_sorted_groups)

        m_groups_holdout_fold <- m_groups
        m_groups_other_folds <- m_groups
      } else {
        fold_i = levels(sorted_folds)[i]
        filtered_sorted_groups <- sorted_groups[sorted_folds != fold_i]
        filtered_sorted_pvalues <- sorted_pvalues[sorted_folds != fold_i]
        filtered_split_sorted_pvalues <-
          split(filtered_sorted_pvalues, filtered_sorted_groups)

        if (!folds_prespecified || length(dim(m_groups)) == 1) {
          # TODO: add check to see if for common use case the first term is 0
          m_groups_holdout_fold <- (m_groups - m_groups_available) / nfolds +
            m_groups_available - sapply(filtered_split_sorted_pvalues, length)
          m_groups_other_folds <- m_groups - m_groups_holdout_fold
          m_groups_other_folds_uncollapsed <- m_groups_other_folds
        } else {
          m_groups_holdout_fold <- m_groups[, colnames(m_groups) == fold_i]
          m_groups_other_folds_uncollapsed <-
            m_groups[, colnames(m_groups) != fold_i, drop = FALSE]
          m_groups_other_folds <-
            rowSums(m_groups_other_folds_uncollapsed)
        }
      }
      ws_mat_list[[i]] <-
        ihw_convex_weights(
          filtered_split_sorted_pvalues,
          alpha,
          m_groups_holdout_fold,
          m_groups_other_folds,
          lambdas,
          penalty = penalty,
          adjustment_type = adjustment_type,
          quiet = quiet,
          ...
        )

      if (!return_num_rejection_list && length(lambdas) > 1) {
        nsplits_internal <- 1L
        rjs_mat  <- matrix(NA, nrow = n_lambdas, ncol = nsplits_internal)
        for (s in seq_len(nsplits_internal)) {
          if (nfolds > 2) {
            sorted_folds_internal = droplevels(sorted_folds[sorted_folds != fold_i])
          } else {
            sorted_folds_internal = NULL
          }
          rjs_int <-
            ihw_internal(
              filtered_sorted_groups,
              filtered_sorted_pvalues,
              alpha,
              lambdas,
              m_groups_other_folds_uncollapsed,
              return_num_rejection_list = TRUE,
              seed = NULL,
              quiet = quiet,
              sorted_folds = sorted_folds_internal,
              nfolds = nfolds_internal,
              adjustment_type = adjustment_type,
              ...
            )
          rjs_mat[, s] <- rjs_int
        }
        rjs_mat_list[[i]] <- rjs_mat
      }
    }
    # within each fold do iterations to also find lambda
    if (return_num_rejection_list) {
      rjs <- rep(NA, n_lambdas)
      for (l in seq_len(n_lambdas)) {
        for (i in seq_along(levels(sorted_folds))) {
          fold_i = levels(sorted_folds)[i]
          sorted_weights[sorted_folds == fold_i] <-
            ws_mat_list[[i]][sorted_groups[sorted_folds == fold_i], l]
        }
        sorted_weighted_pvalues <-
          mydiv(sorted_pvalues, sorted_weights)
        sorted_adj_p <-
          p.adjust(sorted_weighted_pvalues,
                   method = adjustment_type,
                   n = sum(m_groups))
        rjs[l]   <- sum(sorted_adj_p <= alpha)
      }
      return(rjs)
    }

    if (!return_num_rejection_list) {
      for (i in seq_along(levels(sorted_folds))) {
        fold_i = levels(sorted_folds)[i]
        if (length(lambdas) > 1) {
          lambda_idx <- which.max(apply(rjs_mat_list[[i]], 1, mean))
        } else {
          lambda_idx <- 1
        }
        weight_matrix[, i] <- ws_mat_list[[i]][, lambda_idx]
        sorted_weights[sorted_folds == fold_i] <-
          ws_mat_list[[i]][sorted_groups[sorted_folds == fold_i], lambda_idx]
        fold_lambdas[i] <- lambdas[lambda_idx]
      }
      sorted_weighted_pvalues <-
        mydiv(sorted_pvalues, sorted_weights)
      sorted_adj_p <-
        p.adjust(sorted_weighted_pvalues,
                 method = adjustment_type,
                 n = sum(m_groups))
      rjs   <- sum(sorted_adj_p <= alpha)
      lst <-
        list(
          lambda = 0.0,
          fold_lambdas = fold_lambdas,
          rjs = rjs,
          sorted_pvalues = sorted_pvalues,
          sorted_weighted_pvalues = sorted_weighted_pvalues,
          sorted_adj_p = sorted_adj_p,
          sorted_weights = sorted_weights,
          sorted_groups = sorted_groups,
          sorted_folds = sorted_folds,
          weight_matrix = weight_matrix
        )
      return(lst)
    }


    #	lambda_idx <- which.max(apply(rjs,1, mean))
    # fold_lambdas[i] <- lambda



    #if (null_proportion){
    #	pi0_est <- weighted_storey_pi0(sorted_pvalues[sorted_folds==i],
    #		                  sorted_weights[sorted_folds == i],
    #		                  tau= null_proportion_level,
    #		                  m=sum(m_groups_holdout_fold))

    #	sorted_weights[sorted_folds == i] <- sorted_weights[sorted_folds == i]/pi0_est
    #	weight_matrix[,i] <- weight_matrix[,i]/pi0_est

    #}
  }


#' @rdname ihw.default
#' @param formula \code{\link{formula}}, specified in the form pvalue~covariate (only 1D covariate supported)
#' @param data data.frame from which the variables in formula should be taken
#' @export
ihw.formula <- function(formula, data = parent.frame(), ...) {
  if (length(formula) != 3) {
    stop("expecting formula of the form pvalue~covariate")
  }
  pv_name <- formula[[2]]
  cov_name <- formula[[3]]
  pvalues <- eval(pv_name, data, parent.frame())
  covariates <- eval(cov_name, data, parent.frame())
  ihw(pvalues, covariates, ...)
}


#' ihw.DESeqResults: IHW method dispatching on DESeqResults objects
#'
#' @param deseq_res "DESeqResults" object
#' @param filter Vector of length equal to number of rows of deseq_res object. This is used
#'         for the covariates in the call to ihw. Can also be a character,
#'         in which case deseq_res[[filter]] is used as the covariate
#' @param alpha   Numeric, sets the nominal level for FDR control.
#' @param adjustment_type Character ("BH" or "bonferroni") depending on whether you want to control FDR or FWER.
#' @param ... Other optional keyword arguments passed to ihw.
#'
#' @return A "DESeqResults" object, which includes weights and adjusted p-values returned
#'         	by IHW. In addition, includes a metadata slot with an "ihwResult" object.
#' @seealso ihw, ihwResult
#'
#' @examples \dontrun{
#'    library("DESeq2")
#'    library("airway")
#'    data("airway")
#'    dds <- DESeqDataSet(se = airway, design = ~ cell + dex)
#'    dds <- DESeq(dds)
#'    deseq_res <- results(dds)
#'    deseq_res <- ihw(deseq_res, alpha=0.1)
#'    #equivalent: deseq_res2 <- results(dds, filterFun = ihw)
#' }
#'
#' @export
#'
ihw.DESeqResults <- function(deseq_res,
                             filter = "baseMean",
                             alpha = 0.1,
                             adjustment_type = "BH",
                             ...) {
  if (missing(filter)) {
    filter <- deseq_res$baseMean
  } else if (is.character(filter)) {
    stopifnot(filter %in% colnames(deseq_res))
    filter <- deseq_res[[filter]]
  }

  stopifnot(length(filter) == nrow(deseq_res))

  ihw_res <- ihw(deseq_res$pvalue,
                 filter,
                 alpha = alpha,
                 adjustment_type = adjustment_type,
                 ...)

  deseq_res$padj <- adj_pvalues(ihw_res)
  deseq_res$weight <- weights(ihw_res)

  mcols(deseq_res)$type[names(deseq_res) == "padj"] <- "results"
  mcols(deseq_res)$description[names(deseq_res) == "padj"] <-
    paste("Weighted", adjustment_type, "adjusted p-values")

  mcols(deseq_res)$type[names(deseq_res) == "weight"] <- "results"
  mcols(deseq_res)$description[names(deseq_res) == "weight"] <-
    "IHW weights"

  metadata(deseq_res)[["alpha"]] <- alpha
  metadata(deseq_res)[["ihwResult"]] <- ihw_res
  deseq_res
}

ihw_convex_weights <-
  function(split_sorted_pvalues,
           alpha,
           m_groups,
           m_groups_grenander,
           penalty = "total variation",
           lambdas,
           adjustment_type = "BH",
           grenander_binsize = 1,
           quiet = quiet) {
    # m_groups used for balancing (weight budget)
    # m_groups_grenander (used for grenander estimator)
    # Practically always it will hold: m_groups \approx m_groups_grenander

    nbins <- length(split_sorted_pvalues)

    if (length(lambdas) == 1 && lambdas[1] == Inf) {
      ws <- matrix(1, nbins, 1)
      return(ws)
    }

    zero_included <- any(lambdas == 0)
    infty_included <- any(is.infinite(lambdas))

    lambdas_filt <- lambdas[!(lambdas %in% c(0, Inf))]
    # also special case if lambdas == 0

    # also there temporarily to help with numerical stability of downstream optimization
    max_pval <- 0.9 * max(sapply(split_sorted_pvalues, max))
    if (max_pval < 0.1) {
      split_sorted_pvalues <-
        lapply(split_sorted_pvalues, function(x)
          x / max_pval)
      alpha <- alpha / max_pval
    }
    # preprocessing:  Set very low p-values to 10^(-20) to help with numerical stability
    # Note: This only affects the internal optimization, the higher level functions will
    # still return adjusted pvalues based on the original p-values
    split_sorted_pvalues <-
      lapply(split_sorted_pvalues, function(x)
        ifelse(x > 10 ^ (-20), x, 10 ^ (-20)))
    m <- sum(m_groups)

    if (nbins != length(m_groups)) {
      stop("length of m_groups should be equal to number of bins")
    }

    #lapply grenander...
    if (!quiet)
      message("Applying Grenander estimator within each bin.")
    grenander_list <- mapply(
      presorted_grenander,
      split_sorted_pvalues,
      m_groups_grenander,
      grenander_binsize = grenander_binsize,
      quiet = quiet,
      SIMPLIFY = FALSE
    )

    if (length(lambdas_filt) == 0 && zero_included) {
      ts <-  unregularized_thresholds_bh(grenander_list, m_groups, alpha)
      if (infty_included) {
        ts <- cbind(1, ts)
      } else {
        ts <- matrix(ts, ncol = 1)
      }
    } else {
      ts <-  optimal_ts(grenander_list, m_groups, alpha, lambdas_filt)

      if (!zero_included) {
        ts <- ts[,-ncol(ts), drop = FALSE]
      }
      if (!infty_included) {
        ts <- ts[,-1, drop = FALSE]
      }
      ts <- pmax(ts, 0)
    }

    ws <- apply(ts, 2, function(x)
      thresholds_to_weights(x, m_groups))

    # return matrix as follows: nbins x nlambdas)
    return(ws)
  }


#' @importFrom fdrtool gcmlcm
presorted_grenander <-
  function(sorted_pvalues,
           m_total = length(sorted_pvalues),
           grenander_binsize = 1,
           quiet = TRUE) {
    unique_pvalues <- unique(sorted_pvalues)
    ecdf_values <-
      cumsum(tabulate(match(sorted_pvalues, unique_pvalues))) / m_total

    # I think fdrtool neglects this borderline case and this causes returned object
    # to be discontinuous hence also not concave
    unique_pvalues <- c(0, unique_pvalues)
    ecdf_values   <- c(0, ecdf_values)


    if (grenander_binsize != 1) {
      nmax = length(unique_pvalues)
      idx_thin <- c(seq(1, nmax - 1, by = grenander_binsize), nmax)
      unique_pvalues <- unique_pvalues[idx_thin]
      ecdf_values <- ecdf_values[idx_thin]
    }

    ll <- fdrtool::gcmlcm(unique_pvalues, ecdf_values, type = "lcm")

    gren <-
      list(
        locs = ll$x.knots,
        Fs = ll$y.knots,
        fs = ll$slope.knots
      )
    if (!quiet)
      message(paste("Grenander fit with", length(ll$slope.knots), "knots."))
    gren
  }
