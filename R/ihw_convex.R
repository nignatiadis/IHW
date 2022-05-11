#' @export
ihw <- function(...)
{
    UseMethod( "ihw" )
}




#' ihw: Main function for Independent Hypothesis Weighting
#'
#' Given a vector of p-values, a vector of covariates which are independent of the p-values under the null hypothesis and
#' a nominal significance level alpha, IHW learns multiple testing weights and then applies the weighted Benjamini Hochberg
#' (or Bonferroni) procedure.
#'
#'
#' @param pvalues  Numeric vector of unadjusted p-values.
#' @param covariates  Vector or matrix which contains the covariates (independent under the H0 of the p-value)
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
#' @param nsplits_internal  Integer, how many times to repeat the nfolds_internal splitting. Can lead to better regularization
#'              parameter selection but makes ihw a lot slower.
#' @param lambdas  Numeric vector which defines the grid of possible regularization parameters.
#'				Use "auto" for automatic selection.
#' @param seed Integer or NULL. Split of hypotheses into folds is done randomly. To have the output of the function be reproducible, 
#'	the seed of the random number generator is set to this value at the start of the function. Use NULL if you don't want to set the seed.
#' @param distrib_estimator  Character ("grenander" or "ECDF"). Only use this if you know what you are doing. ECDF with nfolds > 1
#'              or lp_solver == "lpsymphony" will in general be excessively slow, except for very small problems.
#' @param lp_solver  Character ("lpsymphony" or "gurobi"). Internally, IHW solves a sequence of linear programs, which
#'        can be solved with either of these solvers.
#' @param adjustment_type Character ("BH" or "bonferroni") depending on whether you want to control FDR or FWER.
#' @param null_proportion Boolean, if True (default is False), a modified version of Storey's estimator
#'        is used within each bin to estimate the proportion of null hypotheses.
#' @param null_proportion_level Numeric, threshold for Storey's pi0 estimation procedure, defaults to 0.5
#' @param return_internal Returns a lower level representation of the output (only useful for debugging purposes).
#' @param stratification_method Character ("quantiles" or "forest" or "none")
#' @param ntrees Integer, see same parameter in \code{\link[randomForestSRC]{rfsrc}} used when \code{stratification_method=="forest"}
#' @param n_censor_thres Integer, number of censoring thresholds tau to be considered for stratification method "forest"
#' @param nodedepth Integer, see same parameter in \code{\link[randomForestSRC]{rfsrc}} used when \code{stratification_method=="forest"}
#' @param nodesize Integer, see same parameter in \code{\link[randomForestSRC]{rfsrc}} used when \code{stratification_method=="forest"}
#' @param mtry Integer, see same parameter in \code{\link[randomForestSRC]{rfsrc}} used when \code{stratification_method=="forest"}
#'				Use "auto" for automatic selection.
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

ihw.default <- function(pvalues, covariates, alpha,
						covariate_type = "ordinal",
						nbins = "auto",
						m_groups = NULL,
						folds = NULL,
						quiet =TRUE ,
						nfolds = 5L,
						nfolds_internal = 5L,
						nsplits_internal=1L,
						lambdas = "auto",
						seed = 1L,
						distrib_estimator = "grenander",
						lp_solver="lpsymphony",
						adjustment_type = "BH",
						null_proportion = FALSE,
						null_proportion_level = 0.5,
						return_internal = FALSE,
						stratification_method = "quantiles",
						ntrees = 5L,
						n_censor_thres = 10L,
						nodedepth = 4,
						nodesize = 1000,
						mtry = "auto",
						...){

	# This function essentially wraps the lower level function ihw_internal
	# e.g. takes care of NAs, sorts pvalues and then
	# returns a nice ihw object



	if (!adjustment_type %in% c("BH","bonferroni")){
		stop("IHW currently only works with BH or bonferroni types of multiple testing corrections")
	}

	if (adjustment_type == "bonferroni" & distrib_estimator != "grenander"){
		stop("For Bonferroni-like FWER adjustment currently only the grenander estimator is supported.")
	}

  if (stratification_method == "forest" & nbins > 1) {
    res <- ihw_forest(
      pvalues, covariates, alpha,
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
      n_censor_thres,
      nodedepth,
      nodesize,
      mtry,
      ...
    )
    return(res)
  }
  
	if (is.null(folds)){
		nfolds <- as.integer(nfolds)
	} else {
		nfolds <- length(unique(folds))
	}

	nfolds_internal <- as.integer(nfolds_internal)


	if (nfolds==1){
		message("Using only 1 fold! Only use this if you want to learn the weights, but NEVER for testing!")
	}


	# gracefully handle NAs
   	nna <- !is.na(pvalues)

   	if (all(nna)){
        nna <- TRUE
   	}

   	# filter our p-values
    pvalues <- pvalues[nna]
    
    if (is.matrix(covariates)) {
      covariates <- covariates[nna, , drop = FALSE]
    } else {
      covariates <- covariates[nna]
    }
    
    weights <- rep(NA, length(pvalues))

    if (any(is.na(covariates))){
    	stop("Covariates corresponding to non-NA p-values should never be NA. Aborting.")
    }


	if (covariate_type =="ordinal" & is.numeric(covariates)){
	  covariates <- as.matrix(covariates)
	  nvar <- ncol(covariates)
	  if (nvar == 1) {
	    colnames(covariates) <- "covariate"
	  } else {
	    colnames(covariates) <- seq_len(nvar)
	  }
	  
		if (!is.null(m_groups)){
			stop("m_groups should only be specified when the covariates are a factor")
		}

		if (nbins == "auto"){
			nbins <- max(1,min(40, floor(length(pvalues)/1500))) # rule of thumb..
		}
	  
	  if (stratification_method == "quantiles") {
	    groups <- groups_by_filter_multivariate(covariates, nbins, seed = seed)
	  } else {
	    stop("Unknown stratification method")
	  }
	  
		penalty <- "total variation"

	} else if (is.factor(covariates)){
	  stratification_method <- "none"
		groups <- covariates
		if (nbins != "auto" & nbins != nlevels(groups)){
			warning("Overwriting manually specified nbins, since it has to equal the number
				of levels for categorical covariates")
		}

		nbins <- nlevels(groups)

		if (covariate_type == "nominal"){
			penalty <- "uniform deviation"
		} else if (covariate_type == "ordinal") {
			penalty <- "total variation"
		}
	} else {
		stop("Covariates are not of the appropriate type")
	}

	if ((length(lambdas)== 1) & (lambdas[1] == "auto")){
		# just a few for now until I get warm starts of the LP solvers to work
		lambdas <- c(0, 1, nbins/8, nbins/4, nbins/2, nbins, Inf)
	}
  
	if (nbins < 1) {
		stop("Cannot have less than one bin.")
	}


	group_levels <- levels(droplevels(groups))
	# sort pvalues globally
	order_pvalues <- order(pvalues) #TODO BREAK ties by filter statistic rank
	reorder_pvalues <- order(order_pvalues)

	sorted_groups <- groups[order_pvalues]
	sorted_pvalues <- pvalues[order_pvalues]

	if (!is.null(folds)){
		sorted_folds <- folds[order_pvalues]
	} else {
		sorted_folds <- NULL
	}

	if (is.null(m_groups)) {
		if (is.null(folds)){
			m_groups <- table(sorted_groups)
		} else {
			m_groups <- table(sorted_groups, sorted_folds)
		}
	} else {
		# TODO: also check if m_groups entered by user is plausible based on observed m_groups and max p_value
	}

    #--- check if m_groups is correctly specified----------------------
    if (is.null(folds)){
		if (length(group_levels) != length(m_groups)){
			stop("Number of unique levels should be equal to length of m_groups.")
		}
	} else {
		if (any(dim(m_groups) != c(length(group_levels), nfolds))){
			stop("m_groups should be of dimension nbins*nfolds")
		}
	}
	# once we have groups, check whether they include enough p-values
	if(nbins > 1 & any(m_groups < 2)){
	  stop("Bins of size < 2 are currently not allowed for hypothesis weighting. Please lower the number of bins or try a different stratification method.")
	} else if (nbins > 1 & any(m_groups < 1000)){
		message("We recommend that you supply (many) more than 1000 p-values for meaningful data-driven hypothesis weighting results.")
	} 

	# start by making sure seed is correctly specified
	if (!is.null(seed)){
		#http://stackoverflow.com/questions/14324096/setting-seed-locally-not-globally-in-r?rq=1
		tmp <- runif(1)
		old <- .Random.seed
  		on.exit( { .Random.seed <<- old } )
  		set.seed(as.integer(seed)) #seed
	}

	if (nbins == 1){
		nfolds <- 1L
		nfolds_internal <- 1L
		nsplits_internal <- 1L
		message("Only 1 bin; IHW reduces to Benjamini Hochberg (uniform weights)")
		sorted_adj_p <- p.adjust(sorted_pvalues, method=adjustment_type, n=sum(m_groups))
		rjs <- sum(sorted_adj_p <= alpha)
		res <- list(lambda=0, fold_lambdas=0, rjs=rjs, sorted_pvalues=sorted_pvalues,
					sorted_weighted_pvalues = sorted_pvalues,
					sorted_adj_p=sorted_adj_p, sorted_weights=rep(1, length(sorted_pvalues)),
					sorted_groups=as.factor(rep(1, length(sorted_pvalues))),
					sorted_folds=as.factor(rep(1, length(sorted_pvalues))),
					weight_matrix = matrix(1))
	} else if (nbins > 1) {
		res <- ihw_internal(sorted_groups, sorted_pvalues, alpha, lambdas,
						m_groups,
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
						...)
	}

	if (return_internal){
		return(res)
	}

	# resort back to original form
	weights <- fill_nas_reorder(res$sorted_weights, nna, reorder_pvalues)
	pvalues <- fill_nas_reorder(res$sorted_pvalues, nna, reorder_pvalues)
	weighted_pvalues <- fill_nas_reorder(res$sorted_weighted_pvalues, nna, reorder_pvalues)
	adj_pvalues <- fill_nas_reorder(res$sorted_adj_p, nna, reorder_pvalues)
	groups     <- factor(fill_nas_reorder(res$sorted_groups, nna, reorder_pvalues),levels=group_levels)
	folds      <- factor(fill_nas_reorder(res$sorted_folds, nna, reorder_pvalues),levels=1:nfolds)
	
	if (is.matrix(covariates)) {
	  covariates <- fill_nas_reorder_dataframe(covariates, nna, seq_len(nrow(covariates)))
	} else {
	  covariates <- fill_nas_reorder(covariates, nna, seq_along(covariates))
	}


	df <- data.frame(pvalue = pvalues,
					 adj_pvalue = adj_pvalues,
				     weight = weights,
				     weighted_pvalue = weighted_pvalues,
					 group= groups,
					 covariate = covariates,
					 fold=as.factor(folds))

	ihw_obj <- new("ihwResult",
		 			df = df,
		 			weights = res$weight_matrix,
		 			alpha = alpha,
					nbins = as.integer(nbins),
					nfolds = nfolds,
					regularization_term = as.matrix(res$fold_lambdas),
					m_groups = as.integer(m_groups),
					penalty = penalty,
					covariate_type = covariate_type,
					adjustment_type = adjustment_type,
					reg_path_information = data.frame(),
		 			solver_information = list(),
					stratification_method = stratification_method,
					weight_matrices_forest = list(),
					n_censor_thres = integer(0),
					ntrees = integer(0))

	ihw_obj

}

# operate on pvalues which have already been sorted and without any NAs
#' @importFrom stats na.exclude p.adjust runif
ihw_internal <- function(sorted_groups, sorted_pvalues, alpha, lambdas,
								m_groups,
							    penalty="total variation",
							    seed=NULL,
								quiet=TRUE,
								sorted_folds=NULL,
								nfolds = 10L,
								nfolds_internal = nfolds,
								nsplits_internal = 1L,
								distrib_estimator = "distrib_estimator",
								lp_solver="lpsymphony",
								adjustment_type="BH",
								null_proportion = FALSE,
								null_proportion_level = 0.5,
								debug_flag=FALSE,
								folds_subset = NULL,
								...){

	folds_prespecified <- !is.null(sorted_folds)

	m <- length(sorted_pvalues)
	split_sorted_pvalues <- split(sorted_pvalues, sorted_groups)
	m_groups_available <- sapply(split_sorted_pvalues, length)


	#  do the k-fold strategy
	if (!is.null(seed)){
		#http://stackoverflow.com/questions/14324096/setting-seed-locally-not-globally-in-r?rq=1
		old <- .Random.seed
  		on.exit( { .Random.seed <<- old } )
  		set.seed(as.integer(seed)) #seed
	}

	if (!folds_prespecified){
		sorted_folds <- sample(1:nfolds, m, replace = TRUE)
	}

	sorted_weights <- rep(NA, m)
	fold_lambdas <- rep(NA, nfolds)
	weight_matrix <- matrix(NA, nlevels(sorted_groups), nfolds)

	if (debug_flag==TRUE){
		rjs_path_mat <- matrix(NA, nrow=nfolds, ncol=length(lambdas))
	}
	
	if (is.null(folds_subset)) folds_subset <- seq_len(nfolds)
	for (i in folds_subset) {

		if (!quiet) message(paste("Estimating weights for fold", i))
		# don't worry about inefficencies right now, they are only O(n) and sorting dominates this
		if (nfolds == 1){
			filtered_sorted_groups <- sorted_groups
			filtered_sorted_pvalues <- sorted_pvalues
			filtered_split_sorted_pvalues <- split(filtered_sorted_pvalues, filtered_sorted_groups)

			m_groups_holdout_fold <- m_groups
			m_groups_other_folds <- m_groups
		} else {
			filtered_sorted_groups <- sorted_groups[sorted_folds!=i]
			filtered_sorted_pvalues <- sorted_pvalues[sorted_folds!=i]
			filtered_split_sorted_pvalues <- split(filtered_sorted_pvalues, filtered_sorted_groups)

      if (!folds_prespecified){
				# TODO: add check to see if for common use case the first term is 0
				m_groups_holdout_fold <- (m_groups-m_groups_available)/nfolds +
						 m_groups_available - sapply(filtered_split_sorted_pvalues, length)
				m_groups_other_folds <- m_groups - m_groups_holdout_fold
			} else {
				m_groups_holdout_fold <- m_groups[,i]
				m_groups_other_folds <- rowSums(m_groups[,-i,drop=FALSE])
			}
		}
		# within each fold do iterations to also find lambda

		# TODO replace for loop by single call to ihw_convex with warm starts
		if (length(lambdas) > 1){
			rjs  <- matrix(NA, nrow=length(lambdas), ncol=nsplits_internal)
			for (k in 1:length(lambdas)){
				lambda <- lambdas[k]
				for (l in 1:nsplits_internal){

					rjs[k,l] <- ihw_internal(filtered_sorted_groups, filtered_sorted_pvalues, alpha,
											lambda, m_groups_other_folds,
									    	seed=NULL, quiet=quiet,
									    	nfolds=nfolds_internal,
									    	distrib_estimator = distrib_estimator,
									    	lp_solver=lp_solver,
									    	adjustment_type=adjustment_type,...)$rjs
				}
				if (debug_flag==TRUE){
					rjs_path_mat[i,k] <- rjs[k] #only works for nsplits_internal currently
				}
			}
			lambda <- lambdas[which.max(apply(rjs,1, mean))]
		} else if (length(lambdas) == 1){
			lambda <- lambdas
		} else {
			stop("something went wrong, need at least 1 value for regularization parameter!")
		}
		fold_lambdas[i] <- lambda

		# ok we have finally picked lambda and can proceed

		if (distrib_estimator=="grenander"){
			res <- ihw_convex(filtered_split_sorted_pvalues, alpha,
						   m_groups_holdout_fold, m_groups_other_folds,
						   penalty=penalty, lambda=lambda, lp_solver=lp_solver, 
						   adjustment_type=adjustment_type, quiet=quiet,...)
		} else if (distrib_estimator == "ECDF"){
			res <- ihw_milp(filtered_split_sorted_pvalues, alpha, m_groups_holdout_fold,
				           penalty=penalty, lambda=lambda, lp_solver=lp_solver, ...)
		} else {
			stop("This type of distribution estimator is not available.")
		}
    
		sorted_weights[sorted_folds == i] <- res$ws[sorted_groups[sorted_folds==i]]
		weight_matrix[,i] <- res$ws

		if (null_proportion){
			pi0_est <- weighted_storey_pi0(sorted_pvalues[sorted_folds==i],
				                  sorted_weights[sorted_folds == i],
				                  tau= null_proportion_level,
				                  m=sum(m_groups_holdout_fold))
			
			sorted_weights[sorted_folds == i] <- sorted_weights[sorted_folds == i]/pi0_est
			weight_matrix[,i] <- weight_matrix[,i]/pi0_est

		}

	}

	sorted_weighted_pvalues <- mydiv(sorted_pvalues, sorted_weights)


	sorted_adj_p <- p.adjust(sorted_weighted_pvalues, method = adjustment_type, n = sum(m_groups))
	rjs   <- sum(sorted_adj_p <= alpha)
	lst <- list(lambda=lambda, fold_lambdas=fold_lambdas, rjs=rjs, sorted_pvalues=sorted_pvalues,
					sorted_weighted_pvalues = sorted_weighted_pvalues,
					sorted_adj_p=sorted_adj_p, sorted_weights=sorted_weights,
					sorted_groups=sorted_groups, sorted_folds=sorted_folds,
					weight_matrix = weight_matrix)

	if (debug_flag == TRUE) {
		lst$rjs_path_mat <- rjs_path_mat
	}
	lst
}


#' @rdname ihw.default
#' @param formula \code{\link{formula}}, specified in the form pvalue~covariate1+covariate2+...
#' @param data data.frame from which the variables in formula should be taken
#' @export
ihw.formula <- function(formula, data=parent.frame(), ...){
	if (length(formula) != 3){
		stop("expecting formula of the form pvalue~covariate")
	}
	pv_name <- formula[[2]]
	
	cov_names <- labels(stats::terms(formula))
	cov_names <- lapply(cov_names, function(cov_name) parse(text = cov_name))
	nvar <- length(cov_names)
	
	pvalues <- eval(pv_name, data, parent.frame())
	
	if (nvar == 1) {
	  covariates <- eval(cov_names[[1]], data, parent.frame())
	} else {
	  m <- length(pvalues)
	  covariates <- vapply(cov_names, function(cov_name) eval(cov_name, data, parent.frame()), numeric(m))
	}
	
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
ihw.DESeqResults <- function(deseq_res, filter="baseMean",
								alpha=0.1, adjustment_type="BH",...){

    if (missing(filter)) {
    	filter <- deseq_res$baseMean
  	} else if (is.character(filter)){
  		stopifnot(filter %in% colnames(deseq_res))
  		filter <- deseq_res[[filter]]
  	}

  	stopifnot(length(filter) == nrow(deseq_res))

	ihw_res <- ihw(deseq_res$pvalue, filter, alpha=alpha,
							adjustment_type=adjustment_type,...)

	deseq_res$padj <- adj_pvalues(ihw_res)
 	deseq_res$weight <- weights(ihw_res)

	mcols(deseq_res)$type[names(deseq_res)=="padj"] <- "results"
	mcols(deseq_res)$description[names(deseq_res)=="padj"] <- paste("Weighted", adjustment_type,"adjusted p-values")

	mcols(deseq_res)$type[names(deseq_res)=="weight"] <- "results"
	mcols(deseq_res)$description[names(deseq_res)=="weight"] <- "IHW weights"

	metadata(deseq_res)[["alpha"]] <- alpha
	metadata(deseq_res)[["ihwResult"]] <- ihw_res
 	deseq_res
}

#' @importFrom slam simple_triplet_zero_matrix simple_triplet_matrix
#' @importFrom lpsymphony lpsymphony_solve_LP
ihw_convex <- function(split_sorted_pvalues, alpha, m_groups, m_groups_grenander,
		penalty="total variation", lambda=Inf, lp_solver="gurobi",
		adjustment_type= "BH", grenander_binsize = 1, quiet=quiet){

	# m_groups used for balancing (weight budget)
	# m_groups_grenander (used for grenander estimator)
	# Practically always it will hold: m_groups \approx m_groups_grenander

	nbins <- length(split_sorted_pvalues)

	if (lambda==0){
		ws <- rep(1, nbins)
		return(list(ws=ws))
	}

  
	# preprocessing:  Set very low p-values to 0, otherwise LP solvers have problems due to numerical instability
	# Note: This only affects the internal optimization, the higher level functions will
	# still return adjusted pvalues based on the original p-values

	split_sorted_pvalues <- lapply(split_sorted_pvalues, function(x) ifelse(x > 10^(-20), x, 0))
	m <- sum(m_groups)

	if (nbins != length(m_groups)){
		stop("length of m_groups should be equal to number of bins")
	}

	#lapply grenander...
	if (!quiet) message("Applying Grenander estimator within each bin.")
	grenander_list <- mapply(presorted_grenander,
		                    split_sorted_pvalues,
		                    m_groups_grenander,
		                    grenander_binsize = grenander_binsize,
		                    quiet=quiet,
		                    SIMPLIFY=FALSE)


	#set up LP
	nconstraints_per_bin <- sapply(grenander_list, function(x) x$length)
	nconstraints <- sum(nconstraints_per_bin)
	i_yi <- 1:nconstraints
	j_yi <- rep(1:nbins,  times=nconstraints_per_bin)
	v_yi <- rep(1,  nconstraints)

	i_ti <- 1:nconstraints
	j_ti <- nbins + rep(1:nbins, times=nconstraints_per_bin)
	v_ti <- unlist(lapply(grenander_list, function(x) -x$slope.knots))

	constr_matrix <- slam::simple_triplet_matrix(c(i_yi, i_ti), c(j_yi,j_ti), c(v_yi, v_ti))
	rhs <- unlist(lapply(grenander_list, function(x) x$y.knots-x$slope.knots*x$x.knots))

	obj <- c(m_groups/m*nbins*rep(1,nbins), rep(0,nbins))

	if (lambda < Inf){
		if (penalty == "total variation"){
			# -f + t_g - t_{g-1} <= 0
			i_fi <- rep(1:(nbins-1),3)
			j_fi <-	c((nbins+2):(2*nbins), (nbins+1):(2*nbins-1), (2*nbins+1):(3*nbins-1))
			v_fi <- c(rep(1,nbins-1), rep(-1,nbins-1), rep(-1, nbins-1))
			absolute_val_constr_matrix_1 <- slam::simple_triplet_matrix(i_fi,j_fi,v_fi)

			# -f - t_g + t_{g-1} <= 0
			i_fi <- rep(1:(nbins-1),3)
			j_fi <-	c((nbins+2):(2*nbins), (nbins+1):(2*nbins-1), (2*nbins+1):(3*nbins-1))
			v_fi <- c(rep(-1,nbins-1), rep(1,nbins-1), rep(-1, nbins-1))
			absolute_val_constr_matrix_2 <- slam::simple_triplet_matrix(i_fi,j_fi,v_fi)

			constr_matrix <- rbind(cbind(constr_matrix, slam::simple_triplet_zero_matrix(nconstraints, nbins-1,mode="double")),
							    absolute_val_constr_matrix_1,
								absolute_val_constr_matrix_2)

			obj <- c(obj, rep(0, nbins-1))

			total_variation_constr <- matrix(c(rep(0,nbins), -lambda*m_groups/m, rep(1,nbins-1)),nrow=1)
			constr_matrix <- rbind(constr_matrix, total_variation_constr)

			rhs <- c(rhs,rep(0, 2*(nbins-1)),0) #add RHS for absolute differences and for total variation penalty

		} else if (penalty == "uniform deviation"){
			# -f + m t_g - sum_g m_i t_i <= 0

			absolute_val_constr_matrix_1 <- matrix(rep(-m_groups,nbins), nbins,nbins, byrow=TRUE)
			diag(absolute_val_constr_matrix_1) <- -m_groups + m
			absolute_val_constr_matrix_1 <- cbind( slam::simple_triplet_zero_matrix(nbins, nbins,mode="double"),
												   absolute_val_constr_matrix_1,
										           -diag(nbins))

			# -f - m t_g +  sum_g m_i t_i  <= 0

			absolute_val_constr_matrix_2 <- matrix(rep(+m_groups,nbins), nbins,nbins, byrow=TRUE)
			diag(absolute_val_constr_matrix_2) <- m_groups - m
			absolute_val_constr_matrix_2 <- cbind( slam::simple_triplet_zero_matrix(nbins, nbins,mode="double"),
												   absolute_val_constr_matrix_2,
										           -diag(nbins))

			constr_matrix <- rbind(cbind(constr_matrix, slam::simple_triplet_zero_matrix(nconstraints, nbins,mode="double")),
							    absolute_val_constr_matrix_1,
								absolute_val_constr_matrix_2)

			obj <- c(obj, rep(0, nbins))

			total_variation_constr <- matrix(c(rep(0,nbins), -lambda*m_groups, rep(1,nbins)),nrow=1)
			constr_matrix <- rbind(constr_matrix, total_variation_constr)

			rhs <- c(rhs,rep(0, 2*nbins),0) #add RHS for absolute differences and for total variation penalty

		}
	}

	if (adjustment_type == "BH"){
		# incorporate the FDR constraint
		fdr_constr<- matrix(c(rep(-alpha,nbins)*m_groups,
						      rep(1,nbins)*m_groups,
						      rep(0,ncol(constr_matrix)-2*nbins)), nrow=1)
		constr_matrix <- rbind(constr_matrix, fdr_constr)
		rhs <- c(rhs,0)
	} else if (adjustment_type == "bonferroni"){
		# incorporate the FWER constraint
		fwer_constr<- matrix(c(rep(0,nbins),
						      rep(1,nbins)*m_groups,
						      rep(0,ncol(constr_matrix)-2*nbins)), nrow=1)
		constr_matrix <- rbind(constr_matrix, fwer_constr)
		rhs <- c(rhs,alpha)
	} else {
		stop("Only BH-like and bonferroni-like adjustments available.")
	}
	nvars <- ncol(constr_matrix)

	if (!quiet) message("Starting to solve LP.")

	if (lp_solver == "gurobi"){

		#if (!require("gurobi", quietly=TRUE)){
		#	stop("Gurobi solver appears not be installed. Please use lpsymphony or install gurobi.")
		#}

		if (!requireNamespace("Matrix", quietly=TRUE)){
			stop("Matrix package required to use gurobi in IHW.")
		}

		model <- list()
		model$A <- Matrix::sparseMatrix(i=constr_matrix$i, j=constr_matrix$j, x=constr_matrix$v,
           		dims=c(constr_matrix$nrow, constr_matrix$ncol))
		model$obj <- obj
		model$modelsense <- "max"
		model$rhs        <- rhs
		model$lb         <- 0
		model$ub         <- 2
		model$sense      <- '<'

		params <- list(OutputFlag=1)
		res <- gurobi(model, params)
		sol <- res$x
		solver_status <- res$status

	} else if (lp_solver=="lpsymphony") {

		#rsymphony_bounds <- list(lower=list(ind=1:nvars, val=rep(0,nvars)),
		#					 upper=list(ind=1:nvars, val=rep(1,nvars)))

		#if (is.infinite(time_limit)) time_limit <- -1
		#if (is.infinite(node_limit)) node_limit <- -1
		res <- lpsymphony::lpsymphony_solve_LP(obj, constr_matrix, rep("<=", nrow(constr_matrix)),
			rhs, #bounds= rsymphony_bounds,
			max = TRUE, verbosity = -2, first_feasible = FALSE)
		sol <- res$solution
		solver_status <- res$status
	} else {
		stop("Only gurobi and lpsymphony solvers currently supported.")
	}
  
	# catch negative thresholds due to numerical rounding
	ts <- pmax(sol[(1:nbins)+nbins],0)
	ws <- thresholds_to_weights(ts, m_groups)	#\approx ts/sum(ts)*nbins

	# return matrix as follows: nbins x nlambdas)
	return(list(ws=ws))
}


#' @importFrom fdrtool gcmlcm
presorted_grenander <- function(sorted_pvalues, m_total=length(sorted_pvalues),
							grenander_binsize = 1,
							quiet=TRUE){

  	unique_pvalues <- unique(sorted_pvalues)
  	ecdf_values <- cumsum(tabulate(match(sorted_pvalues, unique_pvalues)))/m_total

  	if (min(unique_pvalues) > 0){
  		# I think fdrtool neglects this borderline case and this causes returned object
  		# to be discontinuous hence also not concave
  		unique_pvalues <- c(0,unique_pvalues)
  		ecdf_values   <- c(0, ecdf_values)
  	}

  	if (max(unique_pvalues) < 1){
  		unique_pvalues <- c(unique_pvalues,1)
  		ecdf_values    <- c(ecdf_values,1)
  	}

  	if (grenander_binsize != 1){
  		nmax = length(unique_pvalues)
  		idx_thin <- c(seq(1, nmax-1, by=grenander_binsize), nmax)
  		unique_pvalues <- unique_pvalues[idx_thin]
  		ecdf_values <- ecdf_values[idx_thin]
  	}

  	ll <- fdrtool::gcmlcm(unique_pvalues, ecdf_values, type="lcm")
	ll$length <- length(ll$slope.knots)
	ll$x.knots <- ll$x.knots[-ll$length]
  	ll$y.knots <- ll$y.knots[-ll$length]
	if (!quiet) message(paste("Grenander fit with", ll$length, "knots."))
	ll
}


# mostly deprecated function (users should not be using it)
# needed for reproducibility of manuscript
# to demonstrate that "naive IHW" algorithm does not control type-I error

ihw_milp <- function(split_sorted_pvalues, alpha, m_groups, lambda=Inf, lp_solver="gurobi", 
		quiet=quiet,
		optim_pval_threshold=1,
		penalty="total variation",
		lp_relaxation=FALSE,
		t_bh=0,
		bh_rejections =0,
		rejections_bound=NULL,
		time_limit=Inf,
		node_limit=Inf,
		threads=0,
		solution_limit=Inf,
		mip_gap_abs=Inf,
		mip_gap = 10^(-4),
		MIPFocus=0){



	nbins <- length(split_sorted_pvalues)

	if (lambda==0){
		ws <- rep(1, nbins)
		return(list(ws=ws))
	}

	m <- sum(m_groups)
	mtests <- m

	if (nbins != length(m_groups)){
		stop("length of m_groups should be equal to number of bins")
	}


	if (optim_pval_threshold < 1){
		# keep at least 2 p-values in each group, so we do not have to add special cases for n=0 or n=1
		split_sorted_pvalues <- lapply(split_sorted_pvalues,
			function(p) p[p <= max(p[2],optim_pval_threshold)])
	}

	pmax <- max(unlist(split_sorted_pvalues))
	split_sorted_pvalues <- lapply(split_sorted_pvalues, function(p) p/pmax)
	#ts_bh <- sapply(split_sorted_pvalues, function(ps) max(0, ps[ps <= t_bh]))


	# create a simple_triplet_matrix (sparse matrix) representation,

	# function to create sparse-representation of group-wise blocks enforcing the constraints
	# z_1 >= z_2 >= .. in each group  i.e. [(m_group -1) X m_group] matrix of the form:
	#
	# 1 -1  0  ......
	# 0  1 -1  ......
	# 0  0  1  ......
 	# ...............
 	#
 	# Also parameters to make it easy to diagonally combine these blocks afterwards

	single_block_triple_representation <- function(start, p_length,prev_block){
		one_block <- cbind(start:(start+p_length-2), (prev_block+start):(prev_block+start+p_length-2), rep(1,p_length-1))
		minus_one_block <- cbind(start:(start+p_length-2), (prev_block+1+start):(prev_block+start+p_length-1), rep(-1, p_length-1))
	    rbind(one_block, minus_one_block)
	}



	# create, combine above blocks and store them as simple triplet matrix
	p_lengths <- sapply(split_sorted_pvalues,length)
	prev_blocks <- 0:(nbins-1)
	starts <- cumsum(c(1,p_lengths[-nbins]-1))

	full_triplet <- do.call(rbind, mapply(single_block_triple_representation, starts,p_lengths, prev_blocks, SIMPLIFY=FALSE))
	full_triplet <- slam::simple_triplet_matrix(full_triplet[,1],full_triplet[,2], full_triplet[,3])


	# add final constraint to matrix which corresponds to control of plug-in FDR estimate
	diff_coefs <- unlist(mapply(function(p, m_group) diff(c(0,p))* m_group , split_sorted_pvalues, m_groups))
	#plugin_fdr_constraint <- (alpha - diff_coefs * mtests/m)/m #alpha - diff_coefs * mtests/m
	plugin_fdr_constraint <- (alpha/pmax - diff_coefs * mtests/m)/m #alpha - diff_coefs * mtests/m

	full_triplet <- rbind(full_triplet, matrix(plugin_fdr_constraint, nrow=1))

	z_vars <- ncol(full_triplet) #number of binary variables

	obj        <- rep(1, z_vars) # maximize number of rejections

	if (is.finite(lambda)){
		# start with new code which introduces regularization

		# first we introduce #nbins new random variables T_g, such that
		# t_{g} <= T_g < t'_{g} [in practice the second < has to be turned into a <=]
		# this reflects the fact that for all these values of T_g, still the same hypotheses will be rejected
		# in stratum g.

		full_triplet <- cbind(full_triplet, matrix(0, nrow(full_triplet), nbins))

		# Below we modify objective a tiny bit by enforcing some minimization on \sum m_g T_g
		# because 1/(4m)* \sum m_g T_g << 1 we ensure that our objective of interest (number of rejections) does not change
		# Motivation for this soft constraint:
		# 1) Make final solution (in terms of threshold T_g) "somewhat" unique
    	# 2) (more importantly): assist solver in preventing possible violation of plugin FDR control due to numerical issues
    	# due to small p-values.

		obj <- c(obj, -.25*m_groups/m)


		# first g constraints representing t_{g} <= T_g can be represented by [-A I_g] with A:
		#
		#   1 1 |          |          |                  # t_1
		#       |  1  1  1 |          |					 # t_2
		#       |          |  1  1  1 |                  # t_3
		#
		#	 (Replace 1s by diff(c(0,p)) values)

		pdiff_list <- lapply(split_sorted_pvalues, function(p) diff(c(0,p)))

		left_ineq_matrix_block <- function(g, p_length, pdiff){
			mat <- 	matrix(0, nbins, p_length)
			mat[g,] <- - pdiff
			mat
 		}


 		# now put above blocks together
 		ineq_matrix <-  rbind( cbind(do.call(cbind, mapply(left_ineq_matrix_block, 1:nbins, p_lengths, pdiff_list, SIMPLIFY=FALSE)),
 							  	diag(1, nbins)))

 		# for these T_g's we also require the plug-in FDR estimate to be controlled
 		# (in fact we could drop the old constraint on the t_g's, but the extra constraint could help with the relaxations)

 		# we also divide the constraint by m (arbitrary) so that the coefficient range of the matrix is not too large
 		# which could cause numerical problems
		#plugin_fdr_constraint2 <- matrix(c(rep(alpha/m,z_vars), -m_groups/m), nrow=1)
		plugin_fdr_constraint2 <- matrix(c(rep(alpha/m/pmax,z_vars), -m_groups/m), nrow=1)

		# put these constraints together..
		full_triplet <- rbind(full_triplet, ineq_matrix, plugin_fdr_constraint2)


		if (penalty=="total variation"){
			# this implements a penalty of the form: \sum_{g=2}^G |w_g - w_{g-1}| <= reg_par

			# |T_g - T_{g-1}| = f_g  introduce absolute value constants, i.e. introduce nbins-1 new continuous variables

			full_triplet <- cbind(full_triplet, matrix(0,nrow(full_triplet),nbins-1))
			obj <- c(obj, rep(0, nbins-1))

			# we now need inequalities: T_g - T_{g-1} + f_g >= 0    and -T_g + T_{g-1} + f_g >= 0
			# start with matrix diff_matrix (nbins-1) X nbins as a scaffold i.e.
			#
			#  -1  1
			#     -1  1
			#        -1  1


			diff_matrix <- diag(-1,nbins-1,nbins) + cbind(rep(0,nbins-1), diag(1,nbins-1))
			abs_constraint_matrix <- rbind(cbind(slam::simple_triplet_zero_matrix(nbins-1, z_vars,mode="double"),
											diff_matrix, slam::simple_triplet_diag_matrix(1, nrow=nbins-1)),
									   cbind(slam::simple_triplet_zero_matrix(nbins-1, z_vars,mode="double"),
											-diff_matrix, slam::simple_triplet_diag_matrix(1, nrow=nbins-1)))


 			# add final regularization inequality:
 			# 			\sum |w_g - w_{g-1}| <= reg_par
 			# <=>       \sum |T_g - T_{g-1}|*m <= reg_par * \sum m_g T_g
 			regularization_constraint <- matrix(c(rep(0,z_vars), lambda/m*m_groups, rep(-1, nbins-1)), nrow=1)

 			# to be used for ub afterwards
 			regularization_ub <- c(rep(1, nbins), rep(2, nbins-1))

 		} else if (penalty=="uniform deviation"){
 			# this implements a penalty of the form: \sum_{g=1}^G |w_g -1| <= reg_par

 			# |m*T_g - \sum_{i=1}^G m_i T_i| = f_g  introduce absolute value constants, i.e. introduce nbins new continuous variables

			full_triplet <- cbind(full_triplet, matrix(0,nrow(full_triplet),nbins))
			obj <- c(obj, rep(0, nbins))

			# we now need inequalities: m*T_g - \sum_{i=1}^G m_i T_i + f_g >= 0   and -m*T_g  +  \sum_{i=1}^G m_iT_i + f_g >= 0
			# start with matrix diff_matrix (nbins X nbins) as a scaffold i.e.
			#
			#  (m-m_1)    -m_2      -m_3   ......
			#    -m_1  (m - m_2)    -m_3   ......
			#    -m_1     -m_2    (m-m_3)  ......
			#     .         .         .    .
			#     .         .         .      .
			#     .         .         .        .



			diff_matrix <- diag(m, nbins, nbins) - matrix(rep(m_groups, each=nbins), nrow=nbins) 

			abs_constraint_matrix <- rbind(cbind(slam::simple_triplet_zero_matrix(nbins, z_vars,mode="double"),
											diff_matrix, slam::simple_triplet_diag_matrix(1, nrow=nbins)),
									   cbind(slam::simple_triplet_zero_matrix(nbins, z_vars,mode="double"),
											-diff_matrix, slam::simple_triplet_diag_matrix(1, nrow=nbins)))


 			# add final regularization inequality:
 			# 			\sum |w_g - 1| <= reg_par
 			# <=>       \sum |m*T_g - \sum m_i T_i| <= reg_par * \sum m_i T_i
 			regularization_constraint <- matrix(c(rep(0,z_vars), lambda*m_groups, rep(-1, nbins)), nrow=1)

 			# to be used for ub afterwards
 			regularization_ub <- c(rep(1, nbins), rep(2*m, nbins))

 		} else {
 			stop("No such regularization penalty currently supported.")
 		}
 		# put all constraints together
 		full_triplet <- rbind(full_triplet, abs_constraint_matrix, regularization_constraint)

 	}

    # finally add an inequality ensuring that we get at least as many rejections as BH
    # added this feature as an ad-hoc solution to a previous problematic regularization implementation
    # but actually seems to also speed up the optimization problem!

    bh_rj_inequality <- matrix(c(rep(1, z_vars), rep(0, ncol(full_triplet)-z_vars)), nrow=1)
    full_triplet <- rbind(full_triplet, bh_rj_inequality)

    if (!is.null(rejections_bound)){
    	full_triplet <- rbind(full_triplet, -bh_rj_inequality)
    }

    nrows <- nrow(full_triplet)
    model_rhs <- c(rep(0,nrows-1), bh_rejections)

    if (!is.null(rejections_bound)){
    	model_rhs <- c(rep(0,nrows-2), -rejections_bound)
    }

    if (lp_relaxation){
    	model_vtype <- rep('C', ncol(full_triplet))
    } else {
    	model_vtype <- c(rep('B', z_vars), rep('C', ncol(full_triplet)-z_vars))
    }
    # model ub will actually depend on which penalty we choose
    # model ub = 1 for binary variables, 1 for thresholds, x for f_g, where
    # x=2 for total variation, x=2*m for uniform deviation
    model_ub <- c(rep(1, z_vars))
    if (z_vars < nrows){
    	model_ub <- c(model_ub, regularization_ub)
    }

	if (lp_solver=="lpsymphony") {

		rsymphony_ub <- list(upper = list(ind = 1:ncol(full_triplet), val = model_ub))

		if (is.infinite(time_limit)) time_limit <- -1
		if (is.infinite(node_limit)) node_limit <- -1
		res<- lpsymphony::lpsymphony_solve_LP(obj, full_triplet, rep(">=", nrows), 
			model_rhs,
    		types = model_vtype, bounds= rsymphony_ub,
			max = TRUE, verbosity = -2, time_limit = time_limit,
			node_limit = node_limit, gap_limit = ifelse(mip_gap <= 10^(-4), -1, mip_gap*100), #now default mip_gap=10^(-4) as in Gurobi while -1 default for Rsymphony
			first_feasible = FALSE)
		sol <- res$solution
		print(str(res))
		solver_status <- res$status

	} else if (lp_solver=="gurobi"){


		#if (!require("gurobi", quietly=TRUE)){
		#	stop("Gurobi solver appears not be installed. Please use lpsymphony or install gurobi.")
		#}

		if (!requireNamespace("Matrix", quietly=TRUE)){
			stop("Matrix package required to use gurobi in IHW.")
		}
		# keep code for commercial solver for now
		#speed up compared to Symphony appears to be at least ~2-10, depending on problem

		# convert simple_triplet_matrix of pkg Slam to sparse matrix of pkg Matrix
		# gurobi has problems with simple_triplet_matrix for an (undocumented) reason (?bug)
		full_triplet <-  Matrix::sparseMatrix(i=full_triplet$i, j=full_triplet$j, x=full_triplet$v,
           dims=c(full_triplet$nrow, full_triplet$ncol))

		model <- list()
		model$A <- full_triplet
		model$obj <- obj
		model$modelsense <- "max"
		model$rhs        <- model_rhs
		model$lb         <- 0
		model$ub         <- model_ub
		model$sense      <- '>'
		model$vtype      <- model_vtype

		# the parameters below make the optimization procedure MUCH slower
		# but they are actually necessary... (need to check what happens with open source solvers)
		# turns out that enforcing T_g >= t_g is numerically very hard in the context of this problem

		params <- list(NumericFocus=3, FeasibilityTol=10^(-9), ScaleFlag=0 , Quad=1, IntFeasTol=10^(-9),
							OutputFlag = 0, #Presolve=0,#10,
							Threads=threads,MIPFocus=MIPFocus,
							MIPGap = mip_gap)

		if (is.finite(time_limit)) params$TimeLimit <- time_limit
		if (is.finite(node_limit)) params$NodeLimit <- node_limit
		if (is.finite(solution_limit)) params$SolutionLimit <- solution_limit
		if (is.finite(mip_gap_abs)) params$MIPGapAbs <- mip_gap_abs

		res <- gurobi(model, params)
		sol <- res$x
		solver_status <- res$status

	} else {
		stop("Only gurobi and lpsymphony solvers are currently supported.")
	}

	# next try to extract the thresholds/weights from solver output
	lengths <- cumsum(sapply(split_sorted_pvalues, length))
	start_lengths <- c(1, lengths[-nbins]+1)

	pidx_list <- mapply(function(start,end) sol[start:end], start_lengths,lengths, SIMPLIFY=FALSE)

	get_threshold <- function(plist, pidx){
		max_idx <- which.max(which(pidx==1))
		ifelse(length(max_idx)==0, .0, plist[max_idx])
	}

	#get_threshold2 <- function(pdiff, pidx){
	#	sum(pdiff*pidx)
	#}

	t_thresholds1 <- mapply(get_threshold, split_sorted_pvalues, pidx_list)     # to be stored in ddhw object
	#t_thresholds2 <- mapply(get_threshold2, pdiff_list, pidx_list)
	t_thresholds <- if (is.infinite(lambda)){
					 	t_thresholds1
					} else {
						sol[(z_vars+1):(z_vars+nbins)]  # to be used for calculating weights
					}

    ####################################################################################################
    #	Possible checks of whether solver did not do (important) numerical errors: #TODO
    #
    # 1) t_thresholds >= t_thresholds1
    # 2) Plugin FDR has to be controlled (with respect to both t_thresholds, t_thresholds2)
    # 3) total variation of weights has to be less than lambda
    #
    #####################################################################################################

	ws <-  if (all(t_thresholds1 == .0) || all(t_thresholds == .0)){
				rep(1,nbins)
			} else {
				t_thresholds*m/sum(m_groups*t_thresholds) #might need adjustment if mtests > m
			}

	return(list(ws=ws))
}
