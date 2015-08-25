#' DDHW: Data-Driven Hypothesis Weights
#'
#' Given a vector of p-values, a vector of filter-statistics which are independent of the p-values under the null hypothesis and
#' a nominal significance level alpha, DDHW learns multiple testing weights and then applies the weighted Benjamini Hochberg 
#' procedure. When the filter-statistic is informative of the power of the individual tests, this procedure can increase
#' power.
#'
#' @param pvalues  Numeric vector of unadjusted p-values.
#' @param filter_statistics  Vector which contains the one-dimensional filter-statistics (covariates, independent under the H0 of the p-value) 
#'                for each test. Can be numeric or a factor. (If numeric it will be converted into factor by binning.)
#' @param alpha   Numeric, sets the nominal level for FDR control.
#' @param filter_statistic_type  "ordinal" or "nominal" (i.e. whether filter statistics can be sorted in increasing order or not)
#' @param nbins  Integer, number of groups into which p-values will be split based on filter_statistic. Use "auto" for
#'             automatic selection of the number of bins. Only applicable when filter_statistics is not a factor.
#' @param quiet  Boolean, if False a lot of messages are printed during the fitting stages.
#' @param nfolds Number of folds into which the p-values will be split for the pre-validation procedure
#' @param nfolds_internal  Within each fold, a second  (nested) layer of cross-validation can be conducted to choose a good
#'              regularization parameter. This parameter controls the number of nested folds.
#' @param lambdas  Numeric vector which defines the grid of possible regularization parameters.
#'				Use "auto" for automatic selection.
#' @param seed Integer or NULL. Split of hypotheses into folds is done randomly. To have output of the function be reproducible, 
#'	we set a seed. Use NULL if you don't want a seed.
#' @param lp_solver  Character ("lpsymphony" or "gurobi"). Internally, DDHW solves a sequence of linear programs, which
#'        can be solved with either of these solvers.
#' @param return_internal Returns a lower level representation of the output (only useful for debugging purposes).
#' @param ... Arguments passed to internal functions.
#'
#' @return A ddhwResult object.
#' @seealso ddhwResult, plot_ddhw
#' 
#' @examples
#'
#'    set.seed(1)
#'    X   <- runif(20000, min=0.5, max=4.5) #covariate
#'    H   <- rbinom(20000,1,0.1)            #hypothesis true or false
#'    Z   <- rnorm(20000, H*X)              #Z-score
#'    pvalue <- 1-pnorm(Z)                  #pvalue
#'    ddhw_res <- ddhw(pvalue, X, .1)
#'
#'
#' @export
ddhw <- function(pvalues, filter_statistics, alpha,
						filter_statistic_type = "ordinal",
						nbins = "auto",
						quiet =T ,
						nfolds = 5L,
						nfolds_internal = 5L,
						lambdas = "auto",
						seed = 1L,
						lp_solver="lpsymphony",
						return_internal=FALSE,
						...){

	# This function essentially wraps the lower level function ddhw_internal
	# e.g. takes care of NAs, sorts pvalues and then
	# returns a nice ddhw object

	nfolds <- as.integer(nfolds)
	nfolds_internal <- as.integer(nfolds_internal)

	if (nfolds==1){
		warning("Using only 1 fold! Only use this if you want to learn the weights, but NEVER for testing!")
	}


	# gracefully handle NAs
   	nna <- !is.na(pvalues)

   	if (all(nna)){
        nna <- TRUE
   	}

   	# filter our p-values
    pvalues <- pvalues[nna]
    filter_statistics <- filter_statistics[nna]
    weights <- rep(NA, length(pvalues))

    if (any(is.na(filter_statistics))){
    	stop("Filter statistics corresponding to non-NA p-values should never be NA. Aborting.")
    }

	if (nbins == "auto"){
		nbins <- min(300, floor(length(pvalues)/1500)) # rule of thumb..
	}

	if ((length(lambdas)== 1) & (lambdas[1] == "auto")){
		# just a few for now until I get warm starts of the LP solvers to work
		lambdas <- c(1, nbins/8, nbins/4, nbins/2, nbins, Inf)
	}

	if (filter_statistic_type =="ordinal" & is.numeric(filter_statistics)){
		groups <- groups_by_filter(filter_statistics, nbins)
	} else if (is.factor(filter_statistics)){
		groups <- filter_statistics
		nbins <- nlevels(groups)
		if (filter_statistic_type == "nominal"){
			lambdas <- Inf
		}
	} else {
		stop("filter_statistics are not of the appropriate type")
	}


	# once we have groups, check whether they include enough p-values
	if (any(table(groups) < 1000)){
		stop("Data-driven choice of weights requires at least 1000 p-values per stratum.")
	}

	group_levels <- levels(groups)
	# sort pvalues globally
	order_pvalues <- order(pvalues) #TODO BREAK ties by filter statistic rank
	reorder_pvalues <- order(order_pvalues)

	sorted_groups <- groups[order_pvalues]
	sorted_pvalues <- pvalues[order_pvalues]

	res <- ddhw_internal(sorted_groups, sorted_pvalues, alpha, lambdas,
						quiet=quiet,
						nfolds=nfolds,
						nfolds_internal=nfolds_internal,
						seed=seed,
						lp_solver=lp_solver,

						...)

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
	filter_statistics <- fill_nas_reorder(filter_statistics, nna, 1:length(filter_statistics))


	df <- data.frame(pvalue = pvalues,
					 padj = adj_pvalues,
				     weights = weights,
				     weighted_pvalues = weighted_pvalues,
					 group= groups,
					 filter_statistic = filter_statistics,
					 folds=as.factor(folds))

	ddhw_obj <- new("ddhwResult",
		 			df = df,
		 			weights = res$weight_matrix,
		 			alpha = alpha,
					nbins = as.integer(nbins),
					nfolds = nfolds,
					regularization_term = res$fold_lambdas,
					penalty = "total variation",
					filter_statistic_type = filter_statistic_type,
					reg_path_information = data.frame(),
		 			solver_information = list())


	ddhw_obj

}

# operate on pvalues which have already been sorted and without any NAs
ddhw_internal <- function(sorted_groups, sorted_pvalues, alpha, lambdas,
							    seed=NULL,
								quiet=TRUE,
								nfolds = 10,
								nfolds_internal = nfolds,
								lp_solver="lpsymphony",
								debug_flag=FALSE){

	# TODO: check if lambdas are sorted the way I want them to

	m <- length(sorted_pvalues)
	split_sorted_pvalues <- split(sorted_pvalues, sorted_groups)
	m_groups <- sapply(split_sorted_pvalues, length)


	#  do the k-fold strategy
	if (!is.null(seed)) set.seed(as.integer(seed)) #seed

	sorted_folds <- sample(1:nfolds, m, replace = TRUE)

	sorted_weights <- rep(NA, m)
	fold_lambdas <- rep(NA, nfolds)
	weight_matrix <- matrix(NA, nlevels(sorted_groups), nfolds)

	if (debug_flag==TRUE){
		rjs_path_mat <- matrix(NA, nrow=nfolds, ncol=length(lambdas))
	}
	for (i in 1:nfolds){

		if (!quiet) message(paste("Estimating weights for fold", i))
		# don't worry about inefficencies right now, they are only O(n) and sorting dominates this
		if (nfolds == 1){
			filtered_sorted_groups <- sorted_groups
			filtered_sorted_pvalues <- sorted_pvalues
			filtered_split_sorted_pvalues <- split(filtered_sorted_pvalues, filtered_sorted_groups)
		} else {
			filtered_sorted_groups <- sorted_groups[sorted_folds!=i]
			filtered_sorted_pvalues <- sorted_pvalues[sorted_folds!=i]
			filtered_split_sorted_pvalues <- split(filtered_sorted_pvalues, filtered_sorted_groups)
		}
		# within each fold do iterations to also find lambda

		# TODO replace for loop by single call to ddhw_convex with warm starts
		if (length(lambdas) > 1){
			rjs  <- rep(NA, length(lambdas))
			for (k in 1:length(lambdas)){
				lambda <- lambdas[k]
				# To consider: internal nfolds does not have to be the same, could be 2 for speedup
				rjs[k] <- ddhw_internal(filtered_sorted_groups, filtered_sorted_pvalues, alpha, lambda,
									    seed=seed, quiet=quiet,
									    nfolds=nfolds_internal, lp_solver=lp_solver)$rjs
				if (debug_flag==TRUE){
					rjs_path_mat[i,k] <- rjs[k]
				}
			}
			lambda <- lambdas[which.max(rjs)]
		} else if (length(lambdas) == 1){
			lambda <- lambdas
		} else {
			stop("something went wrong, need at least 1 value for regularization parameter!")
		}
		fold_lambdas[i] <- lambda

		# ok we have finally picked lambda and can proceed
		if (nfolds==1){
			m_groups_holdout_fold <- m_groups
		} else {
			m_groups_holdout_fold <- m_groups - sapply(filtered_split_sorted_pvalues, length)
		}

		res <- ddhw_convex(filtered_split_sorted_pvalues, alpha, m_groups_holdout_fold,
						   lambda=lambda, lp_solver=lp_solver, quiet=quiet)

		sorted_weights[sorted_folds == i] <- res$ws[sorted_groups[sorted_folds==i]]
		weight_matrix[,i] <- res$ws
	}

	sorted_weighted_pvalues <- mydiv(sorted_pvalues, sorted_weights)
	sorted_adj_p <- p.adjust(sorted_weighted_pvalues, method = "BH")
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

#' @importFrom slam simple_triplet_zero_matrix simple_triplet_matrix
#' @importFrom lpsymphony lpsymphony_solve_LP
ddhw_convex <- function(split_sorted_pvalues, alpha, m_groups, lambda=Inf, lp_solver="gurobi", quiet=quiet){

	# preprocessing:  Set very low p-values to 0, otherwise LP solvers have problems due to numerical instability
	# Note: This only affects the internal optimization, the higher level functions will
	# still return adjusted pvalues based on the original p-values
	split_sorted_pvalues <- lapply(split_sorted_pvalues, function(x) ifelse(x > 10^(-20), x, 0))
	nbins <- length(split_sorted_pvalues)
	m <- sum(m_groups)

	if (nbins != length(m_groups)){
		stop("length of m_groups should be equal to number of bins")
	}

	if (lambda==0){
		ws <- rep(1, m_groups)
		return(list(ws=ws))
	}

	#lapply grenander...
	if (!quiet) message("Applying Grenander estimator within each bin.")
	grenander_list <- lapply(split_sorted_pvalues, presorted_grenander, quiet=quiet)


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

	}

	# incorporate the FDR constraint
	fdr_constr<- matrix(c(rep(-alpha,nbins)*m_groups, rep(1,nbins)*m_groups, rep(0,ncol(constr_matrix)-2*nbins)), nrow=1)
	constr_matrix <- rbind(constr_matrix, fdr_constr)
	nvars <- ncol(constr_matrix)
	rhs <- c(rhs,0)

	if (!quiet) message("Starting to solve LP.")

	if (lp_solver == "gurobi"){
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
		res <- gurobi::gurobi(model, params)
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


	return(list(ws=ws))
}


#' @importFrom fdrtool gcmlcm
presorted_grenander <- function(sorted_pvalues, quiet=TRUE){
  	n  <- length(sorted_pvalues)
  	unique_pvalues <- unique(sorted_pvalues)
  	ecdf_values <- cumsum(tabulate(match(sorted_pvalues, unique_pvalues)))/n
  	if (min(unique_pvalues) > 0){
  		# I think fdrtool neglects this borderline case and this causes returned object
  		# to be discontinuous hence also not concave
  		unique_pvalues <- c(0,unique_pvalues)
  		ecdf_values   <- c(0, ecdf_values)
  	}
  	ll <- fdrtool::gcmlcm(unique_pvalues, ecdf_values, type="lcm")
	ll$length <- length(ll$slope.knots)
	ll$x.knots <- ll$x.knots[-ll$length]
  	ll$y.knots <- ll$y.knots[-ll$length]
	if (!quiet) message(paste("Grenander fit with", ll$length, "knots."))
	ll
}