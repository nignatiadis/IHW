#' ddhw: Data-Driven Hypothesis Weights
#'
#'
#' @param pvalues  Numeric vector of unadjusted p-values.
#' @param filter_statistics  Vector which contains the one-dimensional filter-statistics (covariates, independent under the H0 of the p-value) 
#'                for each test. Can be numeric or a factor. (If numeric it will be converted into factor by binning.)
#' @param alpha   Numeric, sets the nominal level for FDR control.
#' @param filter_statistic_type  "ordinal" or "nominal" (i.e. whether filter statistics can be sorted in increasing order or not)
#' @param nbins  Integer, number of groups into which p-values will be split based on filter_statistic. Use "auto" for
#'             automatic selection of the number of bins.
#' @param mtests   Numeric, true number of comparisons, must be at least ‘length(unadj_p) !NOT tested yet
#' @param optim_method  Optimization method with which optimal weights are determined. Available options are "MILP" (default) and "subplex".
#' @param lp_relaxation If optim_method = "MILP" then one can either solve the NP hard MILP problem (lp_relaxation = F) or the linear relaxation (lp_relaxation=T). The latter is suggested!
#' @param pi0_adaptive  Enhance DDHW to an alpha-exhaustive procedure. !NOT TESTED YET
#' @param solver   Solver to be used for "MILP" optimization. Available options are "Rsymphony" (default, open source) and "gurobi" (commercial, free for academics, faster).
#' @param regularization_term   Numeric
#' @param penalty  "total variation" for ordered data, "uniform deviation" for categorical data
#' @param optim_pval_threshold   Numeric between 0 and 1 or Character "auto" (default). P-values above this threshold get excluded
#'    from the optimization procedure (e.g. because it is known a-priori that P-values above that threshold will
#'    not get rejected), thus making it faster. Defaults to "auto" which will use heuristics to calculate this threshold.
#' @param time_limit   Numeric, sets maximum time limit for MILP solvers. Defaults to Inf.
#' @param node_limit   Integer, sets maximum number of branch and bound nodes to explore for MILP solvers. Defaults to Inf.
#' @param threads   Integer, number of threads for the MILP solver to use. Defaults to 0 (i.e. solver default, often uses all threads available)
#' @param mip_gap_abs  Gurobi specific parameter
#' @param mip_gap      Gurobi specific parameter
#' @param MIPFocus     Gurobi specific parameter
#' @param local_fdr    Combine DDHW ideas with Cai's estimator. !NOT IMPLEMENTED 

ddhw_tmp <- function(pvalues, filter_statistics, alpha,
						filter_statistic_type = "ordinal",
						nbins = 10,
						nfolds = 5,
						nfolds_internal = 4,
						lambdas = "auto",
						lp_solver="lpsymphony",
						return_internal=FALSE,
						...){

	# This function essentially wraps the lower level function ddhw_internal
	# e.g. takes care of NAs, sorts pvalues and then
	# returns a nice ddhw object

	# code from R's p.adjust function to gracefully handle NAs
	nm <- names(pvalues)
    pvalues <- as.numeric(pvalues)

   	p0 <- setNames(pvalues, nm)

   	if (all(nna <- !is.na(pvalues))){
        nna <- TRUE
   	}

   	# filter our p-values
    pvalues <- pvalues[nna]
    filter_statistics <- filter_statistics[nna]
    weights <- rep(NA, length(pvalues))

	if (filter_statistic_type =="ordinal" & is.numeric(filter_statistics)){
		groups <- groups_by_filter(filter_statistics, nbins)
	} else if (is.factor(filter_statistics)){
		groups <- filter_statistics
		nbins <- nlevels(groups)
		if (filter_statistic_type == "nominal"){
			lambdas <- Inf
		}
		# TODO: check if each group has enough p-values
	} else {
		stop("filter_statistics are not of the appropriate type")
	}

	# once we have groups, check whether they include enough p-values

	# sort pvalues globally
	order_pvalues <- order(pvalues)
	reorder_pvalues <- order(order_pvalues)

	sorted_groups <- groups[order_pvalues]
	sorted_pvalues <- pvalues[order_pvalues]

	res <- ddhw_internal(sorted_groups, sorted_pvalues, alpha, lambdas,
						nbins=nbins, nfolds=nfolds, 
						nfolds_internal=nfolds_internal,
						lp_solver=lp_solver,
						...)

	if (return_internal){
		return(res)
	}

	# resort back to original form
	weights <- fill_nas_reorder(res$sorted_weights, nna, reorder_pvalues)
	pvalues <- p0
	weighted_pvalues <- fill_nas_reorder(res$sorted_weighted_pvalues, nna, reorder_pvalues)
	adj_pvalues <- fill_nas_reorder(res$sorted_adj_p, nna, reorder_pvalues)
	groups     <- fill_nas_reorder(res$sorted_groups, nna, reorder_pvalues)
	folds      <- fill_nas_reorder(res$sorted_folds, nna, reorder_pvalues)
	#weights[nna] <- res$sorted_weights[order_pvalues]

	#weighted_pvalues <- setNames(pvalues, nm)
	#weighted_pvalues <- 
	lst <- list(rjs=res$rjs,
	            pvalues=pvalues,
	            weights = weights,
	            weighted_pvalues=weighted_pvalues,
	            adj_pvalues=adj_pvalues,
	            groups=groups,
	            folds=folds)

	lst
	#weighted_pvalues <- mydiv(pvalues, weights)
	#adj_p <- p.adjust( weighted_pvalues, method = "BH")
	#return(list(pvalues=pvalues, adj_p=adj_p,ws=weights,filter_statistics=filter_statistics, groups=groups, folds=folds))
}

# operate on pvalues which have already been sorted and without any NAs
ddhw_internal <- function(sorted_groups, sorted_pvalues, alpha, lambdas,
								nbins = 10,
								nfolds = 10,
								nfolds_internal = nfolds,
								lp_solver="lpsymphony",
								debug_flag=F){

	# TODO: check if lambdas are sorted the way I want them to

	m <- length(sorted_pvalues)
	split_sorted_pvalues <- split(sorted_pvalues, sorted_groups)
	m_groups <- sapply(split_sorted_pvalues, length)

	#recover weights

	# afterwards do the k-fold strategy
	set.seed(1)
	sorted_folds <- sample(1:nfolds, m, replace = TRUE)

	sorted_weights <- rep(NA, m)
	fold_lambdas <- rep(NA, nfolds)

	if (debug_flag==T){
		rjs_path_mat <- matrix(NA, nrow=nfolds, ncol=length(lambdas))
	}
	for (i in 1:nfolds){

		# don't worry about inefficencies right now, they are only O(n) and sorting dominates this
		filtered_sorted_groups <- sorted_groups[sorted_folds!=i]
		filtered_sorted_pvalues <- sorted_pvalues[sorted_folds!=i]
		filtered_split_sorted_pvalues <- split(filtered_sorted_pvalues, filtered_sorted_groups)

		# within each fold do iterations to also find lambda

		# TODO replace for loop by single call to ddhw_convex..
		if (length(lambdas) > 1){
			rjs  <- rep(NA, length(lambdas))
			for (k in 1:length(lambdas)){
				lambda <- lambdas[k]
				# To consider: internal nfolds does not have to be the same, could be 2 for speedup
				rjs[k] <- ddhw_internal(filtered_sorted_groups, filtered_sorted_pvalues, alpha, lambda,
										nbins=nbins, nfolds=nfolds_internal, lp_solver=lp_solver)$rjs
				if (debug_flag==T){
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
		m_groups_holdout_fold <- m_groups - sapply(filtered_split_sorted_pvalues, length)

		res <- ddhw_convex(filtered_split_sorted_pvalues, alpha, m_groups_holdout_fold,
						   lambda=lambda, lp_solver=lp_solver)

		sorted_weights[sorted_folds == i] <- res$ws[sorted_groups[sorted_folds==i]]
	# I do need some kind of new object though...
	}

	sorted_weighted_pvalues <- mydiv(sorted_pvalues, sorted_weights)
	sorted_adj_p <- p.adjust(sorted_weighted_pvalues, method = "BH")
	rjs   <- sum(sorted_adj_p <= alpha)
	lst <- list(lambda=lambda, fold_lambdas=fold_lambdas, rjs=rjs, sorted_pvalues=sorted_pvalues,
					sorted_weighted_pvalues = sorted_weighted_pvalues,
					sorted_adj_p=sorted_adj_p, sorted_weights=sorted_weights,
					sorted_groups=sorted_groups, sorted_folds=sorted_folds)
	if (debug_flag == T) {
		lst$rjs_path_mat <- rjs_path_mat
	}
	lst
}

ddhw_convex <- function(split_sorted_pvalues, alpha, m_groups, lambda=Inf, lp_solver="gurobi"){

	# preprocessing:  Set too low p-values to 0, otherwise LP solvers have problems
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
	grenander_list <- lapply(split_sorted_pvalues, presorted_grenander)


	#set up LP
	# TODO account for different stratum sizes
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
	fdr_constr<- matrix(c(rep(-alpha,nbins), rep(1,nbins), rep(0,ncol(constr_matrix)-2*nbins)), nrow=1)
	constr_matrix <- rbind(constr_matrix, fdr_constr)
	nvars <- ncol(constr_matrix)
	rhs <- c(rhs,0)

		if (lp_solver == "gurobi"){
		model <- list()
		model$A <- sparseMatrix(i=constr_matrix$i, j=constr_matrix$j, x=constr_matrix$v,
           		dims=c(constr_matrix$nrow, constr_matrix$ncol))
		model$obj <- obj
		model$modelsense <- "max"
		model$rhs        <- rhs
		model$lb         <- 0
		model$ub         <- 2
		model$sense      <- '<'

		params <- list(OutputFlag=1, Method=2)
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
			max = T, verbosity = -2, first_feasible = FALSE)
		sol <- res$solution
		solver_status <- res$status
	} else {
		stop("Only gurobi and lpsymphony solvers currently supported.")
	}

	# catch negative weights due to numerical rounding
	ts <- pmax(sol[(1:nbins)+nbins],0)
	ws <- thresholds_to_weights(ts, m_groups)	#ts/sum(ts)*nbins



	return(list(ws=ws))
}


presorted_grenander <- function(sorted_pvalues){
  	n  <- length(sorted_pvalues)
  	unique_pvalues <- unique(sorted_pvalues)
  	ecdf_values <- cumsum(tabulate(match(sorted_pvalues, unique_pvalues)))/n
  	if (min(unique_pvalues) > 0){
  		# I think fdrtool neglects this borderline case and this causes returned object
  		# to be discontinuous hence also not concave
  		unique_pvalues <- c(0,unique_pvalues)
  		ecdf_values   <- c(0, ecdf_values)
  	}
  	ll <- gcmlcm(unique_pvalues, ecdf_values, type="lcm")
	ll$length <- length(ll$slope.knots)
	ll$x.knots <- ll$x.knots[-ll$length]
  	ll$y.knots <- ll$y.knots[-ll$length]
	ll
}