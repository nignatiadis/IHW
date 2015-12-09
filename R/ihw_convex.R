#' ihw: Independent Hypothesis Weighting
#'
#' Given a vector of p-values, a vector of filter-statistics which are independent of the p-values under the null hypothesis and
#' a nominal significance level alpha, ihw learns multiple testing weights and then applies the weighted Benjamini Hochberg 
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
#' @param nsplits_internal  Integer, how many times to repeat the nfolds_internal splitting. Can lead to better regularization
#'              parameter selection but makes ihw a lot slower.
#' @param lambdas  Numeric vector which defines the grid of possible regularization parameters.
#'				Use "auto" for automatic selection.
#' @param seed Integer or NULL. Split of hypotheses into folds is done randomly. To have output of the function be reproducible, 
#'	we set a seed. Use NULL if you don't want a seed.
#' @param distrib_estimator  Character ("grenander" or "ECDF"). Only use this if you know what you are doing. ECDF with nfolds > 1
#'              or lp_solver == "lpsymphony" will in general be excessively slow, except for very small problems.
#' @param lp_solver  Character ("lpsymphony" or "gurobi"). Internally, ihw solves a sequence of linear programs, which
#'        can be solved with either of these solvers.
#' @param return_internal Returns a lower level representation of the output (only useful for debugging purposes).
#' @param ... Arguments passed to internal functions.
#'
#' @return A ihwResult object.
#' @seealso ihwResult, plot_ihw
#'
#' @examples
#'
#' set.seed(1)
#' X   <- runif(20000, min=0.5, max=4.5) #covariate
#' H   <- rbinom(20000,1,0.1)            #hypothesis true or false
#' Z   <- rnorm(20000, H*X)              #Z-score
#' pvalue <- 1-pnorm(Z)                  #pvalue
#' ihw_res <- ihw(pvalue, X, .1)
#'
#'
#' @export
ihw <- function(pvalues, filter_statistics, alpha,
						filter_statistic_type = "ordinal",
						nbins = "auto",
						quiet =TRUE ,
						nfolds = 5L,
						nfolds_internal = 5L,
						nsplits_internal=1L,
						lambdas = "auto",
						seed = 1L,
						distrib_estimator = "grenander",
						lp_solver="lpsymphony",
						return_internal=FALSE,
						...){

	# This function essentially wraps the lower level function ihw_internal
	# e.g. takes care of NAs, sorts pvalues and then
	# returns a nice ihw object

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


	if (filter_statistic_type =="ordinal" & is.numeric(filter_statistics)){
		groups <- as.factor(groups_by_filter(filter_statistics, nbins))
		penalty <- "total variation"
		if (nbins == "auto"){
			nbins <- min(300, floor(length(pvalues)/1500)) # rule of thumb..
		}
	} else if (is.factor(filter_statistics)){
		groups <- filter_statistics
		if (nbins != "auto" & nbins != nlevels(groups)){
			warning("Overwriting manually specified nbins, since it has to equal the number
				of levels for categorical covariates")
		}

		nbins <- nlevels(groups)

		if (filter_statistic_type == "nominal"){
			penalty <- "uniform deviation"
		} else if (filter_statistic_type == "ordinal") {
			penalty <- "total variation"
		}
	} else {
		stop("filter_statistics are not of the appropriate type")
	}

	if ((length(lambdas)== 1) & (lambdas[1] == "auto")){
		# just a few for now until I get warm starts of the LP solvers to work
		lambdas <- c(0, 1, nbins/8, nbins/4, nbins/2, nbins, Inf)
	}
	# once we have groups, check whether they include enough p-values
	if (any(table(groups) < 1000)){
		message("In general, data-driven choice of weights requires at least 1000 p-values per stratum.")
	}

	group_levels <- levels(groups)
	# sort pvalues globally
	order_pvalues <- order(pvalues) #TODO BREAK ties by filter statistic rank
	reorder_pvalues <- order(order_pvalues)

	sorted_groups <- groups[order_pvalues]
	sorted_pvalues <- pvalues[order_pvalues]

	if (!is.null(seed)) set.seed(as.integer(seed)) #seed

	res <- ihw_internal(sorted_groups, sorted_pvalues, alpha, lambdas,
						penalty=penalty,
						quiet=quiet,
						nfolds=nfolds,
						nfolds_internal=nfolds_internal,
						nsplits_internal=nsplits_internal,
						seed=NULL,
						distrib_estimator = distrib_estimator,
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
					 adj_pvalue = adj_pvalues,
				     weight = weights,
				     weighted_pvalue = weighted_pvalues,
					 group= groups,
					 filter_statistic = filter_statistics,
					 fold=as.factor(folds))

	ihw_obj <- new("ihwResult",
		 			df = df,
		 			weights = res$weight_matrix,
		 			alpha = alpha,
					nbins = as.integer(nbins),
					nfolds = nfolds,
					regularization_term = res$fold_lambdas,
					penalty = penalty,
					filter_statistic_type = filter_statistic_type,
					reg_path_information = data.frame(),
		 			solver_information = list())

	ihw_obj

}

# operate on pvalues which have already been sorted and without any NAs
ihw_internal <- function(sorted_groups, sorted_pvalues, alpha, lambdas,
							    penalty="total variation",
							    seed=NULL,
								quiet=TRUE,
								nfolds = 10L,
								nfolds_internal = nfolds,
								nsplits_internal = 1L,
								distrib_estimator = "distrib_estimator",
								lp_solver="lpsymphony",
								debug_flag=FALSE,
								...){

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

		# TODO replace for loop by single call to ihw_convex with warm starts
		if (length(lambdas) > 1){
			rjs  <- matrix(NA, nrow=length(lambdas), ncol=nsplits_internal)
			for (k in 1:length(lambdas)){
				lambda <- lambdas[k]
				for (l in 1:nsplits_internal){

					rjs[k,l] <- ihw_internal(filtered_sorted_groups, filtered_sorted_pvalues, alpha, lambda,
									    	seed=NULL, quiet=quiet,
									    	nfolds=nfolds_internal,
									    	distrib_estimator = distrib_estimator,
									    	lp_solver=lp_solver)$rjs
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
		if (nfolds==1){
			m_groups_holdout_fold <- m_groups
		} else {
			m_groups_holdout_fold <- m_groups - sapply(filtered_split_sorted_pvalues, length)
		}

		if (distrib_estimator=="grenander"){
			res <- ihw_convex(filtered_split_sorted_pvalues, alpha, m_groups_holdout_fold,
						   penalty=penalty, lambda=lambda, lp_solver=lp_solver, quiet=quiet,...)
		} else if (distrib_estimator == "ECDF"){
			res <- ihw_milp(filtered_split_sorted_pvalues, alpha, m_groups_holdout_fold,
				           penalty=penalty, lambda=lambda, lp_solver=lp_solver, ...)
		} else {
			stop("This type of distribution estimator is not available.")
		}

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
ihw_convex <- function(split_sorted_pvalues, alpha, m_groups, 
		penalty="total variation", lambda=Inf, lp_solver="gurobi", quiet=quiet){

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

			absolute_val_constr_matrix_1 <- matrix(rep(-m_groups,nbins), nbins,nbins, byrow=T)
			diag(absolute_val_constr_matrix_1) <- -m_groups + m
			absolute_val_constr_matrix_1 <- cbind( slam::simple_triplet_zero_matrix(nbins, nbins,mode="double"),
												   absolute_val_constr_matrix_1,
										           -diag(nbins))

			# -f - m t_g +  sum_g m_i t_i  <= 0

			absolute_val_constr_matrix_2 <- matrix(rep(+m_groups,nbins), nbins,nbins, byrow=T)
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

	# incorporate the FDR constraint
	fdr_constr<- matrix(c(rep(-alpha,nbins)*m_groups, rep(1,nbins)*m_groups, rep(0,ncol(constr_matrix)-2*nbins)), nrow=1)
	constr_matrix <- rbind(constr_matrix, fdr_constr)
	nvars <- ncol(constr_matrix)
	rhs <- c(rhs,0)

	if (!quiet) message("Starting to solve LP.")

	if (lp_solver == "gurobi"){

		if (!requireNamespace("gurobi", quietly=TRUE)){
			stop("Gurobi solver appears not be installed. Please use lpsymphony or install gurobi.")
		}

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

	# return matrix as follows: nbins x nlambdas)
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

	full_triplet <- do.call(rbind, mapply(single_block_triple_representation, starts,p_lengths, prev_blocks, SIMPLIFY=F))
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
 		ineq_matrix <-  rbind( cbind(do.call(cbind, mapply(left_ineq_matrix_block, 1:nbins, p_lengths, pdiff_list, SIMPLIFY=F)),
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
			max = T, verbosity = -2, time_limit = time_limit,
			node_limit = node_limit, gap_limit = ifelse(mip_gap <= 10^(-4), -1, mip_gap*100), #now default mip_gap=10^(-4) as in Gurobi while -1 default for Rsymphony
			first_feasible = FALSE)
		sol <- res$solution
		print(str(res))
		solver_status <- res$status

	} else if (lp_solver=="gurobi"){


		if (!requireNamespace("gurobi", quietly=TRUE)){
			stop("Gurobi solver appears not be installed. Please use lpsymphony or install gurobi.")
		}

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

		res <- gurobi::gurobi(model, params)
		sol <- res$x
		solver_status <- res$status

	} else {
		stop("Only gurobi and lpsymphony solvers are currently supported.")
	}

	# next try to extract the thresholds/weights from solver output
	lengths <- cumsum(sapply(split_sorted_pvalues, length))
	start_lengths <- c(1, lengths[-nbins]+1)

	pidx_list <- mapply(function(start,end) sol[start:end], start_lengths,lengths, SIMPLIFY=F)

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
