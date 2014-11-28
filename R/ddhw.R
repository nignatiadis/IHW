#' ddhw: Data-Driven Hypothesis Weights
#'
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param groups  Vector which assigns a group to each p-value.
#' @param filter_statistic  Vector which contains the filter-statistic (covariate, independent under the H0) for each p-value
#' @param nbins  Integer, number of groups into which p-values will be split based on filter_statistic.
#' @param alpha   Numeric, sets the nominal level for FDR control.
#' @param mtests   Numeric, true number of comparisons, must be at least ‘length(unadj_p)
#' @param optim_method  Optimization method with which optimal weights are determined. Available options are "MILP" (default) and "subplex".
#' @param solver   Solver to be used for "MILP" optimization. Available options are "Rsymphony" (default, open source) and "gurobi" (commercial).
#' @param regularization_term   Numeric, biases weights to have the same magnitude, only available for optim_method=="subplex".
#' @param optim_pval_threshold   Numeric between 0 and 1 or Character "auto" (default). P-values above this threshold get excluded
#'    from the optimization procedure (e.g. because it is known a-priori that P-values above that threshold will
#'    not get rejected), thus making it faster. 
#' @param time_limit   Numeric, sets maximum time limit for MILP solvers. Defaults to Inf.


ddhw_grouped <- function(unadj_p, groups, alpha, 
		mtests = length(unadj_p),
		optim_method = "MILP",
		solver = "Rsymphony", #only relevant if optim_method == "MILP"
		regularization_term = 0,
		optim_pval_threshold = "auto",
		time_limit=Inf,
		node_limit=Inf) # only relevant if optim_method == "MILP"
{

	groups <- as.factor(groups) #later: check if factor actually has correct levels..

	if (!check_uniform_groups(groups)){
		warning(paste0("Software and algorithm has only been tested with groups",
			" of equal size."))
	}

	nbins <- length(levels(groups))

	if (optim_pval_threshold == "auto"){
		# heuristic choice, need to check it afterwards
		optim_pval_threshold <- min(1, 
			max(get_bh_threshold(unadj_p, alpha, mtests=mtests)*nbins*2,
				min(unadj_p)*nbins*2) + .Machine$double.eps)
	}


	if (nbins == 1){
		message("Only 1 group supplied, this reduces to the classic Benjamini-Hochberg method.")
		ws <- 1
		weighted_pvalues <- unadj_p
		ts <- get_bh_threshold(unadj_p, alpha, mtests=mtests)

	} else if (optim_method == "MILP"){
		res <- ddhw_grouped_milp(unadj_p, groups, alpha, 
					mtests, optim_pval_threshold, solver=solver,
					time_limit = time_limit, node_limit=node_limit)
		solver_information <- res$solver_information
		ws  <- res$ws
		ts  <- res$ts
		weighted_pvalues <- mydiv(unadj_p, ws[groups])

	} else if (optim_method == "subplex"){
		solver_information <- list(solver="subplex")
		ws <- ddhw_grouped_subplex(unadj_p, groups, alpha, 
					mtests, regularization_term, optim_pval_threshold)
		weighted_pvalues <- mydiv(unadj_p, ws[groups])
		ts <- ws*get_bh_threshold(weighted_pvalues, alpha, mtests=mtests)
	}
	
	# check to see if heuristic optim_pval_threshold was too stringent
	# recursively call function again with higher threshold!
	# seems to be faster to start with low threshold and then grow than start with huge MILP problem
	# unfortunately this is not always true..
	if (any(ts >= 1/2 * optim_pval_threshold)){
		optim_pval_threshold <- min(1,optim_pval_threshold*2)
		message(paste0("rerunning with higher optim_pval_threshold (", optim_pval_threshold, ")"))
		ddhw_res <- ddhw_grouped(unadj_p, groups, alpha, mtests=mtests,
			optim_method = optim_method, solver=solver,
			regularization_term = regularization_term,
			optim_pval_threshold <- optim_pval_threshold)
		return(ddhw_res)
	}

	
	adj_p <- p.adjust( weighted_pvalues, method = "BH", n = mtests)

	df <- data.frame(pvalue = unadj_p, adj_p = adj_p, group= groups)
	# maybe introduce more complicated output object later, but for now, just adjp should be fine!
	new("ddhw",
		 	df = df,
		 	weights = ws,
		 	thresholds = ts,
		 	alpha = alpha,
		 	solver_information = solver_information)
}

#' @rdname ddhw_grouped

ddhw <- function(unadj_p, filter_statistic, nbins, alpha,...){
	groups <- groups_by_filter(filter_statistic, nbins)
	ddhw_object <- ddhw_grouped(unadj_p, groups, alpha, ...)
	filter_statistics(ddhw_object) <- filter_statistic
 	ddhw_object
}


# ddhw based on MILP formulation
# (should) detect global optima
ddhw_grouped_milp <- function(unadj_p, groups, alpha, mtests, 
					optim_pval_threshold, solver="Rsymphony",
					 time_limit=Inf, node_limit=Inf){

	m <- length(unadj_p)
	nbins <- length(levels(groups))

	pvals_list <- split(unadj_p, groups)
	m_groups <- sapply(pvals_list,length)


	if (optim_pval_threshold < 1){
		# keep at least 2 p-values in each group, so we do not have to add special cases for n=0 or n=1
		pvals_list <- lapply(pvals_list, function(p) sort(c(p[p <= max(sort(p,partial=2)[2],optim_pval_threshold)])))
	} else {
		pvals_list <- lapply(pvals_list, function(p) sort(p))
	}



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
	p_lengths <- sapply(pvals_list,length)
	prev_blocks <- 0:(nbins-1)
	starts <- cumsum(c(1,p_lengths[-nbins]-1))

	full_triplet <- do.call(rbind, mapply(single_block_triple_representation, starts,p_lengths, prev_blocks, SIMPLIFY=F))
	full_triplet <- slam::simple_triplet_matrix(full_triplet[,1],full_triplet[,2], full_triplet[,3])

	
	# add final constraint to matrix which corresponds to control of plug-in FDR estimate
	plugin_fdr_constraint <- alpha - unlist(mapply(function(p, m_group) diff(c(0,p))* m_group * mtests/m, pvals_list, m_groups))
	full_triplet <- rbind(full_triplet, matrix(plugin_fdr_constraint, nrow=1))


	z_vars <- ncol(full_triplet) #number of binary variables

	obj        <- rep(1, z_vars) # maximize number of rejections

	if (solver=="Rsymphony") {
		if (is.infinite(time_limit)) time_limit <- -1
		if (is.infinite(node_limit)) node_limit <- -1
		res<- Rsymphony::Rsymphony_solve_LP(obj, full_triplet, rep(">=", nrow(full_triplet)), rep(0, nrow(full_triplet)), 
    		types = "B",
			max = T, verbosity = -2, time_limit = time_limit,
			node_limit = node_limit, gap_limit = -1, first_feasible = FALSE)
		sol <- res$solution
		solver_status <- res$status
	} else if (solver=="gurobi"){

		# keep code for commercial solver for now
		#speed up compared to Symphony appears to range about factor ~2-10, depending on problem
		require("gurobi")
		require("Matrix")

		# convert simple_triplet_matrix of pkg Slam to sparse matrix of pkg Matrix
		# gurobi has problems with simple_triplet_matrix for an (undocumented) reason (?bug)
		full_triplet <-  sparseMatrix(i=full_triplet$i, j=full_triplet$j, x=full_triplet$v,
           dims=c(full_triplet$nrow, full_triplet$ncol))

		model <- list()
		model$A <- full_triplet
		model$obj <- obj
		model$modelsense <- "max"
		model$rhs        <- 0
		model$lb         <- 0
		model$ub         <- 1
		model$sense      <- '>'
		model$vtype      <- 'B'

		params <- list(OutputFlag = 0)
		if (is.finite(time_limit)) params$TimeLimit <- time_limit
		if (is.finite(node_limit)) params$NodeLimit <- node_limit
		res <- gurobi(model, params)
		sol <- res$x
		solver_status <- ""
	} else {
		stop("Only gurobi and Rsymphony solvers are currently supported.")
	}

	# next try to extract the thresholds/weights from solver output
	lengths <- cumsum(sapply(pvals_list, length))
	start_lengths <- c(1, lengths[-nbins]+1)

	pidx_list <- mapply(function(start,end) sol[start:end], start_lengths,lengths, SIMPLIFY=F)

	get_threshold <- function(plist, pidx){
		max_idx <- which.max(which(pidx==1))
		ifelse(length(max_idx)==0, .0, plist[max_idx])
	}

	t_thresholds <- mapply(get_threshold, pvals_list, pidx_list)
	ws <-  if (all(t_thresholds == .0)){
				rep(1,nbins)
			} else { 
				t_thresholds*m/sum(m_groups*t_thresholds)
			}

	list(ts = t_thresholds, ws = ws, solver_information = list(solver=solver, solver_status = solver_status))
}




# ddhw with direct optimization of objective with spherical coordinates
# usually detects high quality local optima [but not neccessarily the global optimum]
ddhw_grouped_subplex <- function(unadj_p, groups, alpha, mtests, regularization_term, optim_pval_threshold){
	nbins <- length(levels(groups))
	m <- length(unadj_p)
	top_unadj_p_idx <- which(unadj_p <= optim_pval_threshold) 
	top_unadj_p <- unadj_p[top_unadj_p_idx]
	top_unadj_p_groups <- groups[top_unadj_p_idx]

	if (regularization_term == 0){
		optim_fun <-  function(w){
    		-sum(  p.adjust( mydiv(top_unadj_p, w[top_unadj_p_groups]), 
    			method="BH", n = mtests) < alpha )
    	}
  	} else {
  		#here we introduce graph-lasso-like regularization term
    	#could also try to just measure deviation from uniformity!

  		optim_fun <-  function(w){
    		-sum(  p.adjust( mydiv(top_unadj_p, w[top_unadj_p_groups]), 
    			method="BH", n = mtests) < alpha) + regularization_term*sum(abs(diff(w)))

    	}	
  	}


  	if (nbins == 2){ # this reduces to a one-dimensional problem which can be handled much easier
		optim_fun_1d <- function(w1){
			optim_fun(c(w1,max(2-w1,0))) #max(2-w1,0) avoid negative weights due to numerical instability
		}
		w1 <- optimize(optim_fun_1d, c(0,2))$minimum
		ws <- c(w1,max(2-w1,0))
	} else {  #interesting part..
		opt <-subplex::subplex(
        	par = w2phi(rep(1, nbins)),
        	fn = function(x) optim_fun(phi2w(x)),
            control = list(trace=0))
     	ws <- phi2w(opt$par)
	}

	ws
}