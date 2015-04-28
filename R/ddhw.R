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
		bh_rejections = "auto",
		rejections_bound = NULL,
		local_fdr = F,
		optim_method = "MILP",
		solver = "Rsymphony", #only relevant if optim_method == "MILP"
		regularization_term = Inf,
		penalty="total variation",
		optim_pval_threshold = "auto",
		time_limit=Inf,
		node_limit=Inf,
		threads=0,
		solution_limit=Inf,
		mip_gap_abs=Inf,
		mip_gap = 10^(-4),
		MIPFocus=0,
		lp_relaxation=FALSE) # only relevant if optim_method == "MILP"
{

	groups <- as.factor(groups) #later: check if factor actually has correct levels..

	nbins <- length(levels(groups))

	if (optim_pval_threshold == "auto"){
		# heuristic choice, need to check it afterwards
		t_tmp <- get_bh_threshold(unadj_p, alpha, mtests=mtests)
		if (bh_rejections == "auto")  bh_rejections <- sum(unadj_p <= t_tmp)
		optim_pval_threshold <- min(1,
			max(t_tmp*nbins*2,
				min(unadj_p)*nbins*2) + .Machine$double.eps)
	} else {
		if (bh_rejections == "auto") bh_rejections <- sum(p.adjust(unadj_p, method="BH") < alpha)
	}

	t_bh <- get_bh_threshold(unadj_p, alpha ,mtests=mtests)

	# simple cases in which we can just use BH rather than MILP
	if (nbins == 1 || regularization_term==0){

		if (nbins == 1){
			message("Only 1 group supplied, this reduces to the classic Benjamini-Hochberg method.")
		}
		ws <- 1
		weighted_pvalues <- unadj_p
		ts <- get_bh_threshold(unadj_p, alpha, mtests=mtests)
		solver_information <- list("R p.adjust function")

		if (regularization_term==0){
			message("Regularization term equal to 0, this reduces to the classic Benjamini-Hochberg method.")
			ws <- rep(1, nbins)
			ts <- rep(ts, nbins)
		}

	} else if (optim_method == "MILP"){
		res <- ddhw_grouped_milp(unadj_p, groups, alpha,
					mtests, optim_pval_threshold,
					regularization_term,
					penalty=penalty,
					local_fdr=local_fdr,
					t_bh=t_bh,
					bh_rejections=bh_rejections, rejections_bound=rejections_bound,
					solver=solver,
					time_limit = time_limit, node_limit=node_limit,
					solution_limit=solution_limit,
					mip_gap_abs=mip_gap_abs,
					mip_gap = mip_gap,
					threads=threads,
					MIPFocus=MIPFocus,
					lp_relaxation=lp_relaxation)
		solver_information <- res$solver_information
		ws  <- res$ws
		ts  <- res$ts		

		weighted_pvalues <- mydiv(unadj_p, ws[groups])

		# instead of using the threshold naturally defined by the procedure, rerun BH with weights returned.
		# This ensures that IF a suboptimal solution was found
		# 1) due to numerical instability
		# 2) because the user set node_limit or time_limit or mip_gap options
		# the BH procedure will be able to improve upon this solution!

		ts <- ws*get_bh_threshold(weighted_pvalues, alpha)
		
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
	if ((!optim_pval_threshold==1) && any(ts >= 1/2 * optim_pval_threshold)){
		optim_pval_threshold <- min(1,optim_pval_threshold*2)
		message(paste0("rerunning with higher optim_pval_threshold (", optim_pval_threshold, ")"))
		ddhw_res <- ddhw_grouped(unadj_p, groups, alpha, mtests=mtests,
			local_fdr=local_fdr,
			optim_method = optim_method, solver=solver,
			regularization_term = regularization_term,
			optim_pval_threshold = optim_pval_threshold,
			time_limit = time_limit, node_limit = node_limit,
			threads = threads,
			solution_limit=solution_limit, 
			mip_gap_abs=mip_gap_abs,
			mip_gap = mip_gap,
			MIPFocus=MIPFocus,
			lp_relaxation=lp_relaxation,
			bh_rejections=bh_rejections, rejections_bound=rejections_bound)
		return(ddhw_res)
	}


	names(ws) <- levels(groups)

	adj_p <- p.adjust( weighted_pvalues, method = "BH", n = mtests)

	df <- data.frame(pvalue = unadj_p, adj_p = adj_p, group= groups)
	# maybe introduce more complicated output object later, but for now, just adjp should be fine!
	new("ddhw",
		 	df = df,
		 	weights = ws,
		 	thresholds = ts,
		 	alpha = alpha,
			nbins = as.integer(nbins),
			regularization_term = regularization_term,
			penalty = penalty,
			reg_path_information = data.frame(),
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
					 	optim_pval_threshold, regularization_term,
					 	penalty="total variation",
					 	local_fdr=F,
					 	t_bh=0, bh_rejections =0 , rejections_bound=NULL,
					 	solver="Rsymphony",
					 	time_limit=Inf, node_limit=Inf,
					 	threads=0,
					 	solution_limit=solution_limit,
					 	mip_gap_abs=mip_gap_abs,
					 	mip_gap = mip_gap,
					 	MIPFocus=MIPFocus,
					 	lp_relaxation=lp_relaxation){

	m <- length(unadj_p)
	nbins <- length(levels(groups))

	pvals_list <- split(unadj_p, groups)
	m_groups <- sapply(pvals_list,length) # number of pvals in each group (BEFORE applying optim_pval_threshold)


	if (optim_pval_threshold < 1){
		# keep at least 2 p-values in each group, so we do not have to add special cases for n=0 or n=1
		pvals_list <- lapply(pvals_list, function(p) sort(c(p[p <= max(sort(p,partial=2)[2],optim_pval_threshold)])))
	} else {
		pvals_list <- lapply(pvals_list, function(p) sort(p))
	}


	pmax <- max(unlist(pvals_list))
	pvals_list <- lapply(pvals_list, function(p) p/pmax)
	#ts_bh <- sapply(pvals_list, function(ps) max(0, ps[ps <= t_bh]))


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
	diff_coefs <- unlist(mapply(function(p, m_group) diff(c(0,p))* m_group , pvals_list, m_groups))
	#plugin_fdr_constraint <- (alpha - diff_coefs * mtests/m)/m #alpha - diff_coefs * mtests/m
	plugin_fdr_constraint <- (alpha/pmax - diff_coefs * mtests/m)/m #alpha - diff_coefs * mtests/m

	full_triplet <- rbind(full_triplet, matrix(plugin_fdr_constraint, nrow=1))

	z_vars <- ncol(full_triplet) #number of binary variables

	obj        <- rep(1, z_vars) # maximize number of rejections

	if (is.finite(regularization_term)){
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

		pdiff_list <- lapply(pvals_list, function(p) diff(c(0,p)))

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
 			regularization_constraint <- matrix(c(rep(0,z_vars), regularization_term/m*m_groups, rep(-1, nbins-1)), nrow=1)

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
 			regularization_constraint <- matrix(c(rep(0,z_vars), regularization_term*m_groups, rep(-1, nbins)), nrow=1)

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
   
	if (solver=="Rsymphony") {

		rsymphony_ub <- list(upper = list(ind = 1:ncol(full_triplet), val = model_ub))

		if (is.infinite(time_limit)) time_limit <- -1
		if (is.infinite(node_limit)) node_limit <- -1
		res<- Rsymphony::Rsymphony_solve_LP(obj, full_triplet, rep(">=", nrows), 
			model_rhs,
    		types = model_vtype, bounds= rsymphony_ub,
			max = T, verbosity = -2, time_limit = time_limit,
			node_limit = node_limit, gap_limit = ifelse(mip_gap <= 10^(-4), -1, mip_gap*100), #now default mip_gap=10^(-4) as in Gurobi while -1 default for Rsymphony
			first_feasible = FALSE)
		sol <- res$solution
		print(str(res))
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

	#get_threshold2 <- function(pdiff, pidx){
	#	sum(pdiff*pidx)
	#}

	t_thresholds1 <- mapply(get_threshold, pvals_list, pidx_list)     # to be stored in ddhw object
	#t_thresholds2 <- mapply(get_threshold2, pdiff_list, pidx_list)
	t_thresholds <- if (is.infinite(regularization_term)){
					 	t_thresholds1
					} else {
						sol[(z_vars+1):(z_vars+nbins)]  # to be used for calculating weights
					}

    ####################################################################################################
    # TODO!!!! I *NEED* to check t_thresholds1 vs t_thresholds, easy way to see if the solver screwed up
    # numerically...In particular three simple checks: 
    #
    # 1) t_thresholds >= t_thresholds1
    # 2) Plugin FDR has to be controlled (with respect to both t_thresholds, t_thresholds2)
    # 3) total variation of weights has to be less than lambda
    #
    # Also what to do if the solver did *screw* up? Maybe restart solver with different numerical
    # tolerances? Easier for now: If this happens, return BH procedure... Better: Avoid numerical  
    # issues by better model scaling. Maybe also have a second simpler linear program which uses
    # previous results as starting point to make sure constraints are satisfied.
    #
    # For now deal with this with special casing regularization_term == 0 case by having it use BH, since
    # this leads to speedup anyway, and this case is particularly bad for SYMPHONY it seems like.
    #
	# print(t_thresholds)
	# print(t_thresholds1)
	####################################################################################################

	ws <-  if (all(t_thresholds1 == .0) || all(t_thresholds == .0)){
				rep(1,nbins)
			} else {
				t_thresholds*m/sum(m_groups*t_thresholds) #might need adjustment if mtests > m
			}

	t_thresholds <- t_thresholds*pmax

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

bootstrapped_weight <- function(obj) {
  groups <- groups_factor(obj)
	pvals <- pvalues(obj)
	for (lvl in levels(groups)){
		idx <- which(groups == lvl)
		sampled_idx <- sample(idx, replace=T)
		pvals[idx] <- pvals[sampled_idx]
	}
	weights(ddhw_grouped(pvals, groups, alpha(obj)))
}

per_bin_fdrs <- function(obj) {
	ts <- thresholds(obj)
	groups <- groups_factor(obj)
	pvals <- pvalues(obj)
	pv_list <- split(pvals, groups)
	print(sapply(pv_list,length))
	fdrs <- mapply(function(t,pvec) length(pvec)*t/max(1, sum(pvec <= t)), ts, pv_list)
	return(fdrs)
}

ddhw_reg_path <- function(unadj_p, filterstat, nbins, alpha, reg_pars, nreps, optim_pval_threshold="auto"){
  rjs <- rep(NA, length(reg_pars))
  ws_mat  <- matrix(NA, nrow = nbins, ncol= length(reg_pars))
  cv_df <- list()

  m <- length(unadj_p)
  all_groups <- groups_by_filter(filterstat,nbins)
  train_idxs <-createDataPartition(all_groups, p = .5, list = T, times = nreps)

  for (i in seq_along(reg_pars)){
  	reg_par <- reg_pars[i]
  	ddhw_res <- ddhw(unadj_p, filterstat,nbins, alpha,
  		regularization_term=reg_par, solver="gurobi", optim_pval_threshold=optim_pval_threshold)
  	rjs[i] <- rejections(ddhw_res)
  	ws_mat[,i] <- weights(ddhw_res)


  	rjs_train <- rep(NA,2*nreps)
  	rjs_test  <- rep(NA,2*nreps)
  	fdr_bias  <- rep(NA,2*nreps)
  	fdr_bias_s <- rep(NA,2*nreps)

  	for (n in 1:nreps){
		train_idx <- train_idxs[[n]]
		# previously was doing non-balanced	 splitting, now lets do it balanced!
    	test_idx  <- setdiff(1:m, train_idx)
    	ddhw_obj <- ddhw(unadj_p[train_idx], filterstat[train_idx], nbins, alpha,
      		regularization_term=reg_par, solver="gurobi", optim_pval_threshold=optim_pval_threshold)
    	ws <- weights(ddhw_obj)
    	rjs_train[2*n-1] <- rejections(ddhw_obj)
    	grps <- groups_by_filter(filterstat[test_idx],nbins)
    	rjs_test[2*n - 1] <- sum(p.adjust(mydiv(unadj_p[test_idx],ws[grps]), method="BH") < alpha)
    	fdr_bias[2*n - 1] <- (rjs_train[2*n - 1]-rjs_test[2*n - 1])/max(rjs_test[2*n-1],rjs_train[2*n - 1],1)

    	tmp <- train_idx
    	train_idx <- test_idx
    	test_idx <- tmp
    	ddhw_obj <- ddhw(unadj_p[train_idx], filterstat[train_idx], nbins, alpha,
      		regularization_term=reg_par, solver="gurobi", optim_pval_threshold=optim_pval_threshold)
    	ws <- weights(ddhw_obj)
    	rjs_train[2*n] <- rejections(ddhw_obj)
    	grps <- groups_by_filter(filterstat[test_idx],nbins)
    	rjs_test[2*n] <- sum(p.adjust(mydiv(unadj_p[test_idx],ws[grps]), method="BH") < alpha)
    	fdr_bias[2*n] <- (rjs_train[2*n]-rjs_test[2*n])/max(rjs_test[2*n], rjs_train[2*n],1)

   		fdr_bias_s[2*n - 1] <- (rjs_train[2*n] - rjs_test[2*n-1])/max(rjs_train[2*n],1)
   		fdr_bias_s[2*n]     <- (rjs_train[2*n-1] - rjs_test[2*n])/max(rjs_train[2*n-1],1)
  	}

  	cv_df[[i]] <- data.frame(fdr_bias = fdr_bias, rjs_train = rjs_train, rjs_test= rjs_test, fdr_bias_s = fdr_bias_s)
  }

  reg_path <- list(reg_pars=reg_pars, rjs=rjs, ws=ws_mat, cv_df = cv_df)
  reg_path
}

ddhw_cv <- function(unadj_p, filterstat, nbins, alpha, reg_par, return_weights=F,...){
 	ddhw_res <- ddhw(unadj_p, filterstat,nbins,alpha,  solver="gurobi", regularization_term=reg_par,...)

 	m <- length(unadj_p)
 	o <- order(filterstat)
 	train_idx <- o[seq(1,m, by=2)]
 	test_idx <- setdiff(1:m, train_idx)

 	ddhw_obj <- ddhw(unadj_p[train_idx], filterstat[train_idx], nbins, alpha,
      		regularization_term=reg_par, solver="gurobi", ...)
    ws <- weights(ddhw_obj)
    rjs_train<- rejections(ddhw_obj)
    grps <- groups_by_filter(filterstat[test_idx],nbins)
    rjs_test <- sum(p.adjust(mydiv(unadj_p[test_idx],ws[grps]), method="BH") < alpha)

    tmp <- train_idx
    train_idx <- test_idx
    test_idx <- tmp
    ddhw_obj <- ddhw(unadj_p[train_idx], filterstat[train_idx], nbins, alpha,
      		regularization_term=reg_par, solver="gurobi", ...)
    ws <- weights(ddhw_obj)
    rjs_train2<- rejections(ddhw_obj)
    grps <- groups_by_filter(filterstat[test_idx],nbins)
    rjs_test2 <- sum(p.adjust(mydiv(unadj_p[test_idx],ws[grps]), method="BH") < alpha)

    lst <- data.frame(reg_par=reg_par, rjs=rejections(ddhw_res), rjs_train = rjs_train, rjs_test = rjs_test2,
    				rjs_train2 = rjs_train2, rjs_test2 = rjs_test)

    if (return_weights){
    	ws_df <- data.frame(bin = 1:nbins, weight = weights(ddhw_res), reg_par = reg_par)
    	return(list(main_df=lst, ws_df=ws_df))
    } else {
    	return(lst)
    }
}

# ok how would model selection look like? Also note that by now we would like to also store information about nbins and reg par.

model_selection <- function(ddhw_obj){
	# this should output, a model selection class, don't even bother cross-validating, just one run

}


# determines number of bins and regularization parameter automatically, returns appropriate output
ddhw_auto <- function(unadj_p, filterstat, alpha, holdout_selection = F, 
		nbins="auto", parametric_reps=10, resampling="nonparametric lfdr",
		lfdr_estimate=F, ...){
	# some heuristic rule for determining number of bins
	if (nbins=="auto") nbins <- max(min(20, floor(length(unadj_p)/1000)), 1)
	
	# start by optimizing with what amounts to reg par = Infty
	ddhw_full <- ddhw(unadj_p, filterstat, nbins, alpha, ...)

	# determine maximum necessary regularization
	max_reg <- total_variation(weights(ddhw_full))

	print(paste("max reg:",max_reg))
	# for reg_par = 0 to max_reg     (e.g. 10 values?) do the holdout-experiment
	reg_par_number <- 7
	#reg_pars <- seq(0, max_reg, length=reg_par_number)

	reg_pars <- c(0, rev(sapply(1:reg_par_number, function(x) max_reg/2^(x-1))))
	print(reg_pars)

	m <- length(unadj_p)
 	o <- order(filterstat)
 	train_idx <- o[seq(1,m, by=2)]
 	test_idx <- setdiff(1:m, train_idx)
    grps <- groups_by_filter(filterstat[test_idx],nbins)

 	rjs_test <- rep(NA, reg_par_number)

 	# apply holdout method for decreasing value of the regularization parameter
 	# thus we can keep track of the rejections in the previous step and pass it to 
 	# the solver as a lower bound to speed things up a tiny bit
 	bh_rejections = "auto"

 	if (holdout_selection){
 		for (i in 1:reg_par_number){
			ddhw_obj <- ddhw(unadj_p[train_idx], filterstat[train_idx], nbins, alpha,
      						regularization_term=reg_pars[i], bh_rejections=bh_rejections, ...)
			bh_rejections <- rejections(ddhw_obj) 
			#print(paste("prev rjs:",bh_rejections))
    		ws <- weights(ddhw_obj)
    		#rjs_train<- rejections(ddhw_obj)
    		rjs_test[i] <- sum(p.adjust(mydiv(unadj_p[test_idx],ws[grps]), method="BH") < alpha)
    		print(paste("prev rjs:", rjs_test[i]))
 		}
 		reg_pars <- reg_pars[which.max(rjs_test)]

 		print(reg_pars)
 	} else {
 		reg_pars <- reg_pars[-1]
 	}




 	# prev. idea: 
	# generate (10?20?) simulations with covmod and estimate FDR on reg_par selected from prev. step.
	# while (FDR > 1.01*alpha) and reg_par >= 0.1 try with reg_par/2 instead

	# instead be less greedy: use regularization_term from above example and try it with that
	# if it does not work, revert to lambda = 0

	if ((parametric_reps >0 || lfdr_estimate)){ #&& reg_pars[1]	 > 0){

		if (resampling=="beta-uniform"){
	    	covmod_fit<- covmod(unadj_p, filterstat, nbins)
	    } else if (resampling=="nonparametric lfdr"){
	    	lfdr_fit_obj <- lfdr_fit(unadj_p, groups_factor(ddhw_full))
	    }

	    if (!lfdr_estimate){
	    	sims <- list()

	    	for (i in 1:parametric_reps){
    			if (resampling=="beta-uniform"){
	    			sims[[i]] <- simulate(covmod_fit)
	    		} else if (resampling=="nonparametric lfdr"){
	    			sims[[i]] <- with(lfdr_fit_obj, lfdr_sim(pvalue, lfdr, group))
	    		}
	    	}

    		#FDP <- matrix(NA, parametric_reps, length(reg_pars))

    		for (j in rev(seq_along(reg_pars))){
    			FDP <- rep(NA, parametric_reps)
    			reg_par <- reg_pars[j]
    			print(paste("par",reg_par))

    			for (i in 1:parametric_reps){
    				sim <- sims[[i]]
    				ddhw_param <- ddhw(sim$pvalue, sim$group, nbins, 
    						alpha, regularization_term=reg_par, ...)
    				tmp_idx <- adj_pvalues(ddhw_param) <= alpha
    				print(paste("rejections:", sum(tmp_idx)))
 					FDP[i] <- ifelse(sum(tmp_idx) == 0, 0,
 								mean(1-sim$H[tmp_idx]))
 				}
 				# should not have NAs :(
 				if (any(is.na(FDP))) warning("NAs!!")
 					FDR <- mean(FDP, na.rm=T)
 					print(FDR)
 				if (FDR < alpha){
 					break
 				} else {
 					reg_par <- 0
 				}
 			}
 		} else if (lfdr_estimate) {
 			for (j in rev(seq_along(reg_pars))){
 				reg_par <- reg_pars[j]
 				ddhw_param <- ddhw(unadj_p, filterstat, nbins, 
    						alpha, regularization_term=reg_par, ...)
 				if (rejections(ddhw_param) > 0){
 					FDR <- mean(lfdr_fit_obj$lfdr[adj_pvalues(ddhw_param) <= alpha])
 					if (FDR < alpha) {
 						break
 					} else {
 						reg_par <- 0
 					}
 				}else{
 					break
 				}
 			}	
 		}
 	}


 		#FDR <- colMeans(FDP, na.rm=T) #find out why occasionally get NA..
 		#print(FDR)
 		#FDR_sd <- apply(FDP, 2, sd, na.rm=T)
 		#reg_par <- ifelse(FDR <= alpha, reg_par, 0)
	
	
	#return(list(FDR=FDR,FDR_sd=FDR_sd))
	print(reg_par)

	ddhw(unadj_p, filterstat, nbins, alpha, regularization_term=reg_par, ...)


	# return final result with 0.1..
}

lfdr_fit <- function(unadj_p, group, lfdr_estimation="fdrtool"){
	if (lfdr_estimation == "covmod"){
		lfdrs <- covmod_grouped(unadj_p, group)$cm_fdr
	} else {

  		pvals_list <- split(unadj_p, group)
  		if (lfdr_estimation == "fdrtool"){
  			lfdr_fun <- function(pv) fdrtool(pv, statistic="pvalue",plot=F,verbose=F)$lfdr
  		} else if (lfdr_estimation == "ConcaveFDR"){
			lfdr_fun <- function(pv) ConcaveFDR(pv, statistic="pvalue",plot=F,verbose=F)$lfdr.log
  		} else if (lfdr_estimation == "locfdr"){
  			lfdr_fun <- function(pv) locfdr(qnorm(pv), nulltype=0, plot=0)$fdr
  		} else if (lfdr_estimation == "mixfdr"){
  	  		lfdr_fun <- function(pv) mixFdr(qnorm(pv), theonull=T, plots=F)$fdr
  		}
  	
  		lfdr_list <- lapply(pvals_list, lfdr_fun) 
  		lfdrs <- unsplit(lfdr_list, group)
  	}
  	fit_obj <- data.frame(pvalue=unadj_p, lfdr=lfdrs, group=group)
  	fit_obj
}

cai_bh <- function(unadj_p, filterstat, nbins, alpha, ...){
	grps <- groups_by_filter(filterstat,nbins)
	lfdr_res <- lfdr_fit(unadj_p, grps, ...)
	lfdrs <- lfdr_res$lfdr
	# sort lfdrs, break ties by pvalues so that in the end within each stratum
	# we get monotonic adjusted p-values as a function of the p-values
	# this is mainly needed for grenander based lfdrs
	o <- order(lfdrs, unadj_p)
	lfdrs_sorted <- lfdrs[o]
	fdr_estimate <- cumsum(lfdrs_sorted)/(1:length(unadj_p))
	adj_p <- rev(cummin(rev(fdr_estimate)))
	adj_p[order(o)]
}



lfdr_sim <- function(pvalue, lfdr, group){
	m <- length(pvalue)
	H <- 1 - rbinom(m, 1 , lfdr)
	sim_pvalue <- runif(m)
	sim_pvalue[H==1] <- pvalue[H==1]
	data.frame(pvalue=sim_pvalue, group=group, H=H)
}

reg_path <- function(unadj_p, filterstat, nbins, alpha, reg_pars,...){
	reg_par_number <- length(reg_pars)

	m <- length(unadj_p)
 	o <- order(filterstat)
 	train_idx <- o[seq(1,m, by=2)]
 	test_idx <- setdiff(1:m, train_idx)
    grps <- groups_by_filter(filterstat[test_idx],nbins)

 	rjs_test <- rep(NA, reg_par_number)

 	# apply holdout method for decreasing value of the regularization parameter
 	# thus we can keep track of the rejections in the previous step and pass it to 
 	# the solver as a lower bound to speed things up a tiny bit
 	bh_rejections = "auto"
 	for (i in 1:reg_par_number){
		ddhw_obj <- ddhw(unadj_p[train_idx], filterstat[train_idx], nbins, alpha,
      						regularization_term=reg_pars[i], bh_rejections=bh_rejections, ...)
		bh_rejections <- rejections(ddhw_obj) 
		print(paste("prev rjs:",bh_rejections))
    	ws <- weights(ddhw_obj)
    	#rjs_train<- rejections(ddhw_obj)
    	rjs_test[i] <- sum(p.adjust(mydiv(unadj_p[test_idx],ws[grps]), method="BH") < alpha)
 	}
 	rjs_test
}

crossval_fdr_estimate <- function(unadj_p, groups, alpha, regularization_term=Inf, nreps=10,...){

	m <- length(unadj_p)
	m_train <- floor(m/2)

	fdr_estimate <- rep(NA,nreps)
	
	for (k in 1:nreps){
		  
		train_idx <- sample(1:m, floor(m/2))  
		test_idx  <- setdiff(1:m, train_idx)

		ddhw_obj <- ddhw_grouped(unadj_p[train_idx], groups[train_idx], alpha, 
						regularization_term=regularization_term, ...)
		train_rjs <- rejections(ddhw_obj, method="thresholds", groupwise=T)
		train_ts  <- thresholds(ddhw_obj)
		#train_ts  <- regularize_thresholds(train_ts, stratum_sizes(ddhw_obj), sum(rjs),alpha, m_train)

  		test_pv_list <- split(unadj_p[test_idx], groups[test_idx])
  		test_m_groups <- sapply(test_pv_list,length)
  		
  		test_rjs <- mapply(function(pv,t) sum(pv <= t), test_pv_list, train_ts)
  	
  		test_fdrs <- mapply(my_grenander, test_pv_list, train_ts)

  		if (sum(train_rjs) == 0){
  			fdr_estimate[k] = 0
  		} else {
  			fdr_estimate[k] <- sum(train_rjs*test_fdrs)/sum(train_rjs)
  		}

  		#if (sum(train_rjs) == 0){
  		#	print("lol..")
  		#	fdr_estimate[k] <- 1 # setting this to 1 seems also plausible.
  		#} else {
			#fdr_estimate[k] <- sum(train_rjs*pmin(1,(test_m_groups-test_rjs)*train_ts/test_rjs, na.rm=T))/sum(train_rjs)
		#}
	}

	fdr_estimate
	#list(test_m = test_m_groups, test_rjs=test_rjs,  train_rjs=train_rjs, train_ts = train_ts,
	#	pm = (test_m_groups-test_rjs)*train_ts/test_rjs, fdr_est = fdr_estimate)
}


ddhw_auto2 <- function(unadj_p, filterstat, alpha, pi0="auto",
		nbins="auto", reg_par_number=7, reg_pars=NULL, nreps=10, sd_penalty=2, ...){

	pi0 <- fdrtool(unadj_p, statistic="pvalue", plot=F, verbose=F)$param[3]
	print(pi0)

	# some heuristic rule for determining number of bins
	if (nbins=="auto") nbins <- max(min(20, floor(length(unadj_p)/1000)), 1)
	
	# start by optimizing with what amounts to reg par = Infty

	if (is.factor(filterstat)){
		ddhw_full <- ddhw_grouped(unadj_p, filterstat, alpha, ...)
		groups <- filterstat
		nbins <- length(levels(groups))
		max_reg <- uniform_deviation(weights(ddhw_full))

	} else {
		ddhw_full <- ddhw(unadj_p, filterstat, nbins, alpha, ...)
		groups <- groups_by_filter(filterstat, nbins)
		# determine maximum necessary regularization
		max_reg <- total_variation(weights(ddhw_full))
	}

	print(paste("max reg:",max_reg))
	# for reg_par = 0 to max_reg     (e.g. 10 values?) do the holdout-experiment

	if (is.null(reg_pars)){
		reg_pars <- seq(0, max_reg, length=reg_par_number)
	}

	print(max_reg)
	#reg_pars <- c(0, rev(sapply(1:reg_par_number, function(x) max_reg/2^(x-1))))

	m <- length(unadj_p)
 
	# instead be less greedy: use regularization_term from above example and try it with that
	# if it does not work, revert to lambda = 0


	# generate indices for train/test datasets
	train_idx_list <- replicate(nreps, sample(1:m, floor(m/2)), simplify=F)

    for (j in rev(seq_along(reg_pars))){

    	fdr_estimate_fold1 <- rep(NA, nreps)
    	fdr_estimate_fold2 <- rep(NA, nreps)

    	reg_par <- reg_pars[j]
    	print(paste("par",reg_par))


    	for (k in 1:nreps){

			train_idx <- train_idx_list[[k]]
			test_idx  <- setdiff(1:m, train_idx)

			ddhw_obj <- ddhw_grouped(unadj_p[train_idx], groups[train_idx], alpha, 
						regularization_term=reg_par, ...)

			train_rjs <- rejections(ddhw_obj, method="thresholds", groupwise=T)
			train_ts  <- thresholds(ddhw_obj)

			# merge adjacent groups which have the same weights
			merged_groups <- as.factor(groups)
			merged_levels <- cumsum(c(0,abs(diff(weights(ddhw_obj))) > .Machine$double.eps))
			merged_levels_rle <- rle(merged_levels)$lengths
			levels(merged_groups) <- merged_levels
			print(length(unique(levels(merged_groups))))

			cnt <- 1
			train_ts_merged <- c()
			train_rjs_merged <- c()
			for (i in seq_along(merged_levels_rle)){
				new_cnt <- cnt + merged_levels_rle[i]
				train_ts_merged[i] <- train_ts[cnt]
				train_rjs_merged[i] <- sum(train_rjs[cnt:(new_cnt-1)])
				cnt <- new_cnt
			}

	  		test_pv_list <- split(unadj_p[test_idx], merged_groups[test_idx])
  			test_m_groups <- sapply(test_pv_list,length)
  		
  			#test_rjs <- mapply(function(pv,t) sum(pv <= t), test_pv_list, train_ts)
  	
  			## test_fdrs <- mapply(my_grenander, test_pv_list, train_ts_merged)
			test_fdrs <- mapply(function(pv,t) my_grenander(pv,t,distrib=T), test_pv_list, train_ts_merged, SIMPLIFY=F)
			test_pi0s <- sapply(test_fdrs, function(x) x$pi_zero)
			test_distribs <- sapply(test_fdrs, function(x) x$distrib)
  			
  			if (sum(train_rjs) == 0){
  				fdr_estimate_fold1[k] <- 0
  			} else {
  				##fdr_estimate_fold1[k] <- sum(train_rjs_merged*test_fdrs)/sum(train_rjs_merged)
  				fdr_estimate_fold1[k] <- min(1,sum(test_m_groups*test_pi0s*train_ts_merged)/sum(test_m_groups*test_distribs))
  			}

  			print(fdr_estimate_fold1[k])

			train_idx <- test_idx
			test_idx  <- train_idx_list[[k]]

			ddhw_obj <- ddhw_grouped(unadj_p[train_idx], groups[train_idx], alpha, 
						regularization_term=reg_par, ...)

			train_rjs <- rejections(ddhw_obj, method="thresholds", groupwise=T)
			train_ts  <- thresholds(ddhw_obj)


			# merge adjacent groups which have the same weights
			merged_groups <- as.factor(groups)
			merged_levels <- cumsum(c(0,abs(diff(weights(ddhw_obj))) > .Machine$double.eps))
			merged_levels_rle <- rle(merged_levels)$lengths
			levels(merged_groups) <- merged_levels
			print(length(unique(levels(merged_groups))))

			cnt <- 1
			train_ts_merged <- c()
			train_rjs_merged <- c()
			for (i in seq_along(merged_levels_rle)){
				new_cnt <- cnt + merged_levels_rle[i]
				train_ts_merged[i] <- train_ts[cnt]
				train_rjs_merged[i] <- sum(train_rjs[cnt:(new_cnt-1)])
				cnt <- new_cnt
			}

	  		test_pv_list <- split(unadj_p[test_idx], merged_groups[test_idx])
  			test_m_groups <- sapply(test_pv_list,length)
  		
  			#test_rjs <- mapply(function(pv,t) sum(pv <= t), test_pv_list, train_ts)
  	
  			##test_fdrs <- mapply(my_grenander, test_pv_list, train_ts_merged)
			test_fdrs <- mapply(function(pv,t) my_grenander(pv,t,distrib=T), test_pv_list, train_ts_merged,SIMPLIFY=F)
			test_pi0s <- sapply(test_fdrs, function(x) x$pi_zero)
			test_distribs <- sapply(test_fdrs, function(x) x$distrib)
  			if (sum(train_rjs) == 0){
  				fdr_estimate_fold2[k] = 0
  			} else {
  				##fdr_estimate_fold2[k] <- sum(train_rjs_merged*test_fdrs)/sum(train_rjs_merged)
  				fdr_estimate_fold2[k] <- min(1,sum(test_m_groups*test_pi0s*train_ts_merged)/sum(test_m_groups*test_distribs))

  			}

  			print(fdr_estimate_fold2[k])

  		}
  		
  		fdr_estimate <- c(fdr_estimate_fold1,fdr_estimate_fold2)
  		print(paste("mean", mean(fdr_estimate,na.rm=T)))
  		print(mean(fdr_estimate,na.rm=T) + sd_penalty*sd(fdr_estimate,na.rm=T)/sqrt(nreps) )
  		# hacky way to deal with NA's that ocassionally pop up :( -> need to check this..(thought it was fixed)
  		decision_bool <- mean(fdr_estimate,na.rm=T) + sd_penalty*sd(fdr_estimate,na.rm=T)/sqrt(nreps) <= pi0*alpha

  		if (!is.na(decision_bool) && class(decision_bool) == "logical" && decision_bool){
 			break
 		} else {
 			reg_par <- 0
 		}
	}

    reg_path_df <- data.frame()
	
	#return(list(FDR=FDR,FDR_sd=FDR_sd))
	print(reg_par)

	if (is.factor(filterstat)){
		ddhw_obj <- ddhw_grouped(unadj_p, filterstat, alpha, regularization_term=reg_par, ...)

	} else {
		ddhw_obj <- ddhw(unadj_p, filterstat, nbins, alpha, regularization_term=reg_par, ...)
	}

	reg_path_information(ddhw_obj) <- reg_path_df
	ddhw_obj
	# return final result with 0.1..
}


ddhw_auto3 <- function(unadj_p, filterstat, alpha, 
		nbins="auto", nbins_check="auto", reg_par_number=7, reg_pars=NULL, nreps=10, sd_penalty=2, ...){

	pi0 <- fdrtool(unadj_p, statistic="pvalue", plot=F, verbose=F)$param[3]
	# some heuristic rule for determining number of bins
	if (nbins=="auto") nbins <- max(min(20, floor(length(unadj_p)/2000)), 1)
	if (nbins_check=="auto") nbins_check <- max(min(20, floor(length(unadj_p)/2000)), 1)

	# start by optimizing with what amounts to reg par = Infty
	ddhw_full <- ddhw(unadj_p, filterstat, nbins_check, alpha, ...)
	groups <- groups_by_filter(filterstat, nbins)
	# determine maximum necessary regularization
	max_reg <- 4*total_variation(weights(ddhw_full))

	print(paste("max reg:",max_reg))
	# for reg_par = 0 to max_reg     (e.g. 10 values?) do the holdout-experiment

	if (is.null(reg_pars)){
		reg_pars <- seq(0, max_reg, length=reg_par_number)
	}

	print(max_reg)
	#reg_pars <- c(0, rev(sapply(1:reg_par_number, function(x) max_reg/2^(x-1))))

	m <- length(unadj_p)
 
	# instead be less greedy: use regularization_term from above example and try it with that
	# if it does not work, revert to lambda = 0


	# generate indices for train/test datasets
	train_idx_list <- replicate(nreps, sample(1:m, floor(m/2)), simplify=F)

    for (j in rev(seq_along(reg_pars))){

    	fdr_estimate_fold1 <- rep(NA, nreps)
    	fdr_estimate_fold2 <- rep(NA, nreps)

    	reg_par <- reg_pars[j]
    	print(paste("par",reg_par))


    	for (k in 1:nreps){

			train_idx <- train_idx_list[[k]]
			test_idx  <- setdiff(1:m, train_idx)

			ddhw_obj <- ddhw(unadj_p[train_idx], filterstat[train_idx], nbins, alpha, 
						regularization_term=reg_par, ...)

			train_rjs <- rejections(ddhw_obj, method="thresholds", groupwise=T)
			train_ts  <- thresholds(ddhw_obj, levels_only=F)
			train_ts_fun <- approxfun(rank(filterstat[train_idx]), train_ts, method="linear",rule=2)


  			test_pv_list <- split(unadj_p[test_idx], groups_by_filter(filterstat[test_idx],nbins_check))
  			test_filter_list <- split(rank(filterstat[test_idx]), groups_by_filter(filterstat[test_idx],nbins_check))
  			test_m_groups <- sapply(test_pv_list,length)
  		
  			#test_rjs <- mapply(function(pv,t) sum(pv <= t), test_pv_list, train_ts)
  			##test_fdrs <- mapply(my_grenander, test_pv_list, train_ts_merged)
			test_fdrs <- mapply(function(pv,t) my_grenander(pv,t,distrib=T), test_pv_list, 1,SIMPLIFY=F)
			test_pi0s <- sapply(test_fdrs, function(x) x$pi_zero)
			test_distrib_funs <- sapply(test_fdrs, function(x) x$distrib_fun)

			summed_ts <- sapply(test_filter_list, function(fs) sum(train_ts_fun(fs)))
			summed_Fts <- mapply(function(fs,fun) sum(fun(train_ts_fun(fs))), test_filter_list, test_distrib_funs)

  			if (sum(train_rjs) == 0){
  				fdr_estimate_fold1[k] = 0
  			} else {
  				fdr_estimate_fold1[k] <- min(1,sum(test_pi0s*summed_ts)/sum(summed_Fts))
  			}



  			print(fdr_estimate_fold1[k])

			train_idx <- test_idx
			test_idx  <- train_idx_list[[k]]

			ddhw_obj <- ddhw(unadj_p[train_idx], filterstat[train_idx], nbins,alpha, 
						regularization_term=reg_par, ...)

			
			train_rjs <- rejections(ddhw_obj, method="thresholds", groupwise=T)
			train_ts  <- thresholds(ddhw_obj, levels_only=F)
			train_ts_fun <- approxfun(rank(filterstat[train_idx]), train_ts, method="linear",rule=2)


  			test_pv_list <- split(unadj_p[test_idx], groups_by_filter(filterstat[test_idx],nbins_check))
  			test_filter_list <- split(rank(filterstat[test_idx]), groups_by_filter(filterstat[test_idx],nbins_check))
  			test_m_groups <- sapply(test_pv_list,length)
  		
  			#test_rjs <- mapply(function(pv,t) sum(pv <= t), test_pv_list, train_ts)
  			##test_fdrs <- mapply(my_grenander, test_pv_list, train_ts_merged)
			test_fdrs <- mapply(function(pv,t) my_grenander(pv,t,distrib=T), test_pv_list, 1,SIMPLIFY=F)
			test_pi0s <- sapply(test_fdrs, function(x) x$pi_zero)
			test_distrib_funs <- sapply(test_fdrs, function(x) x$distrib_fun)

			summed_ts <- sapply(test_filter_list, function(fs) sum(train_ts_fun(fs)))
			summed_Fts <- mapply(function(fs,fun) sum(fun(train_ts_fun(fs))), test_filter_list, test_distrib_funs)

  			if (sum(train_rjs) == 0){
  				fdr_estimate_fold2[k] = 0
  			} else {
  				fdr_estimate_fold2[k] <- min(1,sum(test_pi0s*summed_ts)/sum(summed_Fts))
  			}

  			print(fdr_estimate_fold2[k])

  		}
  		
  		fdr_estimate <- (fdr_estimate_fold1 + fdr_estimate_fold2)/2
  		print(paste("mean", mean(fdr_estimate,na.rm=T)))
  		print(mean(fdr_estimate,na.rm=T) + sd_penalty*sd(fdr_estimate,na.rm=T)/sqrt(nreps) )
  		# hacky way to deal with NA's that ocassionally pop up :( -> need to check this..(thought it was fixed)
  		decision_bool <- mean(fdr_estimate,na.rm=T) + sd_penalty*sd(fdr_estimate,na.rm=T)/sqrt(nreps) <= pi0

  		if (!is.na(decision_bool) && class(decision_bool) == "logical" && decision_bool){
 			break
 		} else {
 			reg_par <- 0
 		}
	}

    reg_path_df <- data.frame()
	
	#return(list(FDR=FDR,FDR_sd=FDR_sd))
	print(reg_par)

	ddhw_obj <- ddhw(unadj_p, filterstat, nbins, alpha, regularization_term=reg_par, ...)
	reg_path_information(ddhw_obj) <- reg_path_df
	ddhw_obj
	# return final result with 0.1..
}



ddhw_unstratified <- function(unadj_p, filter_statistic, lambda, alpha, optim_pval_threshold=0.1){
	mfull <- length(unadj_p)


	o <- order(filter_statistic)
    ro <- order(o)


    pvals <- unadj_p[o]
    filterstats <- filter_statistic[o]

    #idx <- which(pvals > optim_pval_threshold)
    #rj_idx <- which(pvals <= alpha/mfull) #make sure bonferroni rejected p-values get rejected globally
  


    # start constructing LP

    rle_enc <- rle(ifelse(pvals > optim_pval_threshold, -1, pvals))
    rle_idx <- cumsum(rle_enc$lengths)
    m <- length(rle_enc$values)
    ignored_idx <- which(rle_enc$values == -1) #pvals which will not get rejected anyway
    rj_idx      <- which(rle_enc$values <= alpha/mfull) #pvals which will get rejected for sure
    # first only the main problem
    row_idx <- c(1:m, 1:m)
    col_idx <- c(1:(2*m))
    #coef_val <- c(-pvals, rep(1,m))
    coef_val <- c(-rle_enc$values, rep(1,m))

    #print(paste(length(row_idx), length(col_idx), length(coef_val)))
	coef_mat <- slam::simple_triplet_matrix(row_idx,col_idx,coef_val)
	#fdr_constr <- matrix(c(rep(alpha,m), rep(-1,m)),nrow=1)
	fdr_constr <- matrix(c(rep(alpha,m), -1*rle_enc$lengths),nrow=1)

	coef_mat <- rbind(coef_mat,fdr_constr)


	 # Ok now we also unfortunately need the regularization, i.e. introduce m-1 new terms
	coef_mat <- cbind(coef_mat,
	 				slam::simple_triplet_zero_matrix(m+1, m-1 ,mode="double"))


	 				# constraints f_i >= t(i+1) - t_i
	 				# first set: t_i
	 				# second set: t_(i+1)
	 				# third set: f_i
	new_row_idx <- c(1:(m-1), 1:(m-1), 1:(m-1))
	new_col_idx <- c((m+1):(2*m-1), (m+2):(2*m), (2*m+1):(3*m-1)) 
	new_coef_val <- c(rep(-1, m-1), rep(1,m-1), rep(1,m-1))

	#print("2")
	coef_half <- slam::simple_triplet_matrix(new_row_idx, new_col_idx, new_coef_val)
 
	new_coef_val2 <- c(rep(1, m-1), rep(-1,m-1), rep(1,m-1))
	coef_half2 <- slam::simple_triplet_matrix(new_row_idx, new_col_idx, new_coef_val2)


	final_constr <- matrix(c(rep(0,m), rep(lambda/m, m), rep(-1,m-1)),nrow=1)

	coef_mat <- rbind(coef_mat, coef_half, coef_half2, final_constr)

	obj <- c(rep(1,m), rep(0, 2*m-1))
	model_ub <- c(rep(1, 2*m), rep(2,m-1))
	#model_ub[idx] <- 0
	model_ub[ignored_idx] <- 0

	model_lb <- rep(0, 3*m-1)
	#model_lb[rj_idx] <- 1
	#require("gurobi")
	#require("Matrix")

	# convert simple_triplet_matrix of pkg Slam to sparse matrix of pkg Matrix
	# gurobi has problems with simple_triplet_matrix for an (undocumented) reason (?bug)
	#coef_mat <-  sparseMatrix(i=coef_mat$i, j=coef_mat$j, x=coef_mat$v,
     #      dims=c(coef_mat$nrow, coef_mat$ncol))

	#model <- list()
	#model$A <- coef_mat
	#model$obj <- obj
	#model$modelsense <- "max"
	#model$rhs        <- 0
	#model$lb         <- model_lb
	#model$ub         <- model_ub
	#model$sense      <- '>'
	#model$vtype      <- model_vtype


		#params <- list(ScaleFlag=0 ,
		#					OutputFlag = 0)#, TimeLimit=10)
		
		#res <- gurobi(model, params)
		#print( str(res ))

		rsymphony_ub <- list(upper = list(ind = 1:ncol(coef_mat), val = model_ub))

		#if (is.infinite(time_limit)) time_limit <- -1
		#if (is.infinite(node_limit)) node_limit <- -1
		res<- Rsymphony::Rsymphony_solve_LP(obj, coef_mat, rep(">=", nrow(coef_mat)), 
				rep(0, nrow(coef_mat)),
    			bounds= rsymphony_ub,
				max = T, verbosity = 2)
		sol <- res$solution
		#solver_status <- res$status

		#ts <- res$x[(m+1):(2*m)]
		ts <- sol[(m+1):(2*m)]
		rle_enc$values <- ts
		ts <- inverse.rle(rle_enc)
		ws <- thresholds_to_weights_full(ts)

		ts <- ts[ro]
		ws <- ws[ro]

		# linear interpolation
		# outside observed range use same values

		weighted_pvalues <- mydiv(unadj_p, ws)
		ts_corrected <- ws*get_bh_threshold(weighted_pvalues, alpha)

		ts_fun <- approxfun(filterstats, ts_corrected, method="linear",rule=2)

		adj_pv <- p.adjust(weighted_pvalues,method="BH")
		return(list(lambda=lambda, res=res,ts=ts,ws=ws, ts_corrected=ts_corrected,
			filterstats=filterstats, rjs=sum(adj_pv <= alpha), ts_fun=ts_fun))
		#sol <- res$x		#solver_status <- res$status

}

ddhw_unstratified_path <- function(unadj_p, filter_statistic, reg_pars, alpha, nreps=10, ...){
# generate indices for train/test datasets

	m <- length(unadj_p)
	nbins <- floor(m/1000)
	train_idx_list <- replicate(nreps, sample(1:m, floor(m/2)), simplify=F)

    for (j in rev(seq_along(reg_pars))){

    	fdr_estimate_fold1 <- rep(NA, nreps)
    	fdr_estimate_fold2 <- rep(NA, nreps)

    	reg_par <- reg_pars[j]
    	print(paste("par",reg_par))


    	for (k in 1:nreps){

			train_idx <- train_idx_list[[k]]
			test_idx  <- setdiff(1:m, train_idx)

			# probably could pre-sort to save some time
			ddhw_obj <- ddhw_unstratified(unadj_p[train_idx], filter_statistic[train_idx], reg_par, alpha, ...)

			train_rjs <- ddhw_obj$rjs
			train_ts_fun <- ddhw_obj$ts_fun


  			test_pv_list <- split(unadj_p[test_idx], groups_by_filter(filter_statistic[test_idx],nbins))
  			test_filter_list <- split(filter_statistic[test_idx], groups_by_filter(filter_statistic[test_idx],nbins))
  			test_m_groups <- sapply(test_pv_list,length)
  		
  			#test_rjs <- mapply(function(pv,t) sum(pv <= t), test_pv_list, train_ts)
  			##test_fdrs <- mapply(my_grenander, test_pv_list, train_ts_merged)
			test_fdrs <- mapply(function(pv,t) my_grenander(pv,t,distrib=T), test_pv_list, 1,SIMPLIFY=F)
			test_pi0s <- sapply(test_fdrs, function(x) x$pi_zero)
			test_distrib_funs <- sapply(test_fdrs, function(x) x$distrib_fun)

			summed_ts <- sapply(test_filter_list, function(fs) sum(train_ts_fun(fs)))
			summed_Fts <- mapply(function(fs,fun) sum(fun(train_ts_fun(fs))), test_filter_list, test_distrib_funs)

  			if (sum(train_rjs) == 0){
  				fdr_estimate_fold1[k] = 0
  			} else {
  				fdr_estimate_fold1[k] <- min(1,sum(test_pi0s*summed_ts)/sum(summed_Fts))
  			}


  			#test_ts <- train_ts_fun(filter_statistic[test_idx])
  			#test_rjs <- sum(unadj_p[test_idx] <= test_ts)
  	
  			#print(paste("train_rjs:",train_rjs,"   test_rjs:",test_rjs))
  			
  			#if (train_rjs == 0){
  			#	fdr_estimate_fold1[k] <- 0
  			#} else {
  			#	fdr_estimate_fold1[k] <- min(1, sum(test_ts)/test_rjs)
  			#}

  			print(fdr_estimate_fold1[k])

			train_idx <- test_idx
			test_idx  <- train_idx_list[[k]]

			# probably could pre-sort to save some time
			ddhw_obj <- ddhw_unstratified(unadj_p[train_idx], filter_statistic[train_idx], reg_par, alpha, ...)

			train_rjs <- ddhw_obj$rjs
			train_ts_fun <- ddhw_obj$ts_fun


  			#test_ts <- train_ts_fun(filter_statistic[test_idx])
  			#test_rjs <- sum(unadj_p[test_idx] <= test_ts)
  	

  			
  			#if (train_rjs == 0){
  			#	fdr_estimate_fold2[k] <- 0
  			#} else {
  			#	fdr_estimate_fold2[k] <- min(1, sum(test_ts)/test_rjs)
  			#}


  			test_pv_list <- split(unadj_p[test_idx], groups_by_filter(filter_statistic[test_idx],nbins))
  			test_filter_list <- split(filter_statistic[test_idx], groups_by_filter(filter_statistic[test_idx],nbins))
  			test_m_groups <- sapply(test_pv_list,length)
  		
  			#test_rjs <- mapply(function(pv,t) sum(pv <= t), test_pv_list, train_ts)
  	
  			##test_fdrs <- mapply(my_grenander, test_pv_list, train_ts_merged)
			test_fdrs <- mapply(function(pv,t) my_grenander(pv,t,distrib=T), test_pv_list, 1,SIMPLIFY=F)
			test_pi0s <- sapply(test_fdrs, function(x) x$pi_zero)
			test_distrib_funs <- sapply(test_fdrs, function(x) x$distrib_fun)

			summed_ts <- sapply(test_filter_list, function(fs) sum(train_ts_fun(fs)))
			summed_Fts <- mapply(function(fs,fun) sum(fun(train_ts_fun(fs))), test_filter_list, test_distrib_funs)
  			if (sum(train_rjs) == 0){
  				fdr_estimate_fold2[k] = 0
  			} else {
  				##fdr_estimate_fold2[k] <- sum(train_rjs_merged*test_fdrs)/sum(train_rjs_merged)
  				fdr_estimate_fold2[k] <- min(1,sum(test_pi0s*summed_ts)/sum(summed_Fts))

  			}


  			print(fdr_estimate_fold2[k])

  		}
  		
  		fdr_estimate <- (fdr_estimate_fold1 + fdr_estimate_fold2)/2
  		print(paste("mean", mean(fdr_estimate,na.rm=T)))
  		#print(mean(fdr_estimate,na.rm=T) + sd_penalty*sd(fdr_estimate,na.rm=T)/sqrt(nreps) )
  		# hacky way to deal with NA's that ocassionally pop up :( -> need to check this..(thought it was fixed)
  		decision_bool <- mean(fdr_estimate,na.rm=T)  <= alpha#+ sd_penalty*sd(fdr_estimate,na.rm=T)/sqrt(nreps) <= alpha

  		if (!is.na(decision_bool) && class(decision_bool) == "logical" && decision_bool){
 			break
 		} else {
 			reg_par <- 0
 		}
	}

    #reg_path_df <- data.frame()
	
	#return(list(FDR=FDR,FDR_sd=FDR_sd))
	print(reg_par)

	ddhw_unstratified(unadj_p, filter_statistic, reg_par, alpha, ...)

	#ddhw_obj <- ddhw(unadj_p, filterstat, nbins, alpha, regularization_term=reg_par, ...)
	#reg_path_information(ddhw_obj) <- reg_path_df
	#ddhw_obj
	# return final result with 0.1..
}
