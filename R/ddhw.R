


#' ddhw: Data-Driven Hypothesis Weights
#'
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param groups  Vector which assigns a group to each p-value.
#' @param filter_statistic  Vector which contains the filter-statistic (covariate, independent under the H0) for each p-value
#' @param nbins  Integer, number of groups into which p-values will be split based on filter_statistic.
#' @param alpha   Numeric, sets the nominal level for FDR control.
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
		lp_relaxation=TRUE,
		pi0_adaptive=FALSE) # only relevant if optim_method == "MILP"
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

		if (pi0_adaptive){
			pvals_list <- split(unadj_p, groups)
			pi0_groups <- sapply(pvals_list, lsl_pi0_est)
			pi0 <- 0.9*mean(pi0)+0.1
			ws <- 1/pi0
		} else {
			ws <- 1
		}

		weighted_pvalues <- unadj_p / ws 
		ts <- get_bh_threshold(weighted_pvalues, alpha, mtests=mtests)

		solver_information <- list("R p.adjust function")


		if (nbins == 1){
			message("Only 1 group supplied, this reduces to the classic Benjamini-Hochberg method.")
		} else if (regularization_term==0){
			message("Regularization term equal to 0, this reduces to the classic Benjamini-Hochberg method.")
			ws <- rep(ws, nbins)
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
					lp_relaxation=lp_relaxation,
					pi0_adaptive=pi0_adaptive)

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
			bh_rejections=bh_rejections, rejections_bound=rejections_bound, pi0_adaptive=pi0_adaptive)
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
					 	lp_relaxation=lp_relaxation,
					 	pi0_adaptive = F){

	m <- length(unadj_p)
	nbins <- length(levels(groups))

	pvals_list <- split(unadj_p, groups)
	m_groups <- sapply(pvals_list,length) # number of pvals in each group (BEFORE applying optim_pval_threshold)

	if (pi0_adaptive){
		#pi0_groups <- sapply(pvals_list, function(pv) fdrtool(pv, statistic="pvalue", cutoff.method="pct0", pct0=0.5, plot=F, verbose=F)$param[3])
		pi0_groups <- sapply(pvals_list, lsl_pi0_est)
		pi0_groups <- pmin(0.9*pi0_groups + 0.1,1)
		#pi0_groups <- pmin((1-alpha)*pi0_groups + alpha,1)

		#pi0_groups <- 0.99*pi0_groups + 0.01 # add conservative bias
	}

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

	if (pi0_adaptive){
		diff_coefs <- unlist(mapply(function(p, m_group) diff(c(0,p))* m_group , pvals_list, m_groups*pi0_groups))
	}
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

		if (pi0_adaptive){
			plugin_fdr_constraint2 <- matrix(c(rep(alpha/m/pmax,z_vars), -pi0_groups*m_groups/m), nrow=1)
		}
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

 			if (pi0_adaptive){
 				regularization_constraint <- matrix(c(rep(0,z_vars), regularization_term/m*m_groups*pi0_groups, rep(-1, nbins-1)), nrow=1)
			}
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
		res<- lpsymphony::lpsymphony_solve_LP(obj, full_triplet, rep(">=", nrows), 
			model_rhs,
    		types = model_vtype, bounds= rsymphony_ub,
			max = T, verbosity = -2, time_limit = time_limit,
			node_limit = node_limit, gap_limit = ifelse(mip_gap <= 10^(-4), -1, mip_gap*100), #now default mip_gap=10^(-4) as in Gurobi while -1 default for Rsymphony
			first_feasible = FALSE)
		sol <- res$solution
		#print(str(res))
		solver_status <- res$status
		
	} else if (solver=="gurobi"){

		# keep code for commercial solver for now
		#speed up compared to Symphony appears to range about factor ~2-10, depending on problem (even more for hard problems)
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
		# turns out that enforcing T_g >= t_g is numerically hard in the context of this problem

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
	####################################################################################################

	# further note to self: What happens when we solve the linear relaxation of the unregularized procedure

	ws <-  if (all(t_thresholds1 == .0) || all(t_thresholds == .0)){
				rep(1,nbins)
			} else {
				t_thresholds*m/sum(m_groups*t_thresholds) #might need adjustment if mtests > m
			}

	if (pi0_adaptive){
		ws <-  if (all(t_thresholds1 == .0) || all(t_thresholds == .0)){
				rep(1,nbins)
			} else {
				t_thresholds*m/sum(m_groups*pi0_groups*t_thresholds) #might need adjustment if mtests > m
			}
	}

	t_thresholds <- t_thresholds*pmax

	list(ts = t_thresholds, ws = ws, solver_information = list(solver=solver, solver_status = solver_status))
}



#' ddhw_auto: Main function to do Data-Driven Hypothesis Weights
#'
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param filterstat Vector which contains the filter-statistic (covariate, independent under the H0) for each p-value. This can be a numeric vector or a factor.
#' @param alpha   Numeric, sets the nominal level for FDR control.
#' @param nbins  Integer, number of groups into which p-values will be split based on filter_statistic (default = "auto").
#' @param reg_par_number	Integer, number of regularization parameters to try out (length of grid)
#' @param reg_pars  "linear" if sequence of regularization parameters should form a grid, "power" if it should be geometrically decreasing, a numeric vector if the user wants to explicitly prespecify the regularization parameter values to test
#' @param test_distrib  Method to estimate the distribution function in the test set (after holdout split). Available "grenander" (Strimmer's modified Grenander estimator in fdrtools), "ecdf" (the ECDF), "both" (the maximum of "ecdf" and "grenander")
#' @param nreps         Number of replications of the 2-fold CV procedure (i.e. 2*nreps total holdout runs)
#' @param sd_penalty    Additional penalization before accepting a regularization parameter value
#' @param verbosity_level  Verbosity level
#' @return  A ddhw object
#'
#' @seealso the parameters passed to ddhw / ddhw_grouped \code{\link{ddhw}}

ddhw_auto <- function(unadj_p, filterstat, alpha,
		nbins="auto", reg_par_number=7, reg_pars="linear", test_distrib="grenander",  nreps=10, sd_penalty=0,
		verbosity_level = 0, ...){

	pi0 <- 1

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

	if (verbosity_level > 0) print(paste("max reg:",max_reg))

	# for reg_par = 0 to max_reg     (e.g. 10 values?) do the holdout-experiment

	if (reg_pars == "linear"){
		reg_pars <- seq(0, max_reg, length=reg_par_number)
	} else if (reg_pars == "power"){
		reg_pars <- max_reg*0.5^(0:reg_par_number)
	}


	m <- length(unadj_p)
 

	# generate indices for train/test datasets
	train_idx_list <- replicate(nreps, sample(1:m, floor(m/2)), simplify=F)

    for (j in rev(seq_along(reg_pars))){

    	fdr_estimate_fold1 <- rep(NA, nreps)
    	fdr_estimate_fold2 <- rep(NA, nreps)

    	reg_par <- reg_pars[j]
    	if (verbosity_level > 0) print(paste("par",reg_par))


    	for (k in 1:nreps){

			train_idx <- train_idx_list[[k]]
			test_idx  <- setdiff(1:m, train_idx)

			ddhw_obj <- ddhw_grouped(unadj_p[train_idx], groups[train_idx], alpha, 
						regularization_term=reg_par, ...)

			train_rjs <- rejections(ddhw_obj, method="thresholds", groupwise=T)
			train_ts  <- thresholds(ddhw_obj)

			train_ts_merged <- train_ts
			train_rjs_merged <- train_rjs
			merged_groups <- groups

	  		test_pv_list <- split(unadj_p[test_idx], merged_groups[test_idx])
  			test_m_groups <- sapply(test_pv_list,length)


			test_fdrs <- mapply(function(pv,t) my_grenander(pv,t,distrib=T), test_pv_list, train_ts_merged, SIMPLIFY=F)
			test_pi0s <- sapply(test_fdrs, function(x) x$pi_zero)
			test_distribs <- sapply(test_fdrs, function(x) x$distrib)
  			
  			if (test_distrib == "ecdf" || test_distrib == "both"){
  				test_distribs_ecdf <- mapply(function(pv,t,mg) sum(pv <= t)/mg, test_pv_list, train_ts_merged, test_m_groups, SIMPLIFY=T)
  				if (test_distrib == "both"){
  					test_distribs <- pmax(test_distribs_ecdf, test_distribs)
  				} else {
  					test_distribs <- test_distribs_ecdf
  				}
  			}

  			
  			if (sum(train_rjs) == 0){
  				fdr_estimate_fold1[k] <- 1
  			} else {
  				fdr_estimate_fold1[k] <- min(1,sum(test_m_groups*test_pi0s*train_ts_merged)/sum(test_m_groups*test_distribs))
  			}

  			if (verbosity_level > 0) print(fdr_estimate_fold1[k])

			train_idx <- test_idx
			test_idx  <- train_idx_list[[k]]

			ddhw_obj <- ddhw_grouped(unadj_p[train_idx], groups[train_idx], alpha, 
						regularization_term=reg_par, ...)

			train_rjs <- rejections(ddhw_obj, method="thresholds", groupwise=T)
			train_ts  <- thresholds(ddhw_obj)


			train_ts_merged <- train_ts
			train_rjs_merged <- train_rjs
			merged_groups <- groups

	  		test_pv_list <- split(unadj_p[test_idx], merged_groups[test_idx])
  			test_m_groups <- sapply(test_pv_list,length)
  		
  	
			test_fdrs <- mapply(function(pv,t) my_grenander(pv,t,distrib=T), test_pv_list, train_ts_merged,SIMPLIFY=F)
			test_pi0s <- sapply(test_fdrs, function(x) x$pi_zero)
			test_distribs <- sapply(test_fdrs, function(x) x$distrib)


  			if (test_distrib == "ecdf" || test_distrib == "both"){
  				test_distribs_ecdf <- mapply(function(pv,t,mg) sum(pv <= t)/mg, test_pv_list, train_ts_merged, test_m_groups, SIMPLIFY=T)
  				if (test_distrib == "both"){
  					test_distribs <- pmax(test_distribs_ecdf, test_distribs)
  				} else {
  					test_distribs <- test_distribs_ecdf
  				}

  			}

  			if (sum(train_rjs) == 0){
  				fdr_estimate_fold2[k] = 1
  			} else {
  				fdr_estimate_fold2[k] <- min(1,sum(test_m_groups*test_pi0s*train_ts_merged)/sum(test_m_groups*test_distribs))
  			}

  			if (verbosity_level > 0) print(fdr_estimate_fold2[k])

  		}
  		
  		fdr_estimate <-  (fdr_estimate_fold1 + fdr_estimate_fold2)/2

  		if (verbosity_level > 0) print(mean(fdr_estimate,na.rm=T) + sd_penalty*sd(fdr_estimate,na.rm=T))
  		# hacky way to deal with NA's that ocassionally pop up 
  		decision_bool <- mean(fdr_estimate,na.rm=T) + sd_penalty*sd(fdr_estimate,na.rm=T)/sqrt(nreps) <= pi0*alpha

  		#fdr_estimate <- c(fdr_estimate_fold1,fdr_estimate_fold2)
  		#print(paste("mean", mean(fdr_estimate,na.rm=T)))
  		#print(mean(fdr_estimate,na.rm=T) + sd_penalty*sd(fdr_estimate,na.rm=T)/sqrt(nreps) )
  		# hacky way to deal with NA's that ocassionally pop up :( -> need to check this..(thought it was fixed)
  		#decision_bool <- mean(fdr_estimate,na.rm=T) + sd_penalty*sd(fdr_estimate,na.rm=T)/sqrt(nreps) <= pi0*alpha

  		if (!is.na(decision_bool) && class(decision_bool) == "logical" && decision_bool){
 			break
 		} else {
 			reg_par <- 0
 		}
	}

    reg_path_df <- data.frame()
	
	#return(list(FDR=FDR,FDR_sd=FDR_sd))
	if (verbosity_level > 0) print(reg_par)

	if (is.factor(filterstat)){
		ddhw_obj <- ddhw_grouped(unadj_p, filterstat, alpha, regularization_term=reg_par, ...)

	} else {
		ddhw_obj <- ddhw(unadj_p, filterstat, nbins, alpha, regularization_term=reg_par, ...)
	}

	reg_path_information(ddhw_obj) <- reg_path_df

	if (rejections(ddhw_obj) < sum(p.adjust(unadj_p, method="BH") <= alpha)){ # actually should be less than the adaptive version

				ddhw_obj <- ddhw_grouped(unadj_p, filterstat, alpha, regularization_term=0, ...)

	}

	ddhw_obj
}

