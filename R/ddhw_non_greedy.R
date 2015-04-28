ddhw_grouped_nongreedy <- function(unadj_p, groups, rjs,
		mtests = length(unadj_p),
		solver = "Rsymphony", #only relevant if optim_method == "MILP"
		regularization_term = Inf,
		penalty="total variation",
		optim_pval_threshold = 1,
		time_limit=Inf,
		node_limit=Inf,
		threads=0,
		solution_limit=Inf,
		mip_gap_abs=Inf,
		mip_gap = 10^(-4),
		MIPFocus=0) # only relevant if optim_method == "MILP"
{

	groups <- as.factor(groups) #later: check if factor actually has correct levels..

	if (!check_uniform_groups(groups)){
		warning(paste0("Software and algorithm has only been tested with groups",
			" of equal size."))	
	}

	nbins <- length(levels(groups))

	


	res <- ddhw_grouped_milp_nongreedy(unadj_p, groups, rjs,
					mtests, optim_pval_threshold,
					solver=solver,
					time_limit = time_limit, node_limit=node_limit,
					solution_limit=solution_limit,
					mip_gap_abs=mip_gap_abs,
					mip_gap = mip_gap,
					threads=threads,
					MIPFocus=MIPFocus)

	solver_information <- res$solver_information
	ws  <- res$ws
	ts  <- res$ts

		# idea to explore later:
		# instead of using the threshold naturally defined by the procedure, rerun BH with weights returned.
		# This ensures that IF a suboptimal solution was found
		# 1) due to numerical instability
		# 2) because the user set node_limit or time_limit options
		# the BH procedure will be able to improve upon this solution!


	weighted_pvalues <- mydiv(unadj_p, ws[groups])

	adj_p <- p.adjust( weighted_pvalues, method = "BH", n = mtests)

	df <- data.frame(pvalue = unadj_p, adj_p = adj_p, group= groups)
	# maybe introduce more complicated output object later, but for now, just adjp should be fine!
	new("ddhw",
		 	df = df,
		 	weights = ws,
		 	thresholds = ts,
		 	alpha = rjs,
			nbins = as.integer(nbins),
			regularization_term = regularization_term,
			penalty = penalty,
		 	solver_information = solver_information)
}

ddhw_nongreedy <- function(unadj_p, filter_statistic, nbins, rjs,...){
	groups <- groups_by_filter(filter_statistic, nbins)
	ddhw_object <- ddhw_grouped_nongreedy(unadj_p, groups, rjs, ...)
	filter_statistics(ddhw_object) <- filter_statistic
 	ddhw_object
}

ddhw_grouped_milp_nongreedy <- function(unadj_p, groups, rjs, mtests,
					 	optim_pval_threshold, 
					 	solver="Rsymphony",
					 	time_limit=Inf, node_limit=Inf,
					 	threads=0,
					 	solution_limit=solution_limit,
					 	mip_gap_abs=mip_gap_abs,
					 	mip_gap = mip_gap,
					 	MIPFocus=MIPFocus){

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


	diff_coefs <- unlist(mapply(function(p, m_group) diff(c(0,p))* m_group , pvals_list, m_groups))

	z_vars <- ncol(full_triplet) #number of binary variables

	full_triplet <- rbind(full_triplet, matrix(c(1,-1), nrow=2, ncol=z_vars))

	obj        <- -diff_coefs # maximize number of rejections


    nrows <- nrow(full_triplet)
    model_rhs <- c(rep(0,nrows-2), rjs,-rjs)
    model_vtype <- c(rep('B', z_vars))

    # model ub will actually depend on which penalty we choose
    # model ub = 1 for binary variables, 1 for thresholds, x for f_g, where
    # x=2 for total variation, x=2*m for uniform deviation
    model_ub <- c(rep(1, z_vars))
    
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
							OutputFlag = 0, 
							Threads=threads,MIPFocus=MIPFocus,
							MIPGap = mip_gap,  ResultFile='buggy_model.mps')

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
	t_thresholds <- t_thresholds1
				

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

	list(ts = t_thresholds1, ws = ws, solver_information = list(solver=solver, solver_status = solver_status))
}