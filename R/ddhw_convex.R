# actually nevermind let's just work on it directly.

ddhw_tmp <- function(pvalues, filter_statistics, alpha,
						filter_statistic_type = "ordered",
						nbins = 10,
						nfolds = 10,
						lambda = Inf,
						lp_solver="gurobi"){

	# code from R's p.adjust function to gracefully handle NAs
	nm <- names(pvalues)
    pvalues <- as.numeric(pvalues)
    p0 <- setNames(pvalues, nm)
   if (all(nna <- !is.na(pvalues))) 
        nna <- TRUE

    pvalues <- pvalues[nna]
    filter_statistics <- filter_statistics[nna]
    weights <- rep(NA, length(pvalues))

	if (filter_statistic_type =="ordered"){
		groups <- groups_by_filter(filter_statistics, nbins)
	} else {
		groups <- as.factor(filter_statistics)
		nbins <- nlevels(groups)
	}

	# sort pvalues globally
	order_pvalues <- order(pvalues)

	sorted_groups <- groups[order_pvalues]
	sorted_pvalues <- pvalues[order_pvalues]

	split_sorted_pvalues <- split(sorted_pvalues, sorted_groups)
	m_groups <- sapply(split_sorted_pvalues, length)

	#recover weights

	# afterwards do the k-fold strategy
	set.seed(1)
	folds <- sample(1:nfolds, length(pvalues), replace = TRUE)
	sorted_folds <- folds[order_pvalues]

	for (i in 1:nfolds){
		# don't worry about inefficencies right now, they are only O(n) and sorting dominates this
		filtered_sorted_groups <- sorted_groups[sorted_folds!=i]
		filtered_sorted_pvalues <- sorted_pvalues[sorted_folds!=i]
		filtered_split_sorted_pvalues <- split(filtered_sorted_pvalues, filtered_sorted_groups)

		m_groups_holdout_fold <- m_groups - sapply(filtered_split_sorted_pvalues, length)

		res <- ddhw_convex(filtered_split_sorted_pvalues, alpha, m_groups_holdout_fold,
						   lambda=lambda, lp_solver=lp_solver)

		weights[folds == i] <- res$ws[groups[folds==i]]
	# I do need some kind of new object though...
	}

	weighted_pvalues <- mydiv(pvalues, weights)
	adj_p <- p.adjust( weighted_pvalues, method = "BH")
	return(list(adj_p=adj_p,ws=weights,filter_statistics=filter_statistics, groups=groups))
}


ddhw_convex <- function(split_sorted_pvalues, alpha, m_groups, lambda=Inf, lp_solver="gurobi"){

	nbins <- length(split_sorted_pvalues)
	m <- sum(m_groups)

	if (nbins != length(m_groups)){
		stop("length of m_groups should be equal to number of bins")
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

		total_variation_constr <- matrix(c(rep(0,nbins), -lambda*rep(1,nbins), nbins*rep(1,nbins-1)),nrow=1)
		constr_matrix <- rbind(constr_matrix, total_variation_constr)

		rhs <- c(rhs,rep(0, 2*(nbins-1)),0) #add RHS for absolute differences and for total variation penalty

	}

	# incorporate the FDR constraint
	fdr_constr<- matrix(c(rep(-alpha,nbins), rep(1,nbins), rep(0,ncol(constr_matrix)-2*nbins)), nrow=1)
	constr_matrix <- rbind(constr_matrix, fdr_constr)
	rhs <- c(rhs,0)

	model <- list()
	model$A <- sparseMatrix(i=constr_matrix$i, j=constr_matrix$j, x=constr_matrix$v,
           dims=c(constr_matrix$nrow, constr_matrix$ncol))
	model$obj <- obj
	model$modelsense <- "max"
	model$rhs        <- rhs
	model$lb         <- 0
	model$ub         <- 1
	model$sense      <- '<'
	#model$vtype      <- model_vtype

	params <- list()
	res <- gurobi(model, params)
	sol <- res$x
	solver_status <- res$status
	ts <- sol[(1:nbins)+nbins]
	ws <- thresholds_to_weights(ts, m_groups)	#ts/sum(ts)*nbins



	return(list(ws=ws))
}


presorted_grenander <- function(sorted_pvalues){
  	n  <- length(sorted_pvalues)
  	sorted_pvalues <- sorted_pvalues
  	unique_pvalues <- unique(sorted_pvalues)
  	ecdf_values <- cumsum(tabulate(match(sorted_pvalues, unique_pvalues)))/n
  	unique_pvalues <- c(0,unique_pvalues)
  	ecdf_values   <- c(0, ecdf_values)
  	ll <- gcmlcm(unique_pvalues, ecdf_values, type="lcm")
	ll$length <- length(ll$slope.knots)
	ll$x.knots <- ll$x.knots[-ll$length]
	ll$y.knots <- ll$y.knots[-ll$length]
	ll
}