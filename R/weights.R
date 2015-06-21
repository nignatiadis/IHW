# functions to assist in working with weights


total_variation <- function(ws){
    sum(abs(diff(ws)))
}

uniform_deviation <- function(ws){
    sum(abs(ws-1))
}

thresholds_to_weights <- function(ts, m_groups){

    if (length(ts) != length(m_groups)) stop("incosistent number of elements for ts and m_groups")

    nbins <- length(ts)
    m <- sum(m_groups)

    if (all(ts == .0) ){
        rep(1,nbins)
    } else {
        ts*m/sum(m_groups*ts) 
    }
}

thresholds_to_weights_full <- function(ts){

    m <- length(ts)

    if (all(ts == .0) ){
        rep(1,m)
    } else {
        ts*m/sum(ts) 
    }
}
#' Fix numerically unstable MILP thresholds (Helper function)
#'
#' second small linear program to be solved if MILP was numerically unstable
#' then we try to manually check if total variation regularization can really be 
#'  fullfilled
#'
#' @param ts Thresholds returned from MILP (corresponding to highest rejected p-value in each group)
#' @param m_groups  Number of tests in each stratum
#' @param rjs Number of rejections at MILP solver optimum
#' @param alpha Nominal level of multiple testing procedure
#' @param lambda Regularization parameter (so that total_variation(weights) <= lambda)
#' @return Thresholds >= ts, which correspond to weights obeying the regularization condition

regularize_thresholds <- function(ts, m_groups, rjs, alpha, lambda, 
                    penalty="total variation", solver="Rsymphony"){

    if (length(ts) != length(m_groups)) stop("incosistent number of elements for ts and m_groups")

    nbins <- length(ts)
    m <- sum(m_groups)

    # solve the LP
    # Ts >= ts
    # s.t. sum(Ts_reg*m_groups) <= rjs*alpha
    # s.t. sum|m*Ts_{i+1} - m*Ts_reg_{i}| <= lambda * sum(m_g * Ts_reg_{g})


    if (penalty=="total variation"){

            # |T_g - T_{g-1}| = f_g  introduce absolute value constants, i.e. introduce nbins-1 new continuous variables


            # we now need inequalities: T_g - T_{g-1} + f_g >= 0    and -T_g + T_{g-1} + f_g >= 0
            # start with matrix diff_matrix (nbins-1) X nbins as a scaffold i.e.
            #
            #  -1  1
            #     -1  1
            #        -1  1


            diff_matrix <- diag(-1,nbins-1,nbins) + cbind(rep(0,nbins-1), diag(1,nbins-1))

            constraint_matrix <- rbind(cbind(diff_matrix, diag(1, nrow=nbins-1)),
                                       cbind(-diff_matrix, diag(1, nrow=nbins-1)))

            plugin_fdr_constraint <- matrix(c(-m_groups, rep(0, nbins-1)), nrow=1) # >= - rjs*alpha


            # add final regularization inequality:
            #           \sum |w_g - w_{g-1}| <= reg_par
            # <=>       \sum |T_g - T_{g-1}|*m <= reg_par * \sum m_g T_g
            regularization_constraint <- matrix(c(lambda*m_groups, m*rep(-1, nbins-1)), nrow=1)

            constraint_matrix <- rbind(plugin_fdr_constraint, regularization_constraint, constraint_matrix)

            # to be used for ub afterwards
            model_ub  <- c(rep(1, nbins), rep(2, nbins-1))
            model_lb  <- c(ts,            rep(0, nbins-1))
            model_rhs <- c(-rjs*alpha, rep(0,nrow(constraint_matrix)-1))

            #model_obj <- c(m_groups, rep(0, nbins-1)) #instead maybe try to minimize the deviations.
            #model_obj <- c(rep(0, nbins), rep(-1,nbins-1))
            model_obj <- c(-m_groups, rep(0, nbins-1))
    } else if (penalty=="uniform deviation"){
        stop("not yet implemented")
    }


    model_vtype <- rep('C', ncol(constraint_matrix))

   
    if (solver=="Rsymphony") {

        rsymphony_bounds <- list(upper = list(ind = 1:ncol(constraint_matrix), val = model_ub),
                                 lower = list(ind = 1:ncol(constraint_matrix), val = model_lb))

     
        res<- Rsymphony::Rsymphony_solve_LP(model_obj, constraint_matrix, rep(">=", nrow(constraint_matrix)),
                                                model_rhs,  types = model_vtype, bounds= rsymphony_bounds,
                                                max = F, first_feasible = FALSE)
        sol <- res$solution
        print(str(res))
        solver_status <- res$status
    }
    sol[1:nbins]
}

# normalize weights so that their sum will be equal to length(ws)
# also make sure no negative weights appear
normalize_weights <- function(ws){
	abs(ws)*length(ws)/sum(abs(ws))
}

regularize_weights <- function(ws, lambda){
	if (lambda <0 || lambda>1) stop("Regularization factor lambda should be in the interval [0,1]")
	ws <- normalize_weights(ws)
	ws <- (1-lambda)*ws + lambda*mean(ws)
}

