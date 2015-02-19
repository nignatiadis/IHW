# functions to assist in working with weights


total_variation <- function(ws){
    sum(abs(diff(ws)))
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

fix_thresholds <- function(ts, m_groups, rjs, alpha, lambda){
    # solve the LP
    # ts_reg >= ts
    # s.t. sum(ts_reg*m_groups) <= rjs*alpha
    # s.t. sum|m*ts_reg_{i+1} - m*ts_reg_{i}| <= lambda * sum(m_g * ts_reg_{g})

        

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


# mapping weights from n-simplex to n-sphere with spherical coordinates and vice versa
phi2w <- function(p){
	stopifnot(length(p)>=2)
    w <- numeric(length(p)+1)
    w[1] <- cos(p[1])
    w[-1] <- cumprod(sin(p))
    w[2:(length(w)-1)] <- w[2:(length(w)-1)]*cos(p[2:(length(w)-1)])
    w <- length(w) * w^2
    #check(w)
    w
}

  
  ## inverse of phi2w
w2phi <- function(w) {
    #check(w)
	w = sqrt(w/length(w))
    p = numeric(length(w)-1)
    s = 1
    for(i in seq(1, length(p), by=1)) {
      p[i] = acos(w[i]/s)
      s = s * sin(p[i])
    }
    p   
}


#' Generate random weights.
#'
#' this is equivalent to sampling uniformly from the n-simplex
#' see: http://geomblog.blogspot.de/2005/10/sampling-from-simplex.html
#'
#' @param npars Sampling will be from (npars-1)-simplex (embedded in R^npars)
#' @param N Number of random weight combinations to be generated

random_weights <- function(npars, N=1){
    mat <- matrix(ncol = npars, rexp(N*npars))
    mat <- t(apply(mat, 1, function(x) npars*x/sum(x)))
    if (N==1) as.numeric(mat) else mat
}

