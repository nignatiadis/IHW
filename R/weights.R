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

