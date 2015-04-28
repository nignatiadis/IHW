groups_by_filter <- function(filter_statistic, nbins){
	rfs <- rank(filter_statistic, ties.method="first")/length(filter_statistic)
	as.factor(ceiling( rfs* nbins))
}

# given list of adjusted p-values and original p-values
# calculate actual threshold! 
padj_to_threshold <- function(padj, pvals, alpha){
	filtered_pvals <- pvals[padj <= alpha]
	ifelse(length(filtered_pvals) == 0, 0, max(filtered_pvals))
}


mydiv <- function(x,y) ifelse(x == 0, 0, 
						 ifelse(y==0, 1, pmin(x/y,1)))

# returns a threshold!
filter_pvals_for_optim <- function(pvals, alpha, ntests=length(pvals)){
	nrjs <- sum(p.adjust(pvals, method="BH", n = ntests) < alpha)
	min(1, 10*nrjs/length(pvals))
}

get_bh_threshold <- function(pvals, alpha, mtests = length(pvals)){
  m <- length(pvals)
  pvals <- sort(pvals)
  prejected <- which(pvals <= (1:m)/mtests*alpha)
  ifelse(length(prejected)==0, 0, pvals[prejected[which.max(prejected)]])
}




get_bh_thresholds <- function(unadj_p, filterstat, nbins, alpha){
	t <- get_bh_threshold(unadj_p, alpha)
	grps <- groups_by_filter(filterstat, nbins)
	pv_list <- split(unadj_p, grps)
	sapply(pv_list, function(ps) max(0, ps[ps <=t]))
}

get_wbh_weights<- function(obj){
	weighted_pv <- pvalues(obj)/weights(obj, levels_only=F)
	t <- get_bh_threshold(weighted_pv, alpha(obj))
	m_groups <- table(groups_factor(obj))
	grps <- groups_factor(obj)
	pv_list <- split(pvalues(obj), grps)
	ts <- sapply(pv_list, function(ps) max(0, ps[ps <=t]))
	ts*sum(m_groups)/sum(ts*m_groups)
}



my_grenander <- function(unadj_p,t, distrib=F){

  	x0 <- fndr.cutoff(unadj_p, "pvalue")	
  	cf.out <- censored.fit(x=unadj_p, cutoff=x0, statistic="pvalue")

    scale.param = NULL

	eta0 = cf.out[1,3]

 	# determine cumulative empirical distribution function (pvalues)
 	ee <- ecdf.pval(unadj_p, eta0=eta0)

	g.pval <- grenander(ee)

	# mixture density and CDF  

 	F.pval = approxfun( g.pval$x.knots,  g.pval$F.knots, method="linear", 
           yleft=0, yright=g.pval$F.knots[length(g.pval$F.knots)])
	F0.pval = function(x) return( ifelse(x > 1, 1, ifelse(x < 0, 0, x )) )

  
  	Fdr.pval = function(p) pmin( eta0*p / F.pval(p), 1) # eta0*F0/ F

  if (distrib){
    return(list(distrib=F.pval(t), pi_zero=eta0, distrib_fun = F.pval))
  } else {
    return(  min(1,Fdr.pval(t), na.rm=T))
  }
}

ecdf.pval <- function (x, eta0=1) 
{
    # compute empirical CDF as usual
    x = sort(x)
    n = length(x)
    if (n < 1) 
        stop("'x' must have 1 or more non-missing values")
    vals = sort(unique(x))
    F.raw = cumsum(tabulate(match(x, vals)))/n
    
    # control upper bound of F:
    # make sure that the maximum slope of (Grenander) F is eta0
    F.raw = pmin(F.raw, 1-eta0*(1-vals) ) 

    # control lower bound of F: 
    # make sure that (Grenander F) >= eta0*vals
    F.raw = pmax(F.raw, eta0*vals) 

    # if necessary add an atom at 1 to make it a proper CDF
    if (vals[length(vals)] != 1)
    {
       F.raw = c(F.raw, 1)
       vals = c(vals, 1)
    }

    # if necessary also add an atom at 0 with weight zero to get support [0,1]
    if (vals[1] != 0)
    {
       F.raw = c(0, F.raw)
       vals = c(0, vals)
    }

    # finally, modify F such that the last slope of the Grenander F 
    # is *exactly* eta0
    i = length(vals)-1
    F.raw[i] = 1-eta0*(1-vals[i])
    
    rval <- approxfun(vals, F.raw, 
        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) = c("ecdf", "stepfun", class(rval))
    attr(rval, "call") <- sys.call()
    rval
}

  