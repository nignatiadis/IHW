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

