setClass("ddhw",
         	slots = list(
           		df = "data.frame",
           		weights = "numeric",
           		thresholds = "numeric",
           		alpha = "numeric",
              nbins = "integer",
              regularization_term = "numeric",
              penalty = "character",
              reg_path_information = "data.frame",
           		solver_information= "list"))

#----------------- reg_path_information--------------------------------------#

setGeneric("reg_path_information<-", function(object,value) standardGeneric("reg_path_information<-"))

setReplaceMethod("reg_path_information", signature(object="ddhw", value="data.frame"),
                 function( object, value ) {
                    object@reg_path_information <- value
                    validObject(object)
                    object
                 })




########## filter_statistic ###################################################
filter_statistics.ddhw <- function(object){
  object@df$filter_statistic
}

setGeneric("filter_statistics", function(object,value) standardGeneric("filter_statistics"))


setMethod("filter_statistics", signature(object="ddhw"),
          filter_statistics.ddhw)


setGeneric("filter_statistics<-", function(object,value) standardGeneric("filter_statistics<-"))

setReplaceMethod("filter_statistics", signature(object="ddhw", value="numeric"),
                 function( object, value ) {
                  	object@df$filter_statistic <- value
                  	validObject(object)
                  	object
                 })


##########    groups   ###################################################
groups_factor.ddhw <- function(object){
	object@df$group
}

setGeneric("groups_factor", function(object,value) standardGeneric("groups_factor"))


setMethod("groups_factor", signature(object="ddhw"),
          groups_factor.ddhw)


################### weights #################################################
weights.ddhw <-function(object, levels_only = T){
	if (levels_only) {
		object@weights
	} else {
		object@weights[groups_factor(object)]
	}
}

setMethod("weights", signature(object="ddhw"),
          weights.ddhw)

setGeneric("weights<-", function(object,value) standardGeneric("weights<-"))

setReplaceMethod("weights", signature(object="ddhw", value="numeric"),
                 function( object, value ) {
                  	object@weights <- value
                  	validObject(object)
                  	object
                 })

################### alpha #################################################
alpha.ddhw <-function(object) object@alpha

setGeneric("alpha", function(object,value) standardGeneric("alpha"))

setMethod("alpha", signature(object="ddhw"),
          alpha.ddhw)


################### thresholds #################################################
thresholds.ddhw <-function(object, levels_only = T){
	if (levels_only) {
		object@thresholds
	} else {
		object@thresholds[groups_factor(object)]
	}
}

setGeneric("thresholds", function(object,value,...) standardGeneric("thresholds"))

setMethod("thresholds", signature(object="ddhw"),
          thresholds.ddhw)


##########    pvalues   ###################################################
pvalues.ddhw <- function(object){
	object@df$pvalue
}

setGeneric("pvalues", function(object,value) standardGeneric("pvalues"))


setMethod("pvalues", signature(object="ddhw"),
          pvalues.ddhw)

##########     adjusted pvalues   ###############################################
adj_pvalues.ddhw <- function(object){
	object@df$adj_p
}

setGeneric("adj_pvalues", function(object,value) standardGeneric("adj_pvalues"))


setMethod("adj_pvalues", signature(object="ddhw"),
          adj_pvalues.ddhw)

########## rejections ###################################################
rejections.ddhw <- function(object, method="adjpvals", groupwise=F){
	groups <- groups_factor(object)

	# use two different methods to check for number of rejections
	# as internal test

  if (groupwise){
    rjs1 <- sapply(split(adj_pvalues(object), groups_factor(object)),
           function(adj_p) sum(adj_p <= alpha(object)))
                          
    rjs2 <- mapply(function(pv, t) sum(pv <=t), split(pvalues(object), groups_factor(object))
                          ,thresholds(object))

    } else {

	   rjs1 <- sum(pvalues(object) <= thresholds(object, levels_only=F))
	   rjs2 <- sum(adj_pvalues(object) <= alpha(object))
  }

	
  #max(rjs1,rjs2)
  if (method=="adjpvals"){
    return(rjs2)
  } else if (method=="thresholds"){
    return(rjs1)
  }
}

setGeneric("rejections", function(object,value,...) standardGeneric("rejections"))


setMethod("rejections", signature(object="ddhw"),
          rejections.ddhw)


##### FDR estimate #############################################################
plugin_fdr.ddhw <- function(object) {
  ts <- thresholds(object)
  m_groups <- table(groups_factor(object))
  sum(ts*m_groups)/rejections(object, method="thresholds")
}

setGeneric("plugin_fdr", function(object, value,...) standardGeneric("plugin_fdr"))


setMethod("plugin_fdr", signature(object="ddhw"),
          plugin_fdr.ddhw)


##### #############################################################
stratification_breaks.ddhw <- function(object) {
  ts <- thresholds(object)
  groups <- groups_factor(object)
  filterstat_list <- split(filter_statistics(object), groups)
  filterstats <- sapply(filterstat_list, max)
  filterstats
}

setGeneric("stratification_breaks", function(object, value,...) standardGeneric("stratification_breaks"))


setMethod("stratification_breaks", signature(object="ddhw"),
          stratification_breaks.ddhw)


######## temporary: number of pvals in each stratum #############################
stratum_sizes <- function(object) table(groups_factor(object))

############# validity ##########################################################
setValidity( "ddhw", function( object ) {
	return(TRUE)
} )
