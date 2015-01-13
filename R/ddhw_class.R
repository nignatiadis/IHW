setClass("ddhw",
         	slots = list( 
           		df = "data.frame",
           		weights = "numeric",
           		thresholds = "numeric",
           		alpha = "numeric",
           		solver_information= "list"))






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
rejections.ddhw <- function(object){
	groups <- groups_factor(object)

	# use two different methods to check for number of rejections
	# as internal test

	rjs1 <- sum(pvalues(object) <= thresholds(object, levels_only=F))
	rjs2 <- sum(adj_pvalues(object) <= alpha(object))

	if (rjs1 != rjs2) stop("Incosistent number of rejections")
	rjs1
}

setGeneric("rejections", function(object,value) standardGeneric("rejections"))


setMethod("rejections", signature(object="ddhw"),
          rejections.ddhw)


############# validity ##########################################################
setValidity( "ddhw", function( object ) {
	return(TRUE)
} )
