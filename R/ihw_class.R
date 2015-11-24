#' An S4 class to represent the ihw output.
#'
#' @slot df  A data frame which collects the input data such as the vector of p values and the filter statistics, the group assignment, as well as outputs (weighted p-values, adjusted p-values)
#' @slot weights  A (nbins X nfolds) matrix of the weight assigned to each stratum
#' @slot alpha   Numeric, the nominal significance level at which the FDR is to be controlled
#' @slot nbins  Integer, number of distinct levels into which the hypotheses were stratified
#' @slot nfolds  Integer, number of folds for pre-validation procedure
#' @slot regularization_term Numeric vector, the final value of the regularization parameter within each fold
#' @slot penalty  Character, "uniform deviation" or "total variation"
#' @slot filter_statistic_type  Character, "ordinal" or "nominal"
#' @slot reg_path_information  A data frame, information about the whole regularization path. (Currently empty data frame)
#' @slot solver_information  A list, solver specific output, e.g. were all subproblems solved to optimality? (Currently empty list)
#'
#' @examples
#'
#' set.seed(1)
#' X   <- runif(20000, min=0.5, max=4.5) #covariate
#' H   <- rbinom(20000,1,0.1)            #hypothesis true or false
#' Z   <- rnorm(20000, H*X)              #Z-score
#' pvalue <- 1-pnorm(Z)                  #pvalue
#' ihw_res <- ihw(pvalue, X, .1)
#' rejections(ihw_res)
#' colnames(as.data.frame(ihw_res))
#'
#' @seealso ihw, plot_ihw
#' @import methods

ihwResult <- setClass("ihwResult",
         	      slots = list(
           		     df = "data.frame",
           		     weights = "matrix",
           		     alpha = "numeric",
                   nbins = "integer",
                   nfolds = "integer",
                   regularization_term = "numeric",
                   penalty = "character",
                   filter_statistic_type = "character",
                   reg_path_information = "data.frame",
           		     solver_information= "list"))

#-----------------------------adjusted p-values extraction---------------------------------------------------------#
adj_pvalues.ihwResult <- function(object){
  object@df$adj_pvalue
}

#' @rdname ihwResult-class
setGeneric("adj_pvalues", function(object) standardGeneric("adj_pvalues"))

#' @describeIn ihwResult Extract adjusted pvalues
#' @export
setMethod("adj_pvalues", signature(object="ihwResult"),
          adj_pvalues.ihwResult)

#--------------------------- weights extraction --------------------------------------------------------------------#


weights.ihwResult <-function(object, levels_only = FALSE){
  if (levels_only) {
    object@weights
  } else {
    object@df$weight #TODO: Storing redundant information right now
  }
}

#' @param object,x A ihwResult object as returned by a call to ihw(...)
#' @param levels_only Boolean, if FALSE, return a vector of weights (thresholds) with one weight
#'    (threshold) for each hypothesis, otherwise return a nfolds x nbins matrix of weights/thresholds
#' @param ... Parameters passed in to individual methods
#' @describeIn ihwResult Extract weights
#' @export
setMethod("weights", signature(object="ihwResult"),
          weights.ihwResult)


#--------------------------- threshold extraction -------------------------------------------------------------------#


thresholds.ihwResult <-function(object, levels_only = FALSE){
  t <- get_bh_threshold(na.exclude(weighted_pvalues(object)), alpha(object))
  t*weights(object, levels_only = levels_only)
}

#' @rdname ihwResult-class
setGeneric("thresholds", function(object,...) standardGeneric("thresholds"))

#' @describeIn ihwResult Calculate ihw thresholds
#' @export
setMethod("thresholds", signature(object="ihwResult"),
          thresholds.ihwResult)



#--------------------------- p-value extraction ---------------------------------------------------------------------#

pvalues.ihwResult <- function(object){
  object@df$pvalue
}

#' @rdname ihwResult-class
setGeneric("pvalues", function(object) standardGeneric("pvalues"))

#' @describeIn ihwResult Extract pvalues
#' @export
setMethod("pvalues", signature(object="ihwResult"),
          pvalues.ihwResult)


#--------------------------- weighted p-value extraction -------------------------------------------------------------#

weighted_pvalues.ihwResult <- function(object){
  object@df$weighted_pvalue
}

#' @rdname ihwResult-class
setGeneric("weighted_pvalues", function(object) standardGeneric("weighted_pvalues"))

#' @describeIn ihwResult Extract weighted pvalues
#' @export
setMethod("weighted_pvalues", signature(object="ihwResult"),
          weighted_pvalues.ihwResult)


#---------------------------  filter statistic extraction ----------------------------------------------------------#

filter_statistics.ihwResult <- function(object){
  object@df$filter_statistic
}

#' @rdname ihwResult-class
setGeneric("filter_statistics", function(object) standardGeneric("filter_statistics"))

#' @describeIn ihwResult Extract filter statistics
#' @export
setMethod("filter_statistics", signature(object="ihwResult"),
          filter_statistics.ihwResult)


#----------------- extract stratification variable----------------------------------------------------------------#
groups_factor.ihwResult <- function(object){
	object@df$group
}

#' @rdname ihwResult-class
setGeneric("groups_factor", function(object) standardGeneric("groups_factor"))


#' @describeIn ihwResult Extract factor of stratification (grouping) variable
#' @export
setMethod("groups_factor", signature(object="ihwResult"),
          groups_factor.ihwResult)






#----------------- nominal alpha extraction ----------------------------------------------------------------------#
alpha.ihwResult <-function(object) object@alpha

#' @rdname ihwResult-class
setGeneric("alpha", function(object) standardGeneric("alpha"))

#' @describeIn ihwResult Extract nominal significance (alpha) level
#' @export
setMethod("alpha", signature(object="ihwResult"),
          alpha.ihwResult)

#----------------- rejections ------------------------------------------------------------------------------------#


# TODO: Extend this to groupwise calculation (i.e. rejections per stratum)
rejections.ihwResult <- function(object){
  sum(rejected_hypotheses(object), na.rm=TRUE)
}

#' @rdname ihwResult-class
setGeneric("rejections", function(object,...) standardGeneric("rejections"))

#' @describeIn ihwResult Total number of rejected hypotheses by ihw procedure
#' @export
setMethod("rejections", signature(object="ihwResult"),
          rejections.ihwResult)



rejected_hypotheses.ihwResult <- function(object){
  adj_pvalues(object) <= alpha(object)
}

#' @rdname ihwResult-class
setGeneric("rejected_hypotheses", function(object,...) standardGeneric("rejected_hypotheses"))

#' @describeIn ihwResult Get a boolean vector of the rejected hypotheses
#' @export
setMethod("rejected_hypotheses", signature(object="ihwResult"),
          rejected_hypotheses.ihwResult)


#--------------- convenience methods ------------------------------------------------------------------------------#

#' @rdname ihwResult-class
#' @export
as.data.frame.ihwResult <-function(x,row.names=NULL, optional=FALSE, ...){
        x@df
      }

#' @describeIn ihwResult Coerce ihwResult to data frame
#' @export
setMethod("as.data.frame", "ihwResult",as.data.frame.ihwResult)

#' @describeIn ihwResult Convenience method to show ihwResult object
#' @importFrom methods show
#' @export
setMethod("show", signature(object="ihwResult"), function(object) {
  cat("ihwResult object with", nrow(ihw_res@df),"hypothesis tests \n")
  cat("Nominal FDR control level:", alpha(ihw_res),"\n")
  cat("Split into", ihw_res@nbins,"bins, based on an", ihw_res@filter_statistic_type, "covariate\n")
})


#------------------ not exported stuff ----------------------------------------------------------------------------#

##### FDR estimate #############################################################
plugin_fdr.ihwResult <- function(object) {
  ts <- thresholds(object)
  m_groups <- table(groups_factor(object))
  sum(ts*m_groups)/rejections(object, method="thresholds")
}

setGeneric("plugin_fdr", function(object,...) standardGeneric("plugin_fdr"))


setMethod("plugin_fdr", signature(object="ihwResult"),
          plugin_fdr.ihwResult)


##### #############################################################
stratification_breaks.ihwResult <- function(object) {
  ts <- thresholds(object)
  groups <- groups_factor(object)
  filterstat_list <- split(filter_statistics(object), groups)
  filterstats <- sapply(filterstat_list, max)
  filterstats
}

setGeneric("stratification_breaks", function(object,...) standardGeneric("stratification_breaks"))


setMethod("stratification_breaks", signature(object="ihwResult"),
          stratification_breaks.ihwResult)


######## temporary: number of pvals in each stratum #############################
stratum_sizes <- function(object) table(groups_factor(object))

############# validity ##########################################################
setValidity( "ihwResult", function( object ) {
	return(TRUE)
} )


per_bin_fdrs <- function(obj) {
  ts <- thresholds(obj)
  groups <- groups_factor(obj)
  pvals <- pvalues(obj)
  pv_list <- split(pvals, groups)
  print(sapply(pv_list,length))
  fdrs <- mapply(function(t,pvec) length(pvec)*t/max(1, sum(pvec <= t)), ts, pv_list)
  return(fdrs)
}


