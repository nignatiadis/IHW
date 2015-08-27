#' An S4 class to represent the DDHW output.
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
#'    set.seed(1)
#'    X   <- runif(20000, min=0.5, max=4.5) #covariate
#'    H   <- rbinom(20000,1,0.1)            #hypothesis true or false
#'    Z   <- rnorm(20000, H*X)              #Z-score
#'    pvalue <- 1-pnorm(Z)                  #pvalue
#'    ddhw_res <- ddhw(pvalue, X, .1)
#'    rejections(ddhw_res)
#'    colnames(as.data.frame(ddhw_res))
#'
#' @seealso ddhw, plot_ddhw
#' @import methods

ddhwResult <- setClass("ddhwResult",
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
adj_pvalues.ddhwResult <- function(object){
  object@df$padj
}

#' @rdname ddhwResult-class
setGeneric("adj_pvalues", function(object) standardGeneric("adj_pvalues"))

#' @describeIn ddhwResult Extract adjusted pvalues
#' @export
setMethod("adj_pvalues", signature(object="ddhwResult"),
          adj_pvalues.ddhwResult)

#--------------------------- weights extraction --------------------------------------------------------------------#


weights.ddhwResult <-function(object, levels_only = FALSE){
  if (levels_only) {
    object@weights
  } else {
    object@df$weights #TODO: Storing redundant information right now
  }
}

#' @param object,x A ddhwResult object as returned by a call to ddhw(...)
#' @param levels_only Boolean, if FALSE, return a vector of weights (thresholds) with one weight
#'    (threshold) for each hypothesis, otherwise return a nfolds x nbins matrix of weights/thresholds
#' @param ... Parameters passed in to individual methods
#' @describeIn ddhwResult Extract weights
#' @export
setMethod("weights", signature(object="ddhwResult"),
          weights.ddhwResult)


#--------------------------- threshold extraction -------------------------------------------------------------------#


thresholds.ddhwResult <-function(object, levels_only = FALSE){
  t <- get_bh_threshold(na.exclude(weighted_pvalues(object)), alpha(object))
  t*weights(object, levels_only = levels_only)
}

#' @rdname ddhwResult-class
setGeneric("thresholds", function(object,...) standardGeneric("thresholds"))

#' @describeIn ddhwResult Calculate ddhw thresholds
#' @export
setMethod("thresholds", signature(object="ddhwResult"),
          thresholds.ddhwResult)



#--------------------------- p-value extraction ---------------------------------------------------------------------#

pvalues.ddhwResult <- function(object){
  object@df$pvalue
}

#' @rdname ddhwResult-class
setGeneric("pvalues", function(object) standardGeneric("pvalues"))

#' @describeIn ddhwResult Extract pvalues
#' @export
setMethod("pvalues", signature(object="ddhwResult"),
          pvalues.ddhwResult)


#--------------------------- weighted p-value extraction -------------------------------------------------------------#

weighted_pvalues.ddhwResult <- function(object){
  object@df$weighted_pvalue
}

#' @rdname ddhwResult-class
setGeneric("weighted_pvalues", function(object) standardGeneric("weighted_pvalues"))

#' @describeIn ddhwResult Extract weighted pvalues
#' @export
setMethod("weighted_pvalues", signature(object="ddhwResult"),
          weighted_pvalues.ddhwResult)


#---------------------------  filter statistic extraction ----------------------------------------------------------#

filter_statistics.ddhwResult <- function(object){
  object@df$filter_statistic
}

#' @rdname ddhwResult-class
setGeneric("filter_statistics", function(object) standardGeneric("filter_statistics"))

#' @describeIn ddhwResult Extract filter statistics
#' @export
setMethod("filter_statistics", signature(object="ddhwResult"),
          filter_statistics.ddhwResult)


#----------------- extract stratification variable----------------------------------------------------------------#
groups_factor.ddhwResult <- function(object){
	object@df$group
}

#' @rdname ddhwResult-class
setGeneric("groups_factor", function(object) standardGeneric("groups_factor"))


#' @describeIn ddhwResult Extract factor of stratification (grouping) variable
#' @export
setMethod("groups_factor", signature(object="ddhwResult"),
          groups_factor.ddhwResult)






#----------------- nominal alpha extraction ----------------------------------------------------------------------#
alpha.ddhwResult <-function(object) object@alpha

#' @rdname ddhwResult-class
setGeneric("alpha", function(object) standardGeneric("alpha"))

#' @describeIn ddhwResult Extract nominal significance (alpha) level
#' @export
setMethod("alpha", signature(object="ddhwResult"),
          alpha.ddhwResult)

#----------------- rejections ------------------------------------------------------------------------------------#


# TODO: Extend this to groupwise calculation (i.e. rejections per stratum)
rejections.ddhwResult <- function(object){
  sum(rejected_hypotheses(object), na.rm=TRUE)
}

#' @rdname ddhwResult-class
setGeneric("rejections", function(object,...) standardGeneric("rejections"))

#' @describeIn ddhwResult Total number of rejected hypotheses by DDHW procedure
#' @export
setMethod("rejections", signature(object="ddhwResult"),
          rejections.ddhwResult)



rejected_hypotheses.ddhwResult <- function(object){
  adj_pvalues(object) <= alpha(object)
}


#' Generic function for multiple testing procedures
#' Returns boolean vector of rejected hypotheses
#' @param  object Multiple testing procedure object (e.g. ddhwResult)
#' @param ... Arguments passed to individual methods
#' @return boolean vector of rejected hypotheses
#' @export
rejected_hypotheses <- function(object,...) UseMethod("rejected_hypotheses")

#' @describeIn ddhwResult Get a boolean vector of the rejected hypotheses
#' @export
setMethod("rejected_hypotheses", signature(object="ddhwResult"),
          rejected_hypotheses.ddhwResult)


#--------------- convenience methods ------------------------------------------------------------------------------#

#' @describeIn ddhwResult Convert ddhwResult object to data frame
#' @export
setMethod("as.data.frame", signature=c("ddhwResult"),
      function(x){
        x@df
      })

#' @describeIn ddhwResult Convenience method to show ddhwResult object
#' @importFrom methods show
#' @export
setMethod("show", signature(object="ddhwResult"), function(object) {
  cat("ddhwResult object with", nrow(ddhw_res@df),"hypothesis tests \n")
  cat("Nominal FDR control level:", alpha(ddhw_res),"\n")
  cat("Split into", ddhw_res@nbins,"bins, based on an", ddhw_res@filter_statistic_type, "covariate\n")
})


#------------------ not exported stuff ----------------------------------------------------------------------------#

##### FDR estimate #############################################################
plugin_fdr.ddhwResult <- function(object) {
  ts <- thresholds(object)
  m_groups <- table(groups_factor(object))
  sum(ts*m_groups)/rejections(object, method="thresholds")
}

setGeneric("plugin_fdr", function(object,...) standardGeneric("plugin_fdr"))


setMethod("plugin_fdr", signature(object="ddhwResult"),
          plugin_fdr.ddhwResult)


##### #############################################################
stratification_breaks.ddhwResult <- function(object) {
  ts <- thresholds(object)
  groups <- groups_factor(object)
  filterstat_list <- split(filter_statistics(object), groups)
  filterstats <- sapply(filterstat_list, max)
  filterstats
}

setGeneric("stratification_breaks", function(object,...) standardGeneric("stratification_breaks"))


setMethod("stratification_breaks", signature(object="ddhwResult"),
          stratification_breaks.ddhwResult)


######## temporary: number of pvals in each stratum #############################
stratum_sizes <- function(object) table(groups_factor(object))

############# validity ##########################################################
setValidity( "ddhwResult", function( object ) {
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


