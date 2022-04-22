#' An S4 class to represent the ihw output.
#'
#' @slot df  A data.frame that collects the input data, including the vector of p values and the covariate, the group assignment, as well as outputs (weighted p-values, adjusted p-values)
#' @slot weights  A (nbins X nfolds) matrix of the weight assigned to each stratum
#' @slot alpha   Numeric, the nominal significance level at which the FDR is to be controlled
#' @slot nbins  Integer, number of distinct levels into which the hypotheses were stratified
#' @slot nfolds  Integer, number of folds for pre-validation procedure
#' @slot regularization_term Numeric vector, the final value of the regularization parameter within each fold
#' @slot m_groups Integer vector, number of hypotheses tested in each stratum
#' @slot penalty  Character, "uniform deviation" or "total variation"
#' @slot covariate_type  Character, "ordinal" or "nominal"
#' @slot adjustment_type  Character, "BH" or "bonferroni"
#' @slot reg_path_information  A data.frame, information about the whole regularization path. (Currently not used, thus empty)
#' @slot solver_information  A list, solver specific output, e.g. were all subproblems solved to optimality? (Currently empty list)
#' @slot stratification_method Character, "quantiles" or "forest" or "none"
#' @slot weight_matrices_forest list of list of vectors, equivalent of weights, but for stratification method forest
#' @slot n_censor_thres  Integer, number of censoring thresholds tau to be considered for stratification method "forest"
#' @slot ntrees  Integer, see same parameter in \code{\link[randomForestSRC]{rfsrc}} for stratification method "forest"
#' 
#' @return The different methods applied to an ihwResult object can return the following:
#' @return 1) A vector
#'     of length equal to the number of hypotheses tested (e.g. the adjusted p-value or the weight
#'     of each hypothesis).
#' @return 2) A matrix of dimension equal to nfolds x nbins (e.g. the weight of
#'     each stratum, fold combination, set by specifying levels_only=TRUE).
#' @return 3) A vector of length 1 (usually a parameter of the ihwResult object such as nfolds or
#'     the total number of rejections).
#' @return 4) A data.frame (as.data.frame) or just console output (show) for the extended Base generics.
#' @return See section below for the individual methods.
#' @examples
#'
#' save.seed <- .Random.seed; set.seed(1)
#' X   <- runif(n = 20000, min = 0.5, max = 4.5)       # Covariate
#' # Is the null hypothesis (mean=0) true or false ?
#' H   <- rbinom(n = length(X), size = 1, prob = 0.1)  
#' Z   <- rnorm(n = length(X), mean = H * X)           # Z-score
#' .Random.seed <- save.seed
#'
#' pvalue <- 1 - pnorm(Z)                              # pvalue
#' ihw_res <- ihw(pvalue, covariates = X, alpha = 0.1)
#' rejections(ihw_res)
#' colnames(as.data.frame(ihw_res))
#'
#' @seealso ihw, plot,ihwResult-method
#' @import methods

ihwResult <- setClass("ihwResult",
         	      slots = list(
           		     df = "data.frame",
           		     weights = "matrix",
           		     alpha = "numeric",
                   nbins = "integer",
                   nfolds = "integer",
                   regularization_term = "matrix",
                   m_groups = "integer",
                   penalty = "character",
                   covariate_type = "character",
                   adjustment_type = "character",
                   reg_path_information = "data.frame",
           		     solver_information= "list",
           		     stratification_method = "character",
           		     weight_matrices_forest = "list",
           		     n_censor_thres = "integer",
           		     ntrees = "integer"))

#-----------------------------adjusted p-values extraction---------------------------------------------------------#
adj_pvalues_ihwResult <- function(object){
  object@df$adj_pvalue
}

#' @rdname ihwResult-class
setGeneric("adj_pvalues", function(object) standardGeneric("adj_pvalues"))

#' @describeIn ihwResult Extract adjusted pvalues
#' @export
setMethod("adj_pvalues", signature(object="ihwResult"),
          adj_pvalues_ihwResult)

#--------------------------- weights extraction --------------------------------------------------------------------#


weights_ihwResult <-function(object, levels_only = FALSE){
  if (levels_only) {
    if(object@stratification_method == "forest"){
      object@weight_matrices_forest
    }else{
      object@weights
    }
  } else {
    object@df$weight #TODO: Storing redundant information right now
  }
}

#' @param object,x A ihwResult object as returned by a call to ihw(...)
#' @param levels_only Logical, if FALSE, return a vector of weights (thresholds) with one weight
#'    (threshold) for each hypothesis, otherwise return a nfolds x nbins matrix of weights (thresholds)
#' @param row.names,optional See ?base::as.data.frame for a description of these arguments.
#' @param ... Parameters passed in to individual methods
#' @describeIn ihwResult Extract weights
#' @export
setMethod("weights", signature(object="ihwResult"),
          weights_ihwResult)


#--------------------------- threshold extraction -------------------------------------------------------------------#


thresholds_ihwResult <-function(object, levels_only = FALSE){
  if(object@stratification_method == "forest" & levels_only == TRUE){
    stop("threshold() extraction currently not available for stratification_method = forest and levels_only == TRUE")
  }else if(object@stratification_method == "forest" & levels_only == FALSE){
    mtests <- length(na.exclude(weighted_pvalues(object)))
  }else{
    mtests <- sum(as.numeric(m_groups(object)))
  }
  
  t <- get_bh_threshold(na.exclude(weighted_pvalues(object)), alpha(object),
                        mtests = mtests)
  t*weights(object, levels_only = levels_only)
}

#' @rdname ihwResult-class
setGeneric("thresholds", function(object,...) standardGeneric("thresholds"))

#' @describeIn ihwResult Calculate ihw thresholds
#' @export
setMethod("thresholds", signature(object="ihwResult"),
          thresholds_ihwResult)



#--------------------------- p-value extraction ---------------------------------------------------------------------#

pvalues_ihwResult <- function(object){
  object@df$pvalue
}

#' @rdname ihwResult-class
setGeneric("pvalues", function(object) standardGeneric("pvalues"))

#' @describeIn ihwResult Extract pvalues
#' @export
setMethod("pvalues", signature(object="ihwResult"),
          pvalues_ihwResult)


#--------------------------- weighted p-value extraction -------------------------------------------------------------#

weighted_pvalues_ihwResult <- function(object){
  object@df$weighted_pvalue
}

#' @rdname ihwResult-class
setGeneric("weighted_pvalues", function(object) standardGeneric("weighted_pvalues"))

#' @describeIn ihwResult Extract weighted pvalues
#' @export
setMethod("weighted_pvalues", signature(object="ihwResult"),
          weighted_pvalues_ihwResult)


#---------------------------  covariate extraction ----------------------------------------------------------#

covariates_ihwResult <- function(object){ 
  covariate_names <- grepl("^covariate", colnames(object@df)) #TODO test this
  covariates <- object@df[ ,covariate_names, drop = TRUE]
  covariates
}

#' @rdname ihwResult-class
setGeneric("covariates", function(object) standardGeneric("covariates"))

#' @describeIn ihwResult Extract covariates
#' @export
setMethod("covariates", signature(object="ihwResult"),
          covariates_ihwResult)

#---------------------------- covariate type -----------------------------------------------------------------#
covariate_type_ihwResult <- function(object){
  object@covariate_type
}

#' @rdname ihwResult-class
setGeneric("covariate_type", function(object) standardGeneric("covariate_type"))

#' @describeIn ihwResult Extract type of covariate ("ordinal" or "nominal")
#' @export
setMethod("covariate_type", signature(object="ihwResult"),
          covariate_type_ihwResult)

#----------------- extract stratification variable------------------------------------------------------------#
groups_factor_ihwResult <- function(object){#
  if(object@stratification_method == "forest"){
    groups <- dplyr::select(object@df, dplyr::contains("tree")) 
    
    group_names <- grepl("tree", colnames(object@df)) #TODO test this
    groups <- object@df[ ,group_names, drop = TRUE]
    groups
    
  }else{
    object@df$group
  }
}

#' @rdname ihwResult-class
setGeneric("groups_factor", function(object) standardGeneric("groups_factor"))


#' @describeIn ihwResult Extract factor of stratification (grouping) variable
#' @export
setMethod("groups_factor", signature(object="ihwResult"),
          groups_factor_ihwResult)

#----------------- nfolds extraction ----------------------------------------------------------------------#
nfolds_ihwResult <- function(object) object@nfolds

#' @rdname ihwResult-class
setGeneric("nfolds", function(object) standardGeneric("nfolds"))

#' @describeIn ihwResult Extract number of folds
#' @export
setMethod("nfolds", signature(object="ihwResult"),
          nfolds_ihwResult)


#----------------- nbins extraction ----------------------------------------------------------------------#
nbins_ihwResult <- function(object){
  if(object@stratification_method == "forest"){
    stop("nbins() currently not available for stratification_method = forest")
  }
  object@nbins
} 

#' @rdname ihwResult-class
setGeneric("nbins", function(object) standardGeneric("nbins"))

#' @describeIn ihwResult Extract number of bins
#' @export
setMethod("nbins", signature(object="ihwResult"),
          nbins_ihwResult)


#----------------- nominal alpha extraction ----------------------------------------------------------------------#
alpha_ihwResult <-function(object) object@alpha

#' @rdname ihwResult-class
setGeneric("alpha", function(object) standardGeneric("alpha"))

#' @describeIn ihwResult Extract nominal significance (alpha) level
#' @export
setMethod("alpha", signature(object="ihwResult"),
          alpha_ihwResult)

#----------------- adjustment type extraction ----------------------------------------------------------------------#
adjustment_type_ihwResult <-function(object) object@adjustment_type

setGeneric("adjustment_type", function(object) standardGeneric("adjustment_type"))

setMethod("adjustment_type", signature(object="ihwResult"),
          adjustment_type_ihwResult)


#----------------- rejections ------------------------------------------------------------------------------------#


# TODO: Extend this to groupwise calculation (i.e. rejections per stratum)
rejections_ihwResult <- function(object){
  sum(rejected_hypotheses(object), na.rm=TRUE)
}

#' @rdname ihwResult-class
setGeneric("rejections", function(object,...) standardGeneric("rejections"))

#' @describeIn ihwResult Total number of rejected hypotheses by ihw procedure
#' @export
setMethod("rejections", signature(object="ihwResult"),
          rejections_ihwResult)



rejected_hypotheses_ihwResult <- function(object){
  adj_pvalues(object) <= alpha(object)
}

#' @rdname ihwResult-class
setGeneric("rejected_hypotheses", function(object,...) standardGeneric("rejected_hypotheses"))

#' @describeIn ihwResult Get a boolean vector of the rejected hypotheses
#' @export
setMethod("rejected_hypotheses", signature(object="ihwResult"),
          rejected_hypotheses_ihwResult)


#---------------- regularization terms ---------------------------
regularization_term_ihwResult <-function(object) object@regularization_term

#' @rdname ihwResult-class
setGeneric("regularization_term", function(object) standardGeneric("regularization_term"))

#' @describeIn ihwResult Extract vector of regularization parameters used for each stratum
#' @export
setMethod("regularization_term", signature(object="ihwResult"),
          regularization_term_ihwResult)

#---------------- m_groups --------------------------------------------------
m_groups_ihwResult <-function(object){
  if(object@stratification_method == "forest"){
    warning("For stratification_method = forest, m_groups has the structure of a nested list (fold-tree-group)") 
    groups <- groups_factor(object)
    fold <- object@df$fold
    m_groups <- lapply(groups, function(group) table(group, fold))
    m_groups
  }else{
    object@m_groups 
  }
} 

#' @rdname ihwResult-class
setGeneric("m_groups", function(object) standardGeneric("m_groups"))

#' @describeIn ihwResult Extract total number of hypotheses within each stratum
#' @export
setMethod("m_groups", signature(object="ihwResult"),
          m_groups_ihwResult)

#--------------- convenience methods ------------------------------------------------------------------------------#

#' @rdname ihwResult-class
as.data.frame_ihwResult <-function(x,row.names=NULL, optional=FALSE, ...){
        as.data.frame(x@df, row.names=row.names, optional=optional)
      }

#' @describeIn ihwResult Coerce ihwResult to data frame
#' @importFrom BiocGenerics as.data.frame
#' @export
setMethod("as.data.frame", "ihwResult",as.data.frame_ihwResult)

#' @describeIn ihwResult Return number of p-values
#' @export
setMethod("nrow", "ihwResult", function(x) nrow(x@df))

#' @describeIn ihwResult Convenience method to show ihwResult object
#' @importFrom methods show
#' @export
setMethod("show", signature(object="ihwResult"), function(object) {
  adj_type <- adjustment_type(object)
  if (adj_type == "BH"){
    typeI_error <- "FDR"
  } else if (adj_type == "bonferroni"){
    typeI_error <- "FWER"
  } else {
    stop("Adjustment method appears to be invalid, corrupted IHW object.")
  }
  cat("ihwResult object with", nrow(object),"hypothesis tests \n")
  cat("Nominal",  typeI_error, "control level:", alpha(object),"\n")
  cat("Split into", nbins(object),"bins, based on an", covariate_type(object), "covariate\n")
})



#------------------ not exported stuff ----------------------------------------------------------------------------#

##### FDR estimate #############################################################
plugin_fdr.ihwResult <- function(object) {
  if(object@stratification_method == "forest"){
    stop("plugin_fdr() currently not available for stratification_method = forest")
  }
  
  ts <- thresholds(object)
  m_groups <- table(groups_factor(object))
  sum(ts*m_groups)/rejections(object, method="thresholds")
}

setGeneric("plugin_fdr", function(object,...) standardGeneric("plugin_fdr"))


setMethod("plugin_fdr", signature(object="ihwResult"),
          plugin_fdr.ihwResult)


##### #############################################################
stratification_breaks.ihwResult <- function(object, fold = 1, tree = 1) {
  covariates <- covariates(object)
  if(is.matrix(covariates)){
    if(ncol(covariates) > 1) stop("stratification_breaks() currently not available for multi-dimensional covariates")
  }

  ts <- thresholds(object)
  groups <- groups_factor(object)
  if(object@stratification_method == "forest"){
    groups <- groups[ ,(fold - 1) * nfolds(object) + tree, drop = TRUE]
  } 
  filterstat_list <- split(covariates(object), groups)
  filterstats <- sapply(filterstat_list, max)
  sort(filterstats)
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
  if(obj@stratification_method == "forest"){
    stop("per_bin_fdrs() currently not available for stratification_method = forest")
  }
  
  ts <- thresholds(obj)
  groups <- groups_factor(obj)
  pvals <- pvalues(obj)
  pv_list <- split(pvals, groups)
  print(sapply(pv_list,length))
  fdrs <- mapply(function(t,pvec) length(pvec)*t/max(1, sum(pvec <= t)), ts, pv_list)
  return(fdrs)
}


