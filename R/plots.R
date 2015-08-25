# plot weights VS group (categorical)


#' Convenience function to plot weights learned by DDHW method
#'
#' @param x ddhwResult Object
#' @param ... additional arguments passed on to individual methods
setGeneric("plot_ddhw", function(x, ...) standardGeneric("plot_ddhw"))

#' @describeIn plot_ddhw
#' @param x_axis Character, currently only "group" supported.
#' @param scale  Character, "ordinal" or "nominal"
#'
#' @return A ggplot2 object.
#'
#' @examples
#'
#'    set.seed(1)
#'    X   <- runif(20000, min=0.5, max=4.5) #covariate
#'    H   <- rbinom(20000,1,0.1)            #hypothesis true or false
#'    Z   <- rnorm(20000, H*X)              #Z-score
#'    pvalue <- 1-pnorm(Z)                  #pvalue
#'    ddhw_res <- ddhw(pvalue, X, .1)
#'	  plot_ddhw(ddhw_res)
#'
#' @importFrom graphics plot
#' @export
setMethod("plot_ddhw", signature="ddhwResult", 
	function(x, x_axis="group", scale=x@filter_statistic_type){
		if (requireNamespace("ggplot2", quietly=TRUE)){
			ws <- weights(x, levels_only=TRUE)
			group_levels <- levels(groups_factor(x))
			folds  <- factor(1:x@nfolds)
			df <- expand.grid(stratum=1:nlevels(groups_factor(x)),
							  fold=folds)
			df$group <- group_levels[df$stratum]
			df$weight <- mapply(function(x,y) ws[x,y], df$stratum, df$fold)


			if (scale == "nominal"){
				if (x_axis != "group"){
					warning("x_axis argument ignored since scale is nominal, i.e. we cannot rank")
				}
				plt <- ggplot2::ggplot(df, ggplot2::aes_string(x='fold',y='weight',fill='fold')) +
							ggplot2::geom_bar(stat = "identity", position="dodge") +
							ggplot2::facet_wrap(~group)
			} else if (scale=="ordinal"){
				if (x_axis == "group"){
					plt <-	ggplot2::ggplot(df, ggplot2::aes_string(x='stratum', y='weight', col='fold'))+
								ggplot2::geom_point() +
								ggplot2::geom_line()
				} else if (x_axis == "ranks"){

				}
			}
			return(plt)
		} else {
			stop("ddhwResult plotting method requires ggplot2.")
		}
	})
# plot weights vs Rank (ordinal)
# plot weights vs filter statistic value
# requireNamespace("x", quietly = TRUE)