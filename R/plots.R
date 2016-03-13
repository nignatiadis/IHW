# plot weights VS group (categorical)


#' Convenience function to plot weights learned by ihw method
#'
#' @param x ihwResult Object
#' @param x_axis Character, currently only "group" supported.
#' @param scale  Character, "ordinal" or "nominal"
#'
#' @return A ggplot2 object.
#'
#' @examples
#'
#'    save.seed <- .Random.seed; set.seed(1)
#'    X   <- runif(20000, min=0.5, max=4.5) #covariate
#'    H   <- rbinom(20000,1,0.1)            #hypothesis true or false
#'    Z   <- rnorm(20000, H*X)              #Z-score
#'    .Random.seed <- save.seed
#'    pvalue <- 1-pnorm(Z)                  #pvalue
#'    ihw_res <- ihw(pvalue, X, .1)
#'    plot(ihw_res)
#'
#' @export
setMethod("plot", signature="ihwResult",
	function(x, x_axis="group", scale=covariate_type(x)){
		if (requireNamespace("ggplot2", quietly=TRUE) & (requireNamespace("scales", quietly=TRUE))){
			ws <- weights(x, levels_only=TRUE)
			group_levels <- levels(groups_factor(x))
			folds  <- factor(1:nfolds(x))
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
								ggplot2::geom_line() +
								ggplot2::scale_x_continuous(breaks= scales::pretty_breaks())
				} else if (x_axis == "ranks"){

				}
			}
			return(plt)
		} else {
			stop("ihwResult plotting method requires ggplot2 and scales packages.")
		}
	})
# plot weights vs Rank (ordinal)
# plot weights vs filter statistic value
# requireNamespace("x", quietly = TRUE)