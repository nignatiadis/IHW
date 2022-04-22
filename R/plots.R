# plot weights VS group (categorical)


#' Plot functions for IHW 
#'
#' See the vignette for usage examples. 
#'
#' @param x Object of class \code{ihwResult} 
#' @param what \code{Character}: "weights" or "decisionboundary"
#' @param x_axis \code{Character}: "group" or "covariate". Default is "group" if "what" is "weights", and "covariate" if "what" is "decisionboundary". 
#' @param scale  \code{Character}: "ordinal" or "nominal"
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#'
#'    save.seed <- .Random.seed; set.seed(1)
#'    X   <- runif(20000, min = 0.5, max = 4.5) # covariate
#'    H   <- rbinom(20000, 1, 0.1)              # hypothesis true or false
#'    Z   <- rnorm(20000, H*X)                  # z-score
#'    .Random.seed <- save.seed
#'    pvalue <- 1-pnorm(Z)                      #pvalue
#'    ihw_res <- ihw(pvalue, X, .1)
#'    plot(ihw_res)
#'
#' @export
#' @aliases plot,ihwResult-method

setMethod("plot", signature="ihwResult",
	function(x, x_axis = c(weights = "group", decisionboundary = "covariate")[what], what = "weights", scale = covariate_type(x)) {
	  if(x@stratification_method == "forest") stop("plotting not available for forest stratification method")
	  
	  ws <- weights(x, levels_only=TRUE)
	  group_levels <- levels(groups_factor(x))
	  folds  <- factor(seq_len(nfolds(x)))
	  
	  df <- expand.grid(stratum = seq_len(nlevels(groups_factor(x))),
	                    fold    = folds)
	  df$group   <- group_levels[df$stratum]
	  df$weight  <- mapply(function(x,y) ws[x,y], df$stratum, df$fold)

	  if (!(requireNamespace("ggplot2", quietly=TRUE) && (requireNamespace("scales", quietly=TRUE)))) 
	    stop("Please install the packages 'ggplot2' and 'scales'.")
	  
	  switch(what,
      weights = plot_weights(df, x_axis, scale),
      decisionboundary = plot_decisionboundary(x, df, x_axis, scale),
      stop(sprintf("Invalid what='%s'.", what))
    )
})

# plot weights vs rank (ordinal)
# plot weights vs filter statistic value
plot_weights <- function(df, x_axis, scale) {

  switch(scale,
     nominal = {
      if (x_axis != "group")
             warning("'x_axis' argument ignored since 'scale' is nominal.")
           
      ggplot2::ggplot(df, ggplot2::aes_string(x = 'fold',y = 'weight', fill = 'fold')) +
             ggplot2::geom_bar(stat = "identity", position="dodge") +
             ggplot2::facet_wrap(~group)
     },
     ordinal = {
      switch(x_axis,
        group = {
          ggplot2::ggplot(df, ggplot2::aes_string(x='stratum', y='weight', col='fold'))+
               ggplot2::geom_point() +
               ggplot2::geom_line() +
               ggplot2::scale_x_continuous(breaks= scales::pretty_breaks())
        }, 
        stop(sprintf("x_axis = '%s' not implemented.", x_axis))
      )},
      stop(sprintf("Invalid scale='%s'.", scale))
  ) 
}  


plot_decisionboundary <- function(x, df, x_axis, scale) {

  if (x_axis != "covariate")
    warning(sprintf("x_axis='%s' not implemented.", x_axis))
 
  xdf <- x@df[ !is.na(x@df$pvalue), ]

  minp <- 10 ^ (-ceiling(log10(nrow(xdf))) - 6)
                             
  xdf$pvalue[ xdf$pvalue < minp] <- minp
  xdf$rkcv <- rank(xdf$covariate, ties.method = "min")
    
  minrk      <- by(xdf, xdf$group, function(v) min(v$rkcv))
  pthresh    <- by(xdf, list(xdf$group, xdf$fold), function(v) {
                        wh <- (v$adj_pvalue <= x@alpha)
                        if (!any(wh)) minp else max(v$pvalue[wh], na.rm = TRUE) 
                        })
  
  df$minrk   <- minrk[as.character(df$group)]
  df$pthresh <- pthresh[cbind(as.character(df$group), as.character(df$fold))]

  # add another bin for plotting
  toadd <- subset(df, stratum == nlevels(groups_factor(x)))
  toadd$minrk   <- rep(nrow(xdf), nrow(toadd))  
  df <- rbind(df, toadd)
    
  switch(scale,
         nominal = {
           stop("Not yet implemented.")
         },
         ordinal = {
             ggplot2::ggplot(xdf, ggplot2::aes(x = rank(covariate), y = -log10(pvalue))) +
               ggplot2::geom_hex( bins = 100) +
               ggplot2::geom_step( ggplot2::aes(x = minrk, y = -log10(pthresh), col = fold), data = df, direction = "hv")
         },
         stop(sprintf("Invalid scale='%s'.", scale))
  ) 
}  
