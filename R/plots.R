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
	  
	  #multi-dimensional covariate handling
	  nvar <- sum(grepl("^covariate", names(x@df)))
	  if(nvar == 2 & x@stratification_method == "quantiles" & x_axis == "covariate"& what == "weights" & scale == "ordinal"){
	    return(plot_weights_quantile_2d(x))
	  }else if(nvar > 2){
	    stop(sprintf("Plotting not implemented for %s covariates", nvar))
	  }
	  
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

# multivariate quantile binning 
groups_by_filter_multivariate_eval <- function(covariates_train, covariates_eval, nbins) {
  nvar <- ncol(covariates_train)
  nbin_dim <- max(1, floor(nbins ^ (1 / nvar))) #does not change anything, if nvar = 1
  
  groups <- lapply(seq_len(nvar), function(i) {
    covariate_train_i <- covariates_train[, i, drop = TRUE]
    covariates_eval_i <- covariates_eval[, i, drop = TRUE]
    
    breaks <- quantile(covariate_train_i, probs = seq(0, 1, length.out = nbin_dim+1), na.rm = TRUE)
    #add interval open to the right
    
    covariates_eval_i_binned <- cut(covariates_eval_i, breaks, include.lowest = TRUE, right = FALSE)
  })
  
  if(nvar == 1){
    groups <- unname(unlist(groups))
  }else{
    groups <- do.call(cbind, groups)
    groups <- apply(groups, 1, paste, collapse = "-") # base R equivalent of tidyr::unite
  }
  # convert to factor
  groups <- as.factor(groups)
  
  groups
}

plot_weights_quantile_2d <- function(ihw_quantiles, m_eval = 100,  fold = 1){
  ihw_quantiles_df <- ihw_quantiles@df
  
  #build grid to evaluate for plotting
  data_eval <- expand.grid(
    covariate.1 = seq(min(ihw_quantiles_df$covariate.1), max(ihw_quantiles_df$covariate.1), length.out = m_eval),
    covariate.2 = seq(min(ihw_quantiles_df$covariate.2), max(ihw_quantiles_df$covariate.2), length.out = m_eval)
  )
  
  covariates_train <- ihw_quantiles_df[, c("covariate.1","covariate.2")] 
  
  data_eval$groups <- groups_by_filter_multivariate_eval(
    covariates_train = covariates_train, 
    covariates_eval = data_eval, 
    nbins = ihw_quantiles@nbins
  )
  
  #relevel group labels
  levels(data_eval$groups) <- seq_len(nlevels(data_eval$groups))
  
  #obtain weights
  data_eval$weight <- sapply(data_eval$groups, function(group_i){
    ihw_quantiles@weights[group_i,fold]
  })
  
  #plot everything
  plot <- data_eval %>%
    ggplot2::ggplot(ggplot2::aes(x = covariate.1, y = covariate.2, color = weight))+
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks= scales::pretty_breaks())
  
  plot
}