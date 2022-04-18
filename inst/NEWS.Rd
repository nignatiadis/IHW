\name{NEWS}
\title{News for Package \pkg{IHW}}

\section{Changes in version 1.24.0, Bioconductor 3.15 Release}{ 
\itemize{
    \item IHW now also allows for covariates of class matrix
    \item two stratification methods are available: classical quantile slicing or a new random forest approach
    \item formulas with multiple covariates on the RHS are now supported 
    \item classical quantile slicing only works for a small number of dimensions, while the new random forest approach extends nicely to higher dimensions
    \item Add \code{ihw_forest}, \code{group_by_forest}, \code{groups_by_filter_multivariate} 
    \item add test_forest
}}