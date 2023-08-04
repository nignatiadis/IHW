# Independent Hypothesis Weighting

<p align="center">
<img src="https://github.com/Bioconductor/BiocStickers/blob/devel/IHW/IHW.png?raw=true"  width="200">
</p>





Independent hypothesis weighting (IHW) is a multiple testing procedure that increases power compared to the method of Benjamini and Hochberg by assigning
data-driven weights to each hypothesis. The input to IHW is a two-column
table of p-values and covariates. The covariate can be any continuous-valued
or categorical variable that is thought to be informative on the statistical
properties of each hypothesis test, while it is independent of the p-value
under the null hypothesis. IHW is described in the following paper:

> Nikolaos Ignatiadis, Bernd Klaus, Judith Zaugg, Wolfgang Huber (2016)
*Data-driven hypothesis weighting increases detection power in genome-scale multiple testing.*
Nature methods, 13(7), 577-580.

The following paper describes the theoretical underpinning of the method:

> Nikolaos Ignatiadis, Wolfgang Huber (2021)
*Covariate powered cross-weighted multiple testing.*
Journal of the Royal Statistical Society: Series B (JRSS-B), 83(4), 720-751.


# Software availability

The package is available on  [Bioconductor](https://www.bioconductor.org/packages/devel/bioc/html/IHW.html), and may be installed as follows:

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("IHW")
```

The package can be installed as follows with `devtools` from the Github repository:

```R
devtools::install_github("nignatiadis/IHW")
```




