
# Independent Hypothesis Weighting

Independent hypothesis weighting (IHWold) is a multiple testing procedure that increases power compared to the method of Benjamini and Hochberg by assigning
data-driven weights to each hypothesis. The input to IHWold is a two-column
table of p-values and covariates. The covariate can be any continuous-valued
or categorical variable that is thought to be informative on the statistical
properties of each hypothesis test, while it is independent of the p-value
under the null hypothesis. IHWold is described in the following paper:

> N. Ignatiadis, B. Klaus, J.B. Zaugg, W. Huber. *Data-driven hypothesis weighting increases detection power in genome-scale multiple testing.* Nature methods. 2016 Jul;13(7):577-80.

Also see the following paper for the theoretical underpinning of the method:

> N. Ignatiadis and W. Huber. *Covariate-powered cross weighted multiple testing.*  [[arXiv]](https://arxiv.org/abs/1701.05179)


# Software availability

The package is available on  [Bioconductor](https://www.bioconductor.org/packages/devel/bioc/html/IHWold.html), and may be installed as follows:

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("IHWold")
```

The package can be installed as follows with `devtools` from the Github repository:

```R
devtools::install_github("nignatiadis/IHWold")
```




