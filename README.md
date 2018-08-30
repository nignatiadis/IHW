[![Travis-CI Build Status](https://travis-ci.org/nignatiadis/IHW.svg?branch=master)](https://travis-ci.org/nignatiadis/IHW)
[![codecov.io](https://codecov.io/github/nignatiadis/IHW/coverage.svg?branch=master)](https://codecov.io/github/nignatiadis/IHW?branch=master)

# Independent Hypothesis Weighting
An R package which implements the method described in:

[Data-driven hypothesis weighting increases detection power in big data analytics](http://biorxiv.org/content/early/2015/12/13/034330)



# Software availability

The package can be installed as follows with `devtools` from the Github repositories:

```R
devtools::install_github("vladchimescu/lpsymphony", subdir="lpsymphony")
devtools::install_github("nignatiadis/IHW")
```

In addition, these packages will be available in the next release of Bioconductor. Both [lpsymphony](http://bioconductor.org/packages/3.3/bioc/html/lpsymphony.html) and [IHW](https://www.bioconductor.org/packages/devel/bioc/html/IHW.html) are also already available on the development branch of Bioconductor. Thus IHW can also be installed from there as follows:

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("IHW")
```
