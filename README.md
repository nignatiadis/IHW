[![Travis-CI Build Status](https://travis-ci.org/nignatiadis/IHW.svg?branch=master)](https://travis-ci.org/nignatiadis/IHW)
[![codecov.io](https://codecov.io/github/nignatiadis/IHW/coverage.svg?branch=master)](https://codecov.io/github/nignatiadis/IHW?branch=master)

# Independent Hypothesis Weighting
An R package which implements the method described in:

[Data-driven hypothesis weighting increases detection power in genome-scale multiple testing](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3885.html)

You can also find our preprint under:
[Data-driven hypothesis weighting increases detection power in multiple testing](http://biorxiv.org/content/early/2016/03/30/034330)



# Software availability

The package can be installed as follows with `devtools` from the Github repositories:

```R
devtools::install_github("vladchimescu/lpsymphony", subdir="lpsymphony")
devtools::install_github("nignatiadis/IHW")
```

In addition, the package is available in the current release of Bioconductor ([IHW](http://www.bioconductor.org/packages/IHW)). Thus IHW can also be installed from there as follows:

```R
source("https://bioconductor.org/biocLite.R")
biocLite("IHW")
```
