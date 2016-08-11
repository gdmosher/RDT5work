---
title: Installation of R and Add-on Packages
keywords: 
last_updated: Thu Jul 28 05:13:47 2016
sidebar: mydoc_sidebar
permalink: /mydoc/Rbasics_03/
---

(1.) Install R for your operating system from [CRAN](http://cran.at.r-project.org/).

(2.) Install RStudio from [RStudio](http://www.rstudio.com/ide/download).

(3.) Install CRAN Packages from R console like this:


```r
install.packages(c("pkg1", "pkg2")) 
install.packages("pkg.zip", repos=NULL)
```

(4.) Install Bioconductor packages as follows:


```r
source("http://www.bioconductor.org/biocLite.R")
library(BiocInstaller)
BiocVersion()
biocLite()
biocLite(c("pkg1", "pkg2"))
```

(5.) For more details consult the [Bioc Install page](http://www.bioconductor.org/install/)
and [BiocInstaller](http://www.bioconductor.org/packages/release/bioc/html/BiocInstaller.html) package.

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/Rbasics_04/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
