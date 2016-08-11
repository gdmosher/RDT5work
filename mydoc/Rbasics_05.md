---
title: Basic Syntax
keywords: 
last_updated: Thu Jul 28 05:13:47 2016
sidebar: mydoc_sidebar
permalink: /mydoc/Rbasics_05/
---

General R command syntax


```r
object <- function_name(arguments) 
object <- object[arguments] 
```

Finding help


```r
?function_name
```

Load a library/package


```r
library("my_library") 
```

List functions defined by a library


```r
library(help="my_library")
```

Load library manual (PDF or HTML file)


```r
vignette("my_library") 
```

Execute an R script from within R


```r
source("my_script.R")
```

Execute an R script from command-line (the first of the three options is preferred)


```sh
$ Rscript my_script.R
$ R CMD BATCH my_script.R 
$ R --slave < my_script.R 
```

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/Rbasics_06/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
