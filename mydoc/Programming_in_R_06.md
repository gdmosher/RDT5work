---
title: Running R Scripts
keywords: 
last_updated: Thu Jul 28 05:13:48 2016
sidebar: mydoc_sidebar
permalink: /mydoc/Programming_in_R_06/
---

## Possibilities for Executing R Scripts

### R console

```r
source("my_script.R")
```

### Command-line


```sh
Rscript my_script.R # or just ./myscript.R after making it executable
R CMD BATCH my_script.R # Alternative way 1 
R --slave < my_script.R # Alternative way 2
```
### Passing arguments from command-line to R

Create an R script named `test.R` with the following content:


```sh
myarg <- commandArgs()
print(iris[1:myarg[6], ])
```

Then run it from the command-line like this:

```sh
Rscript test.R 10
```

In the given example the number `10` is passed on from the command-line as an argument to the R script which is used to return to `STDOUT` the first 10 rows of the `iris` sample data. If several arguments are provided, they will be interpreted as one string and need to be split in R with the strsplit function. A more detailed example can be found [here](http://manuals.bioinformatics.ucr.edu/home/ht-seq\#TOC-Quality-Reports-of-FASTQ-Files-).

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/Programming_in_R_07/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
