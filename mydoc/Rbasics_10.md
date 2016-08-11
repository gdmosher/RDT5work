---
title: Reading and Writing External Data
keywords: 
last_updated: Thu Jul 28 05:13:47 2016
sidebar: mydoc_sidebar
permalink: /mydoc/Rbasics_10/
---
## Import of tabular data

Import of a tab-delimited tabular file

```r
myDF <- read.delim("myData.xls", sep="\t")
```

Import of Excel file. Note: working with tab- or comma-delimited files is more flexible and preferred.

```r
library(gdata)
myDF <- read.xls"myData.xls")
```

Import of Google Sheets. The following example imports a sample Google Sheet from [here](https://docs.google.com/spreadsheets/d/1U-32UcwZP1k3saKeaH1mbvEAOfZRdNHNkWK2GI1rpPM/edit#gid=472150521).
Detailed instructions for interacting from R with Google Sheets with the required `googlesheets` package are [here](https://github.com/jennybc/googlesheets).


```r
library("googlesheets"); library("dplyr"); library(knitr)
gs_auth() # Creates authorizaton token (.httr-oauth) in current directory if not present
sheetid <-"1U-32UcwZP1k3saKeaH1mbvEAOfZRdNHNkWK2GI1rpPM"
gap <- gs_key(sheetid)
mysheet <- gs_read(gap, skip=4)
myDF <- as.data.frame(mysheet)
myDF
```

## Export of tabular data

```r
write.table(myDF, file="myfile.xls", sep="\t", quote=FALSE, col.names=NA)
```

## Line-wise import

```r
myDF <- readLines("myData.txt")
```

## Line-wise export

```r
writeLines(month.name, "myData.txt")
```

## Copy and paste into R

On Windows/Linux systems

```r
read.delim("clipboard") 
```
On Mac OS X systems

```r
read.delim(pipe("pbpaste")) 
```

## Copy and paste from R 

On Windows/Linux systems

```r
write.table(iris, "clipboard", sep="\t", col.names=NA, quote=F) 
```

On Mac OS X systems

```r
zz <- pipe('pbcopy', 'w')
write.table(iris, zz, sep="\t", col.names=NA, quote=F)
close(zz) 
```

## Homework 3A 

Homework 3A: [Object Subsetting Routines and Import/Export](http://girke.bioinformatics.ucr.edu/GEN242/mydoc/mydoc_homework_03.html)

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/Rbasics_11/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
