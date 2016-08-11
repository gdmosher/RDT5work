---
title: Getting Around
keywords: 
last_updated: Thu Jul 28 05:13:47 2016
sidebar: mydoc_sidebar
permalink: /mydoc/Rbasics_04/
---

## Startup and Closing Behavior

* __Starting R__:
    The R GUI versions, including RStudio, under Windows and Mac OS X can be
    opened by double-clicking their icons. Alternatively, one can start it by
    typing `R` in a terminal (default under Linux). 

* __Startup/Closing Behavior__:
    The R environment is controlled by hidden files in the startup directory:
    `.RData`, `.Rhistory` and `.Rprofile` (optional). 
	
    
* __Closing R__:


```r
q()  
```
```
Save workspace image? [y/n/c]:
```
        
* __Note__:
    When responding with `y`, then the entire R workspace will be written to
    the `.RData` file which can become very large. Often it is sufficient to just
    save an analysis protocol in an R source file. This way one can quickly
    regenerate all data sets and objects. 


## Navigating directories

Create an object with the assignment operator `<-` or `=`

```r
object <- ...
```

List objects in current R session

```r
ls()
```

Return content of current working directory

```r
dir()
```

Return path of current working directory

```r
getwd()
```

Change current working directory

```r
setwd("/home/user")
```

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/Rbasics_05/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
