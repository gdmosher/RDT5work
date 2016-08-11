---
title: Splitting SD Files
keywords: 
last_updated: Thu Jul 28 04:35:19 2016
sidebar: mydoc_sidebar
permalink: /mydoc/ChemmineR_08/
---

The following `write.SDFsplit` function allows to split
SD Files into any number of smaller SD Files. This can become important
when working with very big SD Files. Users should note that this
function can output many files, thus one should run it in a dedicated
directory!  

Create sample SD File with 100 molecules: 

```r
 write.SDF(sdfset, "test.sdf") 
```


Read in sample SD File. Note: reading file into SDFstr is much faster
than into SDFset: 

```r
 sdfstr <- read.SDFstr("test.sdf") 
```


Run export on `SDFstr` object: 

```r
 write.SDFsplit(x=sdfstr, filetag="myfile", nmol=10) # 'nmol' defines the number of molecules to write to each file 
```


Run export on `SDFset` object: 

```r
 write.SDFsplit(x=sdfset, filetag="myfile", nmol=10) 
```


<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/ChemmineR_09/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
