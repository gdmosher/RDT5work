---
title: Export of Compounds
keywords: 
last_updated: Thu Jul 28 04:35:19 2016
sidebar: mydoc_sidebar
permalink: /mydoc/ChemmineR_06/
---

## SDF Export

Write objects of classes `SDFset/SDFstr/SDF` to SD file:


```r
 write.SDF(sdfset[1:4], file="sub.sdf") 
```


Writing customized `SDFset` to file containing
`ChemmineR` signature, IDs from `SDFset`
and no data block: 

```r
 write.SDF(sdfset[1:4], file="sub.sdf", sig=TRUE, cid=TRUE, db=NULL) 
```


Example for injecting a custom matrix/data frame into the data block of
an `SDFset` and then writing it to an SD file:


```r
 props <- data.frame(MF=MF(sdfset), MW=MW(sdfset), atomcountMA(sdfset)) 
 datablock(sdfset) <- props
 write.SDF(sdfset[1:4], file="sub.sdf", sig=TRUE, cid=TRUE) 
```


Indirect export via `SDFstr` object: 

```r
 sdf2str(sdf=sdfset[[1]], sig=TRUE, cid=TRUE) # Uses default components 
 sdf2str(sdf=sdfset[[1]], head=letters[1:4], db=NULL) # Uses custom components for header and data block 
```


Write `SDF`, `SDFset` or
`SDFstr` classes to file: 

```r
 write.SDF(sdfset[1:4], file="sub.sdf", sig=TRUE, cid=TRUE, db=NULL)
 write.SDF(sdfstr[1:4], file="sub.sdf") 
 cat(unlist(as(sdfstr[1:4], "list")), file="sub.sdf", sep="") 
```


## SMILES Export

Write objects of class `SMIset` to SMILES file with and
without compound identifiers: 

```r
 data(smisample); smiset <- smisample # Sample data set 

 write.SMI(smiset[1:4], file="sub.smi", cid=TRUE) write.SMI(smiset[1:4], file="sub.smi", cid=FALSE) 
```


<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/ChemmineR_07/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
