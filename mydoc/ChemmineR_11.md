---
title: Working with SDF/SDFset Classes
keywords: 
last_updated: Thu Jul 28 04:35:19 2016
sidebar: mydoc_sidebar
permalink: /mydoc/ChemmineR_11/
---

Several methods are available to return the different data components of
`SDF/SDFset` containers in batches. The following
examples list the most important ones. To save space their content is
not printed in the manual. 

```r
 view(sdfset[1:4]) # Summary view of several molecules 

 length(sdfset) # Returns number of molecules 
 sdfset[[1]] # Returns single molecule from SDFset as SDF object 

 sdfset[[1]][[2]] # Returns atom block from first compound as matrix

 sdfset[[1]][[2]][1:4,] 
 c(sdfset[1:4], sdfset[5:8]) # Concatenation of several SDFsets 
```


The `grepSDFset` function allows string
matching/searching on the different data components in
`SDFset`. By default the function returns a SDF summary
of the matching entries. Alternatively, an index of the matches can be
returned with the setting `mode="index"`. 

```r
 grepSDFset("650001", sdfset, field="datablock", mode="subset") # To return index, set mode="index") 
```


Utilities to maintain unique compound IDs: 

```r
 sdfid(sdfset[1:4]) # Retrieves CMP IDs from Molecule Name field in header block. 
 cid(sdfset[1:4]) # Retrieves CMP IDs from ID slot in SDFset. 
 unique_ids <- makeUnique(sdfid(sdfset)) # Creates unique IDs by appending a counter to duplicates. 
 cid(sdfset) <- unique_ids # Assigns uniquified IDs to ID slot 
```


Subsetting by character, index and logical vectors: 

```r
 view(sdfset[c("650001", "650012")])
 view(sdfset[4:1])
 mylog <- cid(sdfset)
 view(sdfset[mylog]) 
```


Accessing `SDF/SDFset` components: header, atom, bond and
data blocks: 

```r
 atomblock(sdf); sdf[[2]];
 sdf[["atomblock"]] # All three methods return the same component

 header(sdfset[1:4]) 
 atomblock(sdfset[1:4])
 bondblock(sdfset[1:4]) 
 datablock(sdfset[1:4])  
 header(sdfset[[1]])
 atomblock(sdfset[[1]]) 
 bondblock(sdfset[[1]]) 
 datablock(sdfset[[1]]) 
```


Replacement Methods: 

```r
 sdfset[[1]][[2]][1,1] <- 999 
 atomblock(sdfset)[1] <- atomblock(sdfset)[2] 
 datablock(sdfset)[1] <- datablock(sdfset)[2] 
```


Assign matrix data to data block: 

```r
 datablock(sdfset) <- as.matrix(iris[1:100,])
 view(sdfset[1:4]) 
```


Class coercions from `SDFstr` to `list`,
`SDF` and `SDFset`: 

```r
 as(sdfstr[1:2], "list") as(sdfstr[[1]], "SDF")
 as(sdfstr[1:2], "SDFset") 
```


Class coercions from `SDF` to `SDFstr`,
`SDFset`, list with SDF sub-components: 

```r
 sdfcomplist <- as(sdf, "list") sdfcomplist <-
 as(sdfset[1:4], "list"); as(sdfcomplist[[1]], "SDF") sdflist <-
 as(sdfset[1:4], "SDF"); as(sdflist, "SDFset") as(sdfset[[1]], "SDFstr")
 as(sdfset[[1]], "SDFset") 
```


Class coercions from `SDFset` to lists with components
consisting of SDF or sub-components: 

```r
 as(sdfset[1:4], "SDF") as(sdfset[1:4], "list") as(sdfset[1:4], "SDFstr")
```


<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/ChemmineR_12/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
