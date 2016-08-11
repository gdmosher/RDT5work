---
title: Bond Matrices
keywords: 
last_updated: Thu Jul 28 04:35:19 2016
sidebar: mydoc_sidebar
permalink: /mydoc/ChemmineR_13/
---

Bond matrices provide an efficient data structure for many basic
computations on small molecules. The function `conMA`
creates this data structure from `SDF` and
`SDFset` objects. The resulting bond matrix contains the
atom labels in the row/column titles and the bond types in the data
part. The labels are defined as follows: 0 is no connection, 1 is a
single bond, 2 is a double bond and 3 is a triple bond. 

```r
 conMA(sdfset[1:2],
 exclude=c("H")) # Create bond matrix for first two molecules in sdfset

 conMA(sdfset[[1]], exclude=c("H")) # Return bond matrix for first molecule 
 plot(sdfset[1], atomnum = TRUE, noHbonds=FALSE , no_print_atoms = "", atomcex=0.8) # Plot its structure with atom numbering 
 rowSums(conMA(sdfset[[1]], exclude=c("H"))) # Return number of non-H bonds for each atom
```


<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/ChemmineR_14/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
