---
title: Charges and Missing Hydrogens
keywords: 
last_updated: Thu Jul 28 04:35:19 2016
sidebar: mydoc_sidebar
permalink: /mydoc/ChemmineR_14/
---

The function `bonds` returns information about the number
of bonds, charges and missing hydrogens in `SDF` and
`SDFset` objects. It is used by many other functions
(*e.g.* `MW`, `MF`,
`atomcount`, `atomcuntMA` and
`plot`) to correct for missing hydrogens that are often
not specified in SD files. 

```r
 bonds(sdfset[[1]], type="bonds")[1:4,]
```

```
##   atom Nbondcount Nbondrule charge
## 1    O          2         2      0
## 2    O          2         2      0
## 3    O          2         2      0
## 4    O          2         2      0
```

```r
 bonds(sdfset[1:2], type="charge")
```

```
## $CMP1
## NULL
## 
## $CMP2
## NULL
```

```r
 bonds(sdfset[1:2], type="addNH") 
```

```
## CMP1 CMP2 
##    0    0
```


<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/ChemmineR_15/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
