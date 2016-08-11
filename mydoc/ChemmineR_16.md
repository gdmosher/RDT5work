---
title: Rendering Chemical Structure Images
keywords: 
last_updated: Thu Jul 28 04:35:19 2016
sidebar: mydoc_sidebar
permalink: /mydoc/ChemmineR_16/
---

## R Graphics Device

A new plotting function for compound structures has been added to the
package recently. This function uses the native R graphics device for
generating compound depictions. At this point this function is still in
an experimental developmental stage but should become stable soon.  

If you have `ChemmineOB` available you can use the `regenCoords`
option to have OpenBabel regenerate the coordinates for the compound.
This can sometimes produce better looking plots.

Plot compound Structures with R's graphics device: 

```r
 data(sdfsample)
 sdfset <- sdfsample
 plot(sdfset[1:4], regenCoords=TRUE,print=FALSE) # 'print=TRUE' returns SDF summaries
```

![](../ChemmineR_files/plotstruct2-1.png)


Customized plots: 

```r
 plot(sdfset[1:4], griddim=c(2,2), print_cid=letters[1:4], print=FALSE, 
		noHbonds=FALSE) 
```


In the following plot, the atom block position numbers in the SDF are
printed next to the atom symbols (`atomnum = TRUE`). For
more details, consult help documentation with
`?plotStruc` or `?plot`. 

```r
 plot(sdfset["CMP1"], atomnum = TRUE, noHbonds=F , no_print_atoms = "",
	  	atomcex=0.8, sub=paste("MW:", MW(sdfsample["CMP1"])), print=FALSE) 
```

![](../ChemmineR_files/plotstruct3-1.png)


Substructure highlighting by atom numbers: 

```r
 plot(sdfset[1], print=FALSE, colbonds=c(22,26,25,3,28,27,2,23,21,18,8,19,20,24)) 
```

![](../ChemmineR_files/plotstruct4-1.png)


## Online with ChemMine Tools

Alternatively, one can visualize compound structures with a standard web
browser using the online ChemMine Tools service.

Plot structures using web service ChemMine Tools: 

```r
 sdf.visualize(sdfset[1:4]) 
```

![Figure: Visualization webpage created by calling `sdf.visualize`.](../ChemmineR_files/visualizescreenshot-small.png)


<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/ChemmineR_17/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
