---
title: Format Interconversions
keywords: 
last_updated: Thu Jul 28 04:35:19 2016
sidebar: mydoc_sidebar
permalink: /mydoc/ChemmineR_07/
---

The `sdf2smiles` and `smiles2sdf`
functions provide format interconversion between SMILES strings
(Simplified Molecular Input Line Entry Specification) and
`SDFset` containers.  

Convert an `SDFset` container to a SMILES
`character` string:



```r
 data(sdfsample);
 sdfset <- sdfsample[1] 
 smiles <- sdf2smiles(sdfset) 
 smiles 
```


Convert a SMILES `character` string to an
`SDFset` container:



```r
 sdf <- smiles2sdf("CC(=O)OC1=CC=CC=C1C(=O)O")
 view(sdf) 
```


When the `ChemineOB` package is installed these
conversions are performed with the OpenBabel Open Source Chemistry
Toolbox. Otherwise the functions will fall back to using the ChemMine
Tools web service for this operation. The latter will require internet
connectivity and is limited to only the first compound given.
`ChemmineOB` provides access to the compound format
conversion functions of OpenBabel. Currently, over 160 formats are
supported by OpenBabel. The functions `convertFormat` and
`convertFormatFile` can be used to convert files or
strings between any two formats supported by OpenBabel. For example, to
convert a SMILES string to an SDF string, one can use the
`convertFormat` function.



```r
 sdfStr <- convertFormat("SMI","SDF","CC(=O)OC1=CC=CC=C1C(=O)O_name") 
```


This will return the given compound as an SDF formatted string. 2D
coordinates are also computed and included in the resulting SDF string.

To convert a file with compounds encoded in one format to another
format, the `convertFormatFile` function can be used
instead. 

```r
 convertFormatFile("SMI","SDF","test.smiles","test.sdf") 
```


To see the whole list of file formats supported by OpenBabel, one can
run from the command-line "obabel -L formats".


<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/ChemmineR_08/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
