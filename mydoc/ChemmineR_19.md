---
title: Searching PubChem
keywords: 
last_updated: Thu Jul 28 04:35:19 2016
sidebar: mydoc_sidebar
permalink: /mydoc/ChemmineR_19/
---

## Get Compounds from PubChem by Id

The function `getIds` accepts one or more numeric PubChem
compound ids and downloads the corresponding compounds from PubChem
Power User Gateway (PUG) returning results in an `SDFset`
container. The ChemMine Tools web service is used as an intermediate, to
translate queries from plain HTTP POST to a PUG SOAP query.  

Fetch 2 compounds from PubChem:



```r
 compounds <- getIds(c(111,123))
 compounds 
```


## Search a SMILES Query in PubChem

The function `searchString` accepts one SMILES string
(Simplified Molecular Input Line Entry Specification) and performs a
\>0.95 similarity PubChem fingerprint search, returning the hits in an
`SDFset` container. The ChemMine Tools web service is
used as an intermediate, to translate queries from plain HTTP POST to a
PubChem Power User Gateway (PUG) query.  

Search a SMILES string on PubChem:



```r
 compounds <- searchString("CC(=O)OC1=CC=CC=C1C(=O)O") compounds 
```


## Search an SDF Query in PubChem

The function `searchSim` performs a PubChem similarity
search just like `searchString`, but accepts a query in
an `SDFset` container. If the query contains more than
one compound, only the first is searched.  

Search an `SDFset` container on PubChem:



```r
 data(sdfsample); 
 sdfset <- sdfsample[1] 
 compounds <- searchSim(sdfset) 
 compounds 
```


<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/ChemmineR_20/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
