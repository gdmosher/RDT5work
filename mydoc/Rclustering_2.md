---
title: Data Preprocessing
keywords: 
last_updated: Thu Jul 28 05:13:48 2016
sidebar: mydoc_sidebar
permalink: /mydoc/Rclustering_2/
---

## Data Transformations

Choice depends on data set!
	
- Center and standardize
    1. Center: subtract from each value the mean of the corresponding vector
	2.  Standardize: devide by standard deviation
	- Result: _Mean = 0_ and _STDEV = 1_

- Center and scale with the `scale()` function
    1. Center: subtract from each value the mean of the corresponding vector
	2. Scale: divide centered vector by their _root mean square_ (_rms_):
    $$ x_{rms} = \sqrt[]{\frac{1}{n-1}\sum_{i=1}^{n}{x_{i}{^2}}} $$
    - Result: _Mean = 0_ and _STDEV = 1_

- Log transformation 
- Rank transformation: replace measured values by ranks 
- No transformation

## Distance Methods

List of most common ones!

- Euclidean distance for two profiles _X_ and _Y_:
  $$ d(X,Y) = \sqrt[]{ \sum_{i=1}^{n}{(x_{i}-y_{i})^2} }$$
    - __Disadvantages__: not scale invariant, not for negative correlations
- Maximum, Manhattan, Canberra, binary, Minowski, ...
- Correlation-based distance: _1-r_
    - Pearson correlation coefficient (PCC):
	  $$r = \frac{n\sum_{i=1}^{n}{x_{i}y_{i}} - \sum_{i=1}^{n}{x_{i}} \sum_{i=1}^{n}{y_{i}}}{ \sqrt[]{(\sum_{i=1}^{n}{x_{i}^2} - (\sum_{i=1}^{n}{x_{i})^2}) (\sum_{i=1}^{n}{y_{i}^2} - (\sum_{i=1}^{n}{y_{i})^2})} }$$
        - __Disadvantage__: outlier sensitive 
	- Spearman correlation coefficient (SCC)
	    - Same calculation as PCC but with ranked values!

There are many more distance measures

- If the distances among items are quantifiable, then clustering is possible.
- Choose the most accurate and meaningful distance measure for a given field of application.
- If uncertain then choose several distance measures and compare the results. 

## Cluster Linkage

<center><img title="cluster linkage" src="../Rclustering_files/linkage.png" ></center>

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/Rclustering_3/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
