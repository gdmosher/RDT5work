---
title: Operators and Calculations
keywords: 
last_updated: Thu Jul 28 05:13:47 2016
sidebar: mydoc_sidebar
permalink: /mydoc/Rbasics_09/
---

## Comparison Operators

Comparison operators: `==`, `!=`, `<`, `>`, `<=`, `>=`

```r
1==1
```

```
## [1] TRUE
```
Logical operators: AND: `&`, OR: `|`, NOT: `!`

```r
x <- 1:10; y <- 10:1
x > y & x > 5
```

```
##  [1] FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
```

## Basic Calculations

To look up math functions, see Function Index [here](http://cran.at.r-project.org/doc/manuals/R-intro.html#Function-and-variable-index)

```r
x + y
```

```
##  [1] 11 11 11 11 11 11 11 11 11 11
```

```r
sum(x)
```

```
## [1] 55
```

```r
mean(x)
```

```
## [1] 5.5
```

```r
apply(iris[1:6,1:3], 1, mean) 
```

```
##        1        2        3        4        5        6 
## 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667
```

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/Rbasics_10/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
