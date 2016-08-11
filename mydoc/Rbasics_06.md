---
title: Data Types 
keywords: 
last_updated: Thu Jul 28 05:13:47 2016
sidebar: mydoc_sidebar
permalink: /mydoc/Rbasics_06/
---

## Numeric data

Example: `1, 2, 3, ...`


```r
x <- c(1, 2, 3)
x
```

```
## [1] 1 2 3
```

```r
is.numeric(x)
```

```
## [1] TRUE
```

```r
as.character(x)
```

```
## [1] "1" "2" "3"
```

## Character data

Example: `"a", "b", "c", ...`


```r
x <- c("1", "2", "3")
x
```

```
## [1] "1" "2" "3"
```

```r
is.character(x)
```

```
## [1] TRUE
```

```r
as.numeric(x)
```

```
## [1] 1 2 3
```

## Complex data

Example: mix of both


```r
c(1, "b", 3)
```

```
## [1] "1" "b" "3"
```

## Logical data

Example: `TRUE` of `FALSE`


```r
x <- 1:10 < 5
x  
```

```
##  [1]  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
```

```r
!x
```

```
##  [1] FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
```

```r
which(x) # Returns index for the 'TRUE' values in logical vector
```

```
## [1] 1 2 3 4
```

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/Rbasics_07/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
