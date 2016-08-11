---
title: Useful R Functions
keywords: 
last_updated: Thu Jul 28 05:13:47 2016
sidebar: mydoc_sidebar
permalink: /mydoc/Rbasics_11/
---

## Unique entries

Make vector entries unique with `unique`


```r
length(iris$Sepal.Length)
```

```
## [1] 150
```

```r
length(unique(iris$Sepal.Length))
```

```
## [1] 35
```

## Count occurrences

Count occurrences of entries with `table`

```r
table(iris$Species)
```

```
## 
##     setosa versicolor  virginica 
##         50         50         50
```

## Aggregate data

Compute aggregate statistics with `aggregate`

```r
aggregate(iris[,1:4], by=list(iris$Species), FUN=mean, na.rm=TRUE)
```

```
##      Group.1 Sepal.Length Sepal.Width Petal.Length Petal.Width
## 1     setosa        5.006       3.428        1.462       0.246
## 2 versicolor        5.936       2.770        4.260       1.326
## 3  virginica        6.588       2.974        5.552       2.026
```

## Intersect data

Compute intersect between two vectors with `%in%`

```r
month.name %in% c("May", "July")
```

```
##  [1] FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE
```

## Merge data frames

Join two data frames by common field entries with `merge` (here row names `by.x=0`). To obtain only the common rows, change `all=TRUE` to `all=FALSE`. To merge on specific columns, refer to them by their position numbers or their column names.

```r
frame1 <- iris[sample(1:length(iris[,1]), 30), ]
frame1[1:2,]
```

```
##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 3          4.7         3.2          1.3         0.2  setosa
## 5          5.0         3.6          1.4         0.2  setosa
```

```r
dim(frame1)
```

```
## [1] 30  5
```

```r
my_result <- merge(frame1, iris, by.x = 0, by.y = 0, all = TRUE)
dim(my_result)
```

```
## [1] 150  11
```

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/Rbasics_12/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
