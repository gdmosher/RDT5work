---
title: Programming Exercises
keywords: 
last_updated: Thu Jul 28 05:13:48 2016
sidebar: mydoc_sidebar
permalink: /mydoc/Programming_in_R_08/
---

## Exercise 1

### `for` loop

__Task 1.1__: Compute the mean of each row in `myMA` by applying the mean function in a `for` loop.


```r
myMA <- matrix(rnorm(500), 100, 5, dimnames=list(1:100, paste("C", 1:5, sep="")))
myve_for <- NULL
for(i in seq(along=myMA[,1])) {
	myve_for <- c(myve_for, mean(as.numeric(myMA[i, ])))
}
myResult <- cbind(myMA, mean_for=myve_for)
myResult[1:4, ]
```

```
##           C1         C2         C3          C4         C5   mean_for
## 1 -0.9832766  0.8446066  0.4196481  0.23814667  0.3493797  0.1737009
## 2  1.2980835  0.6924483  0.5996748 -0.51642965  0.6701380  0.5487830
## 3 -1.1466949  0.6752775 -0.9384848  0.07464206 -0.7651351 -0.4200791
## 4  0.8122978 -1.3107710  0.6664631  0.12316103 -0.1726270  0.0237048
```

### `while` loop

__Task 1.2__: Compute the mean of each row in `myMA` by applying the mean function in a `while` loop.


```r
z <- 1
myve_while <- NULL
while(z <= length(myMA[,1])) {
	myve_while <- c(myve_while, mean(as.numeric(myMA[z, ])))
	z <- z + 1
}
myResult <- cbind(myMA, mean_for=myve_for, mean_while=myve_while)
myResult[1:4, -c(1,2)]
```

```
##           C3          C4         C5   mean_for mean_while
## 1  0.4196481  0.23814667  0.3493797  0.1737009  0.1737009
## 2  0.5996748 -0.51642965  0.6701380  0.5487830  0.5487830
## 3 -0.9384848  0.07464206 -0.7651351 -0.4200791 -0.4200791
## 4  0.6664631  0.12316103 -0.1726270  0.0237048  0.0237048
```
__Task 1.3__: Confirm that the results from both mean calculations are identical

```r
all(myResult[,6] == myResult[,7])
```

```
## [1] TRUE
```

### `apply` loop
	
__Task 1.4__: Compute the mean of each row in myMA by applying the mean function in an `apply` loop

```r
myve_apply <- apply(myMA, 1, mean)
myResult <- cbind(myMA, mean_for=myve_for, mean_while=myve_while, mean_apply=myve_apply)
myResult[1:4, -c(1,2)]
```

```
##           C3          C4         C5   mean_for mean_while mean_apply
## 1  0.4196481  0.23814667  0.3493797  0.1737009  0.1737009  0.1737009
## 2  0.5996748 -0.51642965  0.6701380  0.5487830  0.5487830  0.5487830
## 3 -0.9384848  0.07464206 -0.7651351 -0.4200791 -0.4200791 -0.4200791
## 4  0.6664631  0.12316103 -0.1726270  0.0237048  0.0237048  0.0237048
```

### Avoiding loops

__Task 1.5__: When operating on large data sets it is much faster to use the rowMeans function


```r
mymean <- rowMeans(myMA)
myResult <- cbind(myMA, mean_for=myve_for, mean_while=myve_while, mean_apply=myve_apply, mean_int=mymean)
myResult[1:4, -c(1,2,3)]
```

```
##            C4         C5   mean_for mean_while mean_apply   mean_int
## 1  0.23814667  0.3493797  0.1737009  0.1737009  0.1737009  0.1737009
## 2 -0.51642965  0.6701380  0.5487830  0.5487830  0.5487830  0.5487830
## 3  0.07464206 -0.7651351 -0.4200791 -0.4200791 -0.4200791 -0.4200791
## 4  0.12316103 -0.1726270  0.0237048  0.0237048  0.0237048  0.0237048
```

## Exercise 2 

### Custom functions

__Task 2.1__: Use the following code as basis to implement a function that allows the user to compute the mean for any combination of columns in a matrix or data frame. The first argument of this function should specify the input data set, the second the mathematical function to be passed on (_e.g._ `mean`, `sd`, `max`) and the third one should allow the selection of the columns by providing a grouping vector.


```r
myMA <- matrix(rnorm(100000), 10000, 10, dimnames=list(1:10000, paste("C", 1:10, sep="")))
myMA[1:2,]
```

```
##           C1         C2         C3        C4         C5        C6         C7         C8         C9
## 1 0.17958477 -0.5714262 -0.8866647  2.463907 -0.8126814 -0.933438 -0.8118745 -0.5546591 -0.9026947
## 2 0.04191538  1.2456075  0.4953736 -2.703100  0.2055566 -1.945671 -2.8696620 -1.5043589 -1.8783809
##          C10
## 1 -0.5698361
## 2 -1.6321152
```

```r
myList <- tapply(colnames(myMA), c(1,1,1,2,2,2,3,3,4,4), list) 
names(myList) <- sapply(myList, paste, collapse="_")
myMAmean <- sapply(myList, function(x) apply(myMA[,x], 1, mean))
myMAmean[1:4,] 
```

```
##     C1_C2_C3   C4_C5_C6      C7_C8     C9_C10
## 1 -0.4261687  0.2392626 -0.6832668 -0.7362654
## 2  0.5942988 -1.4810715 -2.1870104 -1.7552480
## 3  0.1217488 -0.7225502 -0.6295343  0.4990018
## 4 -0.9118941 -0.3107419  0.3284317 -0.5693107
```
<!---
Solution

-->


## Exercise 3

### Nested loops to generate similarity matrices

__Task 3.1__: Create a sample list populated with character vectors of different lengths


```r
setlist <- lapply(11:30, function(x) sample(letters, x, replace=TRUE))
names(setlist) <- paste("S", seq(along=setlist), sep="") 
setlist[1:6]
```

```
## $S1
##  [1] "x" "r" "j" "n" "l" "z" "b" "o" "v" "j" "i"
## 
## $S2
##  [1] "k" "b" "p" "c" "z" "f" "v" "u" "e" "d" "c" "f"
## 
## $S3
##  [1] "l" "e" "p" "j" "i" "k" "y" "i" "w" "l" "w" "x" "p"
## 
## $S4
##  [1] "d" "e" "v" "o" "h" "q" "i" "e" "d" "y" "o" "m" "q" "y"
## 
## $S5
##  [1] "s" "q" "r" "j" "o" "z" "q" "g" "s" "v" "w" "j" "l" "r" "d"
## 
## $S6
##  [1] "c" "l" "h" "v" "e" "a" "i" "u" "g" "h" "s" "f" "u" "b" "e" "y"
```

__Task 3.2__: Compute the length for all pairwise intersects of the vectors stored in `setlist`. The intersects can be determined with the `%in%` function like this: `sum(setlist[[1]] %in% setlist[[2]])`


```r
setlist <- sapply(setlist, unique)
olMA <- sapply(names(setlist), function(x) sapply(names(setlist), 
               function(y) sum(setlist[[x]] %in% setlist[[y]])))
olMA[1:12,] 
```

```
##     S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19 S20
## S1  10  3  4  3  6  4  6  6  7   8   7   2   5   8   8   8   8   7   8   7
## S2   3 10  3  3  3  6  7  4  7   2   5   6   4   8   7   6   8   6   7   7
## S3   4  3  9  3  3  4  4  5  4   5   5   3   4   5   7   7   8   5   5   5
## S4   3  3  3  9  4  5  3  4  5   5   5   7   7   5   8   7   7   5   6   5
## S5   6  3  3  4 11  4  3  5  6   6   7   7   8  10   9   8   8   7   8   8
## S6   4  6  4  5  4 13  4  5  7   5   7   8   6   9   9   8   9   7   7   9
## S7   6  7  4  3  3  4 12  6  7   5   9   4   6   9   9   8  10   8  10   8
## S8   6  4  5  4  5  5  6 13  9   9  10   7   7   9  12   9  11   9   9   7
## S9   7  7  4  5  6  7  7  9 14   7   8   7   7  11  11  10  11   8   9  11
## S10  8  2  5  5  6  5  5  9  7  13   8   6   6   6  10  11   9   9   8   9
## S11  7  5  5  5  7  7  9 10  8   8  15   7   9  11  11  10  12  11  12  10
## S12  2  6  3  7  7  8  4  7  7   6   7  14   8  10  11   9  11   9   9   8
```
__Task 3.3__ Plot the resulting intersect matrix as heat map. The `image` or the `heatmap.2` function from the `gplots` library can be used for this.

```r
image(olMA)
```

![](../Programming_in_R_files/nested_loops3-1.png)

## Exercise 4

### Build your own R package

__Task 4.1__: Save one or more of your functions to a file called `script.R` and build the package with the `package.skeleton` function.


```r
package.skeleton(name="mypackage", code_files=c("script1.R"), namespace=TRUE)
```

__Task 4.2__: Build tarball of the package


```r
system("R CMD build mypackage")
```

__Task 4.3__: Install and use package


```r
install.packages("mypackage_1.0.tar.gz", repos=NULL, type="source")
library(mypackage)
?myMAcomp # Opens help for function defined by mypackage
```

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/Programming_in_R_09/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
