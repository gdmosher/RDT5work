---
title: Control Structures
keywords: 
last_updated: Thu Jul 28 05:13:48 2016
sidebar: mydoc_sidebar
permalink: /mydoc/Programming_in_R_02/
---

## Important Operators

### Comparison operators

* `==` (equal)
* `!=` (not equal)
* `>` (greater than)
* `>=` (greater than or equal)
* `<` (less than)
* `<=` (less than or equal)

### Logical operators
		
* `&` (and)
* `|` (or) 
* `!` (not)

## Conditional Executions: `if` Statements

An `if` statement operates on length-one logical vectors.

__Syntax__

```r
if(TRUE) { 
	statements_1 
} else { 
	statements_2 
}
```

__Example__

```r
if(1==0) { 
	print(1) 
} else { 
	print(2) 
}
```

```
## [1] 2
```

## Conditional Executions: `ifelse` Statements

The `ifelse` statement operates on vectors.

__Syntax__

```r
ifelse(test, true_value, false_value)
```
__Example__

```r
x <- 1:10 
ifelse(x<5, x, 0)
```

```
##  [1] 1 2 3 4 0 0 0 0 0 0
```

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/Programming_in_R_03/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
