---
title: Utilities for coverage data
keywords: 
last_updated: Thu Jul 28 05:13:49 2016
sidebar: mydoc_sidebar
permalink: /mydoc/systemPipeChIPseq_05/
---

The following introduces several utilities useful for ChIP-Seq data. They are not part of the actual
workflow.

## Rle object stores coverage information

```r
library(rtracklayer); library(GenomicRanges); library(Rsamtools); library(GenomicAlignments)
aligns <- readGAlignments(outpaths(args)[1])
cov <- coverage(aligns)
cov
```

## Resizing aligned reads

```r
trim(resize(as(aligns, "GRanges"), width = 200))
```

## Naive peak calling

```r
islands <- slice(cov, lower = 15)
islands[[1]]
```

## Plot coverage for defined region

```r
library(ggbio)
myloc <- c("Chr1", 1, 100000)
ga <- readGAlignments(outpaths(args)[1], use.names=TRUE, param=ScanBamParam(which=GRanges(myloc[1], IRanges(as.numeric(myloc[2]), as.numeric(myloc[3])))))
autoplot(ga, aes(color = strand, fill = strand), facets = strand ~ seqnames, stat = "coverage")
```

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/systemPipeChIPseq_06/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
