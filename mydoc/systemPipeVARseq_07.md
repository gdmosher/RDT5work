---
title: Annotate filtered variants
keywords: 
last_updated: Thu Jul 28 05:13:49 2016
sidebar: mydoc_sidebar
permalink: /mydoc/systemPipeVARseq_07/
---

The function `variantReport` generates a variant report using
utilities provided by the `VariantAnnotation` package. The report for
each sample is written to a tabular file containing genomic context annotations
(_e.g._ coding or non-coding SNPs, amino acid changes, IDs of affected
genes, etc.) along with confidence statistics for each variant. The parameter
file `annotate_vars.param` defines the paths to the input and output
files which are stored in a new `SYSargs` instance. 

## Basics of annotating variants

Variants overlapping with common annotation features can be identified with `locateVariants`.

```r
library("GenomicFeatures")
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_gatk_filtered.txt")
txdb <- loadDb("./data/tair10.sqlite")
vcf <- readVcf(infile1(args)[1], "A. thaliana")
locateVariants(vcf, txdb, CodingVariants())
```
Synonymous/non-synonymous variants of coding sequences are computed by the predictCoding function for variants overlapping with coding regions.


```r
fa <- FaFile(systemPipeR::reference(args))
predictCoding(vcf, txdb, seqSource=fa)
```

## Annotate filtered variants called by `GATK`


```r
library("GenomicFeatures")
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_gatk_filtered.txt")
txdb <- loadDb("./data/tair10.sqlite")
fa <- FaFile(systemPipeR::reference(args))
suppressAll(variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana"))
```

## Annotate filtered variants called by `BCFtools`


```r
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_sambcf_filtered.txt")
txdb <- loadDb("./data/tair10.sqlite")
fa <- FaFile(systemPipeR::reference(args))
suppressAll(variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana"))
```

## Annotate filtered variants called by `VariantTools`


```r
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
txdb <- loadDb("./data/tair10.sqlite")
fa <- FaFile(systemPipeR::reference(args))
suppressAll(variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana"))
```

View annotation result for single sample

```r
read.delim(outpaths(args)[1])[38:40,]
```

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/systemPipeVARseq_08/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
