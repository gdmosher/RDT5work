---
title: Summary statistics of variants
keywords: 
last_updated: Thu Jul 28 05:13:49 2016
sidebar: mydoc_sidebar
permalink: /mydoc/systemPipeVARseq_09/
---

The `varSummary` function counts the number of variants for each feature type
included in the anntation reports.

## Summary for `GATK`


```r
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_gatk_filtered.txt")
varSummary(args)
write.table(varSummary(args), "./results/variantStats_gatk.xls", quote=FALSE, col.names = NA, sep="\t")
```

## Summary for `BCFtools`


```r
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_sambcf_filtered.txt")
varSummary(args)
write.table(varSummary(args), "./results/variantStats_sambcf.xls", quote=FALSE, col.names = NA, sep="\t")
```

## Summary for `VariantTools`  


```r
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
varSummary(args)
write.table(varSummary(args), "./results/variantStats_vartools.xls", quote=FALSE, col.names = NA, sep="\t")
```

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/systemPipeVARseq_10/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
