---
title: Clustering and heat maps
keywords: 
last_updated: Thu Jul 28 05:13:49 2016
sidebar: mydoc_sidebar
permalink: /mydoc/systemPipeRNAseq_08/
---

The following example performs hierarchical clustering on the `rlog`
transformed expression matrix subsetted by the DEGs identified in the above
differential expression analysis. It uses a Pearson correlation-based distance
measure and complete linkage for cluster joining.


```r
library(pheatmap)
geneids <- unique(as.character(unlist(DEG_list[[1]])))
y <- assay(rlog(dds))[geneids, ]
pdf("heatmap1.pdf")
pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
dev.off()
```

![](../systemPipeRNAseq_files/heatmap1.png)
<div align="center">Figure 6: Heat Map with Hierarchical Clustering Dendrograms of DEGs</div>

<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/systemPipeRNAseq_09/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
