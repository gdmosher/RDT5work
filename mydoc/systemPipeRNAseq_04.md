---
title: Alignments
keywords: 
last_updated: Thu Jul 28 05:13:49 2016
sidebar: mydoc_sidebar
permalink: /mydoc/systemPipeRNAseq_04/
---

## Read mapping with `Bowtie2/Tophat2` 

The NGS reads of this project will be aligned against the reference
genome sequence using `Bowtie2/TopHat2` (Kim et al., 2013; Langmead et al., 2012). The parameter 
settings of the aligner are defined in the `tophat.param` file.


```r
args <- systemArgs(sysma="param/tophat.param", mytargets="targets.txt")
sysargs(args)[1] # Command-line parameters for first FASTQ file
```


Submission of alignment jobs to compute cluster, here using 72 CPU cores
(18 `qsub` processes each with 4 CPU cores).


```r
moduleload(modules(args))
system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta")
resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="10gb")
reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01", 
                  resourceList=resources)
waitForJobs(reg)
```

Check whether all BAM files have been created


```r
file.exists(outpaths(args))
```

## Read and alignment stats

The following provides an overview of the number of reads in each sample
and how many of them aligned to the reference.


```r
read_statsDF <- alignStats(args=args) 
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
```


```r
read.table(system.file("extdata", "alignStats.xls", package="systemPipeR"), header=TRUE)[1:4,]
```

```
##   FileName Nreads2x Nalign Perc_Aligned Nalign_Primary Perc_Aligned_Primary
## 1      M1A   192918 177961     92.24697         177961             92.24697
## 2      M1B   197484 159378     80.70426         159378             80.70426
## 3      A1A   189870 176055     92.72397         176055             92.72397
## 4      A1B   188854 147768     78.24457         147768             78.24457
```


## Create symbolic links for viewing BAM files in IGV

The `symLink2bam` function creates symbolic links to view the BAM alignment files in a
genome browser such as IGV. The corresponding URLs are written to a file
with a path specified under `urlfile` in the `results` directory.


```r
symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"), 
            urlbase="http://biocluster.ucr.edu/~tgirke/", 
	        urlfile="./results/IGVurl.txt")
```


<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/systemPipeRNAseq_05/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
