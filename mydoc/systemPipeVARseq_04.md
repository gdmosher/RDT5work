---
title: Alignments
keywords: 
last_updated: Thu Jul 28 05:13:49 2016
sidebar: mydoc_sidebar
permalink: /mydoc/systemPipeVARseq_04/
---

## Read mapping with `BWA-MEM` 

The NGS reads of this project are aligned against the reference genome
sequence using the highly variant tolerant short read aligner `BWA-MEM`
(Li , 2013; Li et al., 2009). The parameter settings of the aligner are
defined in the `bwa.param` file.


```r
args <- systemArgs(sysma="param/bwa.param", mytargets="targets.txt")
sysargs(args)[1] # Command-line parameters for first FASTQ file
```


Runs the alignments sequentially (_e.g._ on a single machine)


```r
moduleload(modules(args))
system("bwa index -a bwtsw ./data/tair10.fasta")
bampaths <- runCommandline(args=args)
```

Alternatively, the alignment jobs can be submitted to a compute cluster,
here using 72 CPU cores (18 `qsub` processes each with 4 CPU cores).


```r
moduleload(modules(args))
system("bwa index -a bwtsw ./data/tair10.fasta")
resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="10gb")
reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01", 
                  resourceList=resources)
waitForJobs(reg)
writeTargetsout(x=args, file="targets_bam.txt", overwrite=TRUE)
```

Check whether all BAM files have been created


```r
file.exists(outpaths(args))
```

## Read mapping with `gsnap` 

An alternative variant tolerant aligner is `gsnap` from the `gmapR` package
(Wu et al., 2010). The following code shows how to run this aligner on
multiple nodes of a computer cluster that uses Torque as scheduler.


```r
library(gmapR); library(BiocParallel); library(BatchJobs)
args <- systemArgs(sysma="param/gsnap.param", mytargets="targetsPE.txt")
gmapGenome <- GmapGenome(systemPipeR::reference(args), directory="data", name="gmap_tair10chr", create=TRUE)
f <- function(x) {
    library(gmapR); library(systemPipeR)
    args <- systemArgs(sysma="param/gsnap.param", mytargets="targetsPE.txt")
    gmapGenome <- GmapGenome(reference(args), directory="data", name="gmap_tair10chr", create=FALSE)
    p <- GsnapParam(genome=gmapGenome, unique_only=TRUE, molecule="DNA", max_mismatches=3)
    o <- gsnap(input_a=infile1(args)[x], input_b=infile2(args)[x], params=p, output=outfile1(args)[x])
}
funs <- makeClusterFunctionsTorque("torque.tmpl")
param <- BatchJobsParam(length(args), resources=list(walltime="20:00:00", nodes="1:ppn=1", memory="6gb"), cluster.functions=funs)
register(param)
d <- bplapply(seq(along=args), f)
writeTargetsout(x=args, file="targets_gsnap_bam.txt", overwrite=TRUE)
```


## Read and alignment stats

The following generates a summary table of the number of reads in each
sample and how many of them aligned to the reference.


```r
read_statsDF <- alignStats(args=args) 
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
```


## Create symbolic links for viewing BAM files in IGV

The `symLink2bam` function creates symbolic links to view the BAM alignment files in a
genome browser such as IGV. The corresponding URLs are written to a file
with a path specified under `urlfile`, here `IGVurl.txt`.


```r
symLink2bam(sysargs=args, htmldir=c("~/.html/", "projects/gen242/"), 
            urlbase="http://biocluster.ucr.edu/~tgirke/", 
            urlfile="./results/IGVurl.txt")
```


<div class="tags">
<b>Jump to: </b>
<a href="../../mydoc/systemPipeVARseq_05/" class="btn btn-default navbar-btn cursorNorm" role="button">next_page</a>
</div>
