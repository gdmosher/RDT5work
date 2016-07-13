~~~ ruby
def what?
aStr = "Hello"
fn(x);
  42
end
~~~

* hello
*   world
* goodbye
    * world


- Four major graphics environments
    - Low-level infrastructure
        - R Base Graphics (low- and high-level)

```r
p <- ggplot(dsmall, aes(carat, price)) + geom_point() + 
            geom_smooth(method="lm", se=FALSE) +
    	    theme(panel.background=element_rect(fill = "white", colour = "black"))
print(p) 
```

``` ruby
def what?
aStr = "Hello"
fn(x);
  42
end
```

``` sublime
def what?
aStr = "Hello"
fn(x);
  42
end
```
- Genome broswer concepts 
    - A genome browser is a visulalization tool for plotting different types of genomic data in separate tracks along chromosomes. 
	- The `ggbio` package [@Yin2012-jj] facilitates plotting of complex genome data objects, such as read alignments (SAM/BAM), genomic context/annotation information (gff/txdb), variant calls (VCF/BCF), and more. To easily compare these data sets, it extends the faceting facility of `ggplot2` to genome browser-like tracks.
	- Most of the core object types for handling genomic data with R/Bioconductor are supported: `GRanges`, `GAlignments`, `VCF`, etc. For more details, see Table 1.1 of the `ggbio` vignette [here](http://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf).
    - `ggbio`'s convenience plotting function is `autoplot`. For more customizable plots, one can use the generic `ggplot` function.
	- Apart from the standard `ggplot2` plotting components, `ggbio` defines serval new components useful for genomic data visualization. A detailed list is given in Table 1.2 of the vignette [here](http://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf). 
	- Useful web sites:
        - [ggbio manual](http://www.tengfei.name/ggbio/docs/)
		- [ggbio functions](http://www.tengfei.name/ggbio/docs/man/)
		- [autoplot demo](http://www.tengfei.name/ggbio/docs/man/autoplot-method.html)


~~~ vim
def what?
aStr = "Hello"
fn(x);
  42
end
~~~
