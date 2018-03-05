# kmer


--------------------------------------------------------------------------------
K-mer counting and clustering for biological sequence analysis  

`kmer` is an R package for rapidly computing distance matrices and 
clustering large sequence datasets using fast alignment-free k-mer counting and 
divisive clustering techniques. 


### Installation
`kmer` is currently available as a development version, with a stable
release available on CRAN shortly. To download the package from 
GitHub users will first need to ensure they have a C/C++ compliler and the 
[devtools](https://github.com/hadley/devtools) R package installed. 
Linux users will generally have a compiler such as `gcc` installed by default; 
however Windows users will need to download 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and Mac 
OSX users will need [Xcode](https://developer.apple.com/xcode) 
(note that Rtools and Xcode are not R packages). To download and install 
devtools, run 
```R
install.packages("devtools")
``` 
and then install and load `kmer` by running 
```R
devtools::install_github("shaunpwilkinson/kmer") 
library("kmer")
```

#### Example: clustering a sequence dataset
The function `cluster` builds a tree by divisive clustering.
This is done by counting k-mers and recursively partitioning 
the sequence set using successive k-means clustering steps. 
No alignment is necessary and no distance matrix is computed,
making it possible to rapidly and efficiently cluster 
very large sequence datasets.

This following code demonstrates how to build and plot a divisive 
tree using the `woodmouse` data from the ape package:

```R
library("kmer")
library("ape")
data(woodmouse)
x <- cluster(woodmouse, k = 5, nstart = 10)
op <- par(no.readonly = TRUE)
par(mar = c(4, 4, 4, 5))
plot(x, horiz = TRUE)
par(op)
```

### Help
An overview of the package and it's functions can be found by running
```R
?kmer
```
If you experience a problem using this package please feel free to
raise it as an issue on [GitHub](http://github.com/shaunpwilkinson/kmer/issues).

### Acknowledgements
This software was developed at 
[Victoria University of Wellington](http://www.victoria.ac.nz/) 
with funding from a Rutherford Foundation Postdoctoral Research Fellowship 
award from the Royal Society of New Zealand.
