## ---- echo = FALSE, message = FALSE, warning = FALSE---------------------
#knitr::opts_chunk$set(out.width='750px', dpi=200)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
library(ape)
data(woodmouse)
## view the first few rows and columns 
as.character.DNAbin(woodmouse[1:5, 1:5])

## ------------------------------------------------------------------------
woodmouse <- woodmouse[, apply(woodmouse, 2, function(v) !any(v == 0xf0))]

## ------------------------------------------------------------------------
### Compute the full distance matrix and print the first few rows and columns
library(kmer)
woodmouse.kdist <- kdistance(woodmouse, k = 5)
print(as.matrix(woodmouse.kdist)[1:7, 1:7], digits = 2)

### Compute and print the embedded distance matrix
set.seed(999)
seeds <- sample(1:15, size = 3)
woodmouse.mbed <- mbed(woodmouse, seeds = seeds, k = 5)
### remove the attributes for printing by subsetting the distance matrix
print(woodmouse.mbed[,], digits = 2)

## ---- echo = FALSE-------------------------------------------------------
op <- par(no.readonly = TRUE)

## ---- out.width='700px', out.height='550px', dpi=500---------------------
## set out plotting panes
par(mfrow = c(1, 3), mar = c(1, 2, 3, 3), cex = 0.3)

## (1) neighbor joining tree with Kimura 1980 distance
### compute the n x n K80 distance matrix 
woodmouse.dist <- ape::dist.dna(woodmouse, model = "K80") 
### build the neighbor-joining tree
tree1.phylo <- ape::nj(woodmouse.dist)
### export as Newick text
tree1.newick <- ape::write.tree(tree1.phylo)
### import as "dendrogram" object
tree1.dendro <- phylogram::read.dendrogram(text = tree1.newick)
### sort nodes by size
tree1.dendro <- phylogram::ladder(tree1.dendro)
### plot the nj tree
plot(tree1.dendro, horiz = TRUE, yaxt = "n", 
     main = "Neighbor-joining tree with\nK80 distance matrix")

## (2) neighbor joining tree with k-mer distance
### compute the n x n k-mer distance matrix 
woodmouse.kdist <- kdistance(woodmouse, k = 5) 
### build the neighbor-joining tree
tree2.phylo <- ape::nj(woodmouse.kdist)
### export as Newick text
tree2.newick <- ape::write.tree(tree2.phylo)
### import as "dendrogram" object
tree2.dendro <- phylogram::read.dendrogram(text = tree2.newick)
### sort nodes by size
tree2.dendro <- phylogram::ladder(tree2.dendro)
### plot the nj tree
plot(tree2.dendro, horiz = TRUE, yaxt = "n", 
     main = "Neighbor-joining tree with\nk-mer distance matrix (k=5)")

## (3) topdown tree without distance matrix
set.seed(999)
tree3 <- cluster(woodmouse, k = 5, nstart = 20)
## sort nodes by size
tree3 <- phylogram::ladder(tree3)
### plot the topdown tree
plot(tree3, horiz = TRUE, yaxt = "n", 
     main = "Top-down tree without\ndistance matrix (k=5)")

## ---- echo = FALSE-------------------------------------------------------
## reset plotting parameters
par(op)

## ------------------------------------------------------------------------
set.seed(999)
woodmouse.OTUs <- otu(woodmouse, k = 5, threshold = 0.97, nstart = 20)
woodmouse.OTUs

