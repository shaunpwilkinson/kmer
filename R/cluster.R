#' Top down (divisive) tree-building.
#'
#' Create phylogenetic trees by successively splitting the sequence dataset
#'   into smaller and smaller subsets.
#'
#' @param x a list or matrix of sequences, possibly an object of class
#'   \code{"DNAbin"} or \code{"AAbin"}.
#' @param k integer. The k-mer size required.
#' @param residues either NULL (default; emitted residues are automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Specifying the residue alphabet is therefore
#'   recommended unless the sequence list is a "DNAbin" or "AAbin" object.
#' @param gap the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @param ... further arguments to be passed to \code{kmeans} (not including
#'   \code{centers}).
#' @return Returns an object of class \code{"dendrogram"}.
#'
#' @details This function creates a tree by successively splitting
#'   the dataset into smaller and smaller subsets (recursive
#'   partitioning). This is a divisive, or "top-down" approach to tree-building,
#'   as opposed to agglomerative "bottom-up" methods such as neighbour joining
#'   and UPGMA. It is particularly useful for large large datasets with many sequences
#'   (\emph{n} > 10, 000) since the need to compute a large \emph{n} * \emph{n}
#'   distance matrix is circumvented.
#'   Instead, a matrix of k-mer counts is computed, and split recursively row-wise
#'   using a k-means clustering algorithm (\emph{k} = 2). This effectively reduces
#'   the time and memory complexity from quadratic to linear, while generally
#'   maintaining comparable accuracy.
#'
#'   If a more accurate tree is required, users can increase the value
#'   of \code{nstart} passed to \code{kmeans} \emph{via} the \code{...} argument.
#'   While this can increase computation time, it can improve tree accuracy
#'   considerably.
#'
#'   DNA and amino acid sequences can be passed to the function either as
#'   a list of non-aligned sequences or a matrix of aligned sequences,
#'   preferably in the "DNAbin" or "AAbin" raw-byte format
#'   (Paradis et al 2004, 2012; see the \code{\link[ape]{ape}} package
#'   documentation for more information on these S3 classes).
#'   Character sequences are supported; however ambiguity codes may
#'   not be recognized or treated appropriately, since raw ambiguity
#'   codes are counted according to their underlying residue frequencies
#'   (e.g. the 5-mer "ACRGT" would contribute 0.5 to the tally for "ACAGT"
#'   and 0.5 to that of "ACGGT").
#'
#'   To minimize computation time when counting longer k-mers (k > 3),
#'   amino acid sequences in the raw "AAbin" format are automatically
#'   compressed using the Dayhoff-6 alphabet as detailed in Edgar (2004).
#'   Note that amino acid sequences will not be compressed if they
#'   are supplied as a list of character vectors rather than an "AAbin"
#'   object, in which case the k-mer length should be reduced
#'   (k < 4) to avoid excessive memory use and computation time.
#'
#' @author Shaun Wilkinson
#'
#' @references
#'   Edgar RC (2004) Local homology recognition and distance measures in
#'   linear time using compressed amino acid alphabets.
#'   \emph{Nucleic Acids Research}, \strong{32}, 380-385.
#'
#'   Paradis E, Claude J, Strimmer K, (2004) APE: analyses of phylogenetics
#'   and evolution in R language. \emph{Bioinformatics} \strong{20}, 289-290.
#'
#'   Paradis E (2012) Analysis of Phylogenetics and Evolution with R
#'   (Second Edition). Springer, New York.
#'
#' @seealso \code{\link{kcount}}
#'
#' @examples
#' \dontrun{
#' ## Cluster the woodmouse dataset (ape package)
#' library(ape)
#' data(woodmouse)
#' ## trim gappy ends to subset global alignment
#' woodmouse <- woodmouse[, apply(woodmouse, 2, function(v) !any(v == 0xf0))]
#' ## build tree divisively
#' set.seed(999)
#' woodmouse.tree <- cluster(woodmouse, nstart = 5)
#' ## plot tree
#' op <- par(no.readonly = TRUE)
#' par(mar = c(5, 2, 4, 8) + 0.1)
#' plot(woodmouse.tree, main = "Woodmouse phylogeny", horiz = TRUE)
#' par(op)
#' }
################################################################################
cluster <- function(x, k = 5, residues = NULL, gap = "-", ...){
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- .alphadetect(x, residues = residues, gap = gap)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  if(is.matrix(x)) x <- .unalign(x, gap = gap)
  nseq <- length(x)
  if(nseq == 1){# singleton tree (leaf)
    tree <- 1
    attr(tree, "leaf") <- TRUE
    attr(tree, "height") <- 0
    attr(tree, "midpoint") <- 0
    attr(tree, "label") <- names(x)[1]
    attr(tree, "members") <- 1
    class(tree) <- "dendrogram"
    return(tree)
  }
  # stop("Only a single sequence provided")
  if(is.null(names(x))) names(x) <- paste0("S", 1:nseq)
  catchnames <- names(x)
  hashes <- .digest(x, simplify = TRUE)
  duplicates <- duplicated(hashes)
  if(sum(!duplicates) == 1){
    tree <- vector(mode = "list", length = length(x))
    attr(tree, "height") <- 0
    attr(tree, "midpoint") <- 0.5
    attr(tree, "members") <- length(x)
    for(i in seq_along(x)){
      tree[[i]] <- i
      attr(tree[[i]], "height") <- 0
      attr(tree[[i]], "midpoint") <- 0
      attr(tree[[i]], "members") <- 1
      attr(tree[[i]], "leaf") <- TRUE
      attr(tree[[i]], "label") <- catchnames[i]
    }
    class(tree) <- "dendrogram"
    return(tree)
  }
  nuseq <- sum(!duplicates)
  if(any(duplicates)){
    pointers <- integer(length(x))
    dupehashes <- hashes[duplicates]
    uniquehashes <- hashes[!duplicates]
    pointers[!duplicates] <- seq_along(uniquehashes)
    pd <- integer(length(dupehashes))
    for(i in unique(dupehashes)) pd[dupehashes == i] <- match(i, uniquehashes)
    pointers[duplicates] <- pd
    x <- x[!duplicates]
  }else{
    pointers <- seq_along(x)
  }
  kcounts <- kcount(x, k = k, residues = residues, gap = gap, named = FALSE)
  kfreqs <- kcounts/apply(kcounts, 1, sum)
  tree <- 1
  attr(tree, "leaf") <- TRUE
  attr(tree, "sequences") <- 1:nuseq
  attr(tree, "height") <- 0
  ## define recursive splitting functions
  cluster1 <- function(tree, d, ...){ # d is the kfreq matrix
    tree <- cluster2(tree, d = d, ... = ...)
    if(is.list(tree)) tree[] <- lapply(tree, cluster1, d = d, ... = ...)
    return(tree)
  }
  cluster2 <- function(node, d, ...){
    if(!is.list(node) & length(attr(node, "sequences")) > 1){
      ## fork leaves only
      seqs <- d[attr(node, "sequences"), , drop = FALSE]
      errfun <- function(er){## used when >3 uniq hashes but kmeans throws error
        out <- list()
        nrs <- nrow(seqs)
        cls <- rep(1, nrs)
        cls[sample(1:nrs, 1)] <- 2 ## peel randomly selected one off
        out$cluster <- cls
        out$centers <- rbind(apply(seqs[cls == 1, , drop = FALSE], 2, mean),
                             apply(seqs[cls == 2, , drop = FALSE], 2, mean))
        return(out)
      }
      km <- if(nrow(seqs) > 2) {
        tryCatch(kmeans(seqs, centers = 2, ... = ...), error = errfun, warning = errfun)
      }else{
        list(cluster = 1:2, centers = seqs)
      }
      membership <- km$cluster
      centers <- km$centers
      tmpattr <- attributes(node)
      node <- vector(mode = "list", length = 2)
      attributes(node) <- tmpattr
      attr(node, "leaf") <- NULL
      attr(node, "avdist") <- sqrt(sum(abs(km$centers[1,] - km$centers[2,])^2))
      for(i in 1:2){
        node[[i]] <- 1
        attr(node[[i]], "height") <- attr(node, "height") - 0.0001 ## cleaned up later
        attr(node[[i]], "leaf") <- TRUE
        attr(node[[i]], "sequences") <- attr(node, "sequences")[membership == i]
        attr(node[[i]], "avdist") <- 0 ## recalculated later if subnode is list
      }
    }
    return(node)
  }
  ##  build tree recursively
  tree <- cluster1(tree, d = kfreqs)
  tree <- phylogram::remidpoint(tree)
  class(tree) <- "dendrogram"
  reheight <- function(node){
    if(is.list(node)){
      ## arbitrary minimum branch length of 0.0001
      node1edge <- max(0.0001, (attr(node, "avdist") - attr(node[[1]], "avdist"))/2)
      node2edge <- max(0.0001, (attr(node, "avdist") - attr(node[[2]], "avdist"))/2)
      attr(node[[1]], "height") <- attr(node, "height") - node1edge
      attr(node[[2]], "height") <- attr(node, "height") - node2edge
    }
    attr(node, "avdist") <- NULL ## tidy up by removing attrs
    return(node)
  }
  reheight1 <- function(node){
    node <- reheight(node)
    if(is.list(node)) node[] <- lapply(node, reheight1)
    return(node)
  }
  tree <- reheight1(tree)
  class(tree) <- "dendrogram"
  tree <- phylogram::reposition(tree)
  if(any(duplicates)){
    reduplicate <- function(node, pointers){
      attr(node, "sequences") <- which(pointers %in% attr(node, "sequences"))
      if(is.leaf(node)){
        lams <- length(attr(node, "sequences"))
        if(lams > 1){
          labs <- attr(node, "label")
          hght <- attr(node, "height")
          seqs <- attr(node, "sequences")
          node <- vector(mode = "list", length = lams)
          attr(node, "height") <- hght
          attr(node, "sequences") <- seqs
          for(i in 1:lams){
            node[[i]] <- 1
            attr(node[[i]], "height") <- hght
            attr(node[[i]], "label") <- labs[i]
            attr(node[[i]], "sequences") <- seqs[i]
            attr(node[[i]], "leaf") <- TRUE
          }
        }
      }
      return(node)
    }
    tree <- dendrapply(tree, reduplicate, pointers)
    tree <- remidpoint(tree)
  }
  label <- function(node, labs){
    if(is.leaf(node)) attr(node, "label") <- labs[attr(node, "sequences")]
    return(node)
  }
  tree <- dendrapply(tree, label, labs = catchnames)
  rmseqs <- function(node){
    if(is.leaf(node)){
      tmpattr <- attributes(node)
      node[] <- tmpattr$sequences
      tmpattr$sequences <- NULL
      attributes(node) <- tmpattr
    }else{
      attr(node, "sequences") <- NULL
    }
    return(node)
  }
  tree <- dendrapply(tree, rmseqs)
  return(tree)
}
################################################################################
