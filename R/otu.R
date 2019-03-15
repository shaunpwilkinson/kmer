#' Cluster sequences into operational taxonomic units.
#'
#' This function performs divisive heirarchical clustering on a set of
#'   DNA sequences using sequential k-means partitioning,
#'   returning an integer vector of OTU membership.
#'
#' @param x a "DNAbin" object.
#' @param k integer giving the k-mer size used to generate the input matrix
#'   for k-means clustering.
#' @param threshold numeric between 0 and 1 giving the OTU identity cutoff.
#'   Defaults to 0.97.
#' @param method the maximum distance criterion to use for terminating the
#'   recursive partitioning procedure. Accepted options are "central" (splitting
#'   stops if the similarity between the central sequence
#'   and its farthest neighbor within the cluster is greater than the threshold),
#'   "centroid" (splitting stops if the similarity between the centroid
#'   and its farthest neighbor within the cluster is greater than the threshold),
#'   and "farthest" (splitting
#'   stops if the similarity between the two farthest sequences within the cluster
#'   is greater than the threshold). Defaults to "central".
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
#' @return a named integer vector of cluster membership with values ranging from 1 to
#'   the total number of OTUs. Asterisks indicate the representative sequence within
#'   each cluster.
#' @details This function clusters sequences into OTUs by first
#'   generating a matrix of k-mer counts, and then splitting the matrix
#'   into two subsets (row-wise) using the k-means algorithm (\emph{k} = 2).
#'   The splitting continues recursively until the farthest k-mer distance
#'   in every cluster is below the threshold value.
#'
#'   This is a divisive, or "top-down" approach to OTU clustering,
#'   as opposed to agglomerative "bottom-up" methods.
#'   It is particularly useful for large large datasets with many sequences
#'   (\emph{n} > 10, 000) since the need to compute a large \emph{n} * \emph{n}
#'   distance matrix is circumvented.
#'   This effectively reduces the time and memory complexity from quadratic to linear,
#'   while generally maintaining comparable accuracy.
#'
#'   It is recommended to increase the value
#'   of \code{nstart} passed to \code{kmeans} \emph{via} the \code{...} argument
#'   to at least 20.
#'   While this can increase computation time, it can improve clustering accuracy
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

#' @author Shaun Wilkinson
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
#' @examples
#' \dontrun{
#' ## Cluster the woodmouse dataset (from the ape package) into OTUs
#' library(ape)
#' data(woodmouse)
#' ## trim gappy ends to subset global alignment
#' woodmouse <- woodmouse[, apply(woodmouse, 2, function(v) !any(v == 0xf0))]
#' ## cluster sequences into OTUs at 0.97 threshold with kmer size = 5
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(999)
#' woodmouse.OTUs <- otu(woodmouse, k = 5, threshold = 0.97, nstart = 20)
#' woodmouse.OTUs
#' }
################################################################################
otu <- function(x, k = 5, threshold = 0.97, method = "central", residues = NULL,
                gap = "-", ...){
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- .alphadetect(x, residues = residues, gap = gap)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  if(is.matrix(x)) x <- .unalign(x, gap = gap)
  if(is.null(names(x))) names(x) <- paste0("SEQ", seq_along(x))
  if(any(duplicated(names(x)))) stop("Sequence names must be unique\n")
  catchnames <- names(x)
  dthresh <- 1 - threshold
  hashes <- .digest(x)
  pointers <- .point(hashes)
  x <- x[!duplicated(hashes)]
  xlengths <- sapply(x, length)
  kcounts <- as.data.frame(kcount(x, k = k, residues = residues, gap = gap, named = FALSE))
  tree <- 1
  attr(tree, "leaf") <- TRUE
  attr(tree, "sequences") <- seq_along(x)
  attr(tree, "height") <- 10
  otun <- function(node, kcounts, seqlengths, k, threshold, ...){
    if(!is.list(node)){
      if(length(attr(node, "sequences")) > 1){
        ## fork leaves only
        seqs <- kcounts[attr(node, "sequences"), , drop = FALSE]
        lens <- seqlengths[attr(node, "sequences")]
        if(method == "farthest"){
          fths <- .farthest2(seqs, k = k, seqlengths = lens)
          if(attr(fths, "distance") <= dthresh){
            attr(node, "central") <- .central1(seqs, k = k, seqlengths = lens, maxdist_central = FALSE)
            return(node)
          }
        }else if(method == "central"){
          cseq <- .central1(seqs, k = k, seqlengths = lens, maxdist_central = TRUE)
          if(attr(cseq, "maxdist_central") <= dthresh){
            attr(node, "central") <- cseq[1]
            return(node)
          }
        }else if(method == "centroid"){
          cseq <- .central1(seqs, k = k, seqlengths = lens, maxdist_central = FALSE)
          if(attr(cseq, "maxdist_centroid") <= dthresh){
            attr(node, "central") <- cseq[1]
            return(node)
          }
        }
        km <- if(nrow(seqs) > 2){
          tryCatch(kmeans(seqs, centers = 2, ... = ...),
                   error = function(er) return(NULL),
                   warning = function(wa) return(NULL))
        }else{
          list(cluster = 1:2)
        }
        if(is.null(km)) return(node)
        tmpattr <- attributes(node)
        node <- vector(mode = "list", length = 2)
        attributes(node) <- tmpattr
        attr(node, "leaf") <- NULL
        for(i in 1:2){
          node[[i]] <- 1
          attr(node[[i]], "height") <- attr(node, "height") - 1
          attr(node[[i]], "leaf") <- TRUE
          attr(node[[i]], "sequences") <- attr(node, "sequences")[km$cluster == i]
        }
      }else{
        attr(node, "central") <- 1
      }
    }
    return(node)
  }
  otur <- function(tree, kcounts, seqlengths, k, threshold, ...){
    tree <- otun(tree, kcounts, seqlengths, k, threshold, ...)
    if(is.list(tree)) tree[] <- lapply(tree, otur, kcounts, seqlengths, k, threshold, ...)
    return(tree)
  }
  ##  build tree recursively
  tree <- otur(tree, kcounts = kcounts, seqlengths = xlengths, k = k,
               threshold = threshold, ... = ...)
  class(tree) <- "dendrogram"
  counter <- 1
  fun <- function(node){
    if(is.leaf(node)){
      repl <- attr(node, "sequences")
      cent <- attr(node, "central") # central sequence
      names(repl) <- rep(counter, length(repl))
      names(repl)[cent] <- paste0(names(repl)[cent], "*")
      node <- repl
      counter <<- counter + 1
    }
    return(node)
  }
  tree <- dendrapply(tree, fun)
  tree <- sort(unlist(tree, use.names = TRUE))
  centrals <- grepl("\\*$", names(tree))
  res <- as.integer(gsub("\\*$", "", names(tree)))
  res <- res[pointers]
  cinds <- match(names(x)[centrals], catchnames) # central indices in derep'd x
  catchnames[cinds] <- paste0(catchnames[cinds], "*")
  names(res) <- catchnames
  return(res)
}
################################################################################
