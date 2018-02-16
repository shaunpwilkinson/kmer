#' Cluster sequences into operational taxonomic units.
#'
#' This function performs divisive heirarchical clustering on a set of
#'   DNA sequences using sequential k-means partitioning,
#'   returning an integer vector of OTU membership.
#'
#' @param x a "DNAbin" object.
#' @param k integer giving the k-mer size used to generate the input matrix
#'   for k-means clustering.
#' @param threshold numeric between 0 and 1 giving the OTU similarity cutoff.
#' @param ... further arguments to be passed to \code{kmeans}.
#' @return an integer vector of cluster membership with values from 1 to
#'   the total number of OTUs.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples
#'   ##TBA
#'
################################################################################
otu <- function(x, k = 5, threshold = 0.97, ...){
  nseq <- length(x)
  hashes <- .digest(x, simplify = TRUE)
  duplicates <- duplicated(hashes)
  nuseq <- sum(!duplicates)
  if(any(duplicates)){
    pointers <- .point(hashes)
    x <- x[!duplicates]
  }else{
    pointers <- seq_along(x)
  }
  kcounts <- round(kcount(x, k = k, named = FALSE))
  tree <- 1
  attr(tree, "leaf") <- TRUE
  attr(tree, "sequences") <- 1:nuseq
  attr(tree, "height") <- 10
  otu1 <- function(tree, d, threshold, ...){ # d is the kcount matrix
    tree <- otu2(tree, d = d, threshold = threshold, ... = ...)
    if(is.list(tree)) tree[] <- lapply(tree, otu1, d = d, threshold = threshold, ... = ...)
    return(tree)
  }
  otu2 <- function(node, d, threshold, ...){
    if(!is.list(node) & length(attr(node, "sequences")) > 1){
      ## fork leaves only
      seqs <- d[attr(node, "sequences"), , drop = FALSE]
      conserved <- apply(seqs, 2, function(v) all(v == v[1]))
      propcons <- sum(conserved)/length(conserved)
      if(propcons > threshold) return(node)
      km <- if(nrow(seqs) > 2){
        tryCatch(kmeans(seqs, centers = 2, ... = ...),
                 error = function(er) return(NULL),
                 warning = function(wa) return(NULL))
      }else{
        list(cluster = 1:2)
      }
      if(is.null(km)) return(node)
      membership <- km$cluster
      tmpattr <- attributes(node)
      node <- vector(mode = "list", length = 2)
      attributes(node) <- tmpattr
      attr(node, "leaf") <- NULL
      for(i in 1:2){
        node[[i]] <- 1
        attr(node[[i]], "height") <- attr(node, "height") - 1
        attr(node[[i]], "leaf") <- TRUE
        attr(node[[i]], "sequences") <- attr(node, "sequences")[membership == i]
      }
    }
    return(node)
  }
  ##  build tree recursively
  tree <- otu1(tree, d = kcounts, threshold = threshold)
  #tree <- phylogram::remidpoint(tree)
  class(tree) <- "dendrogram"
  counter <- 1
  fun <- function(node){
    if(is.leaf(node)){
      repl <- attr(node, "sequences")
      names(repl) <- rep(counter, length(repl))
      node <- repl
      counter <<- counter + 1
    }
    return(node)
  }
  tree <- dendrapply(tree, fun)
  tree <- unlist(tree, use.names = TRUE)
  res <- as.integer(names(sort(tree)))
  res <- res[pointers]
  return(res)
}
################################################################################
