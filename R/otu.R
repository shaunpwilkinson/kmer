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
#' @return an integer vector of cluster membership with values from 1 to
#'   the total number of OTUs.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples
#'   ##TBA
#'
################################################################################
otu <- function(x, k = 5, threshold = 0.97){
  dthresh <- 1 - threshold
  hashes <- .digest(x)
  pointers <- .point(hashes)
  x <- x[!duplicated(hashes)]
  xlengths <- sapply(x, length)
  kcounts <- as.data.frame(kcount(x, k = k))
  tree <- 1
  attr(tree, "leaf") <- TRUE
  attr(tree, "sequences") <- seq_along(x)
  attr(tree, "height") <- 10
  farthest2 <- function(y, k, seqlengths){ # y is a kcount (or other) matrix
    y <- as.matrix(y) # in case y is a df
    point1 <- sample(1:nrow(y), size = 1)
    checked <- integer(100)
    checked[1] <- point1
    dists <- .kdist(y, from = 1:nrow(y) - 1, to = point1 - 1,
                    seqlengths = seqlengths, k = k)[, 1]
    point2 <- which.max(dists)
    for(i in 2:100){
      dists <- .kdist(y, from = 1:nrow(y) - 1, to = point2 - 1,
                      seqlengths = seqlengths, k = k)[, 1]
      point3 <- which.max(dists)
      if(point3 %in% checked) break
      point1 <- point2
      checked[i] <- point1
      point2 <- point3
    }
    if(i == 100) stop("Farthest distances not found\n")
    res <- c(point2, point3)
    attr(res, "distance") <- unname(dists[point3])
    return(res)
  }
  otun <- function(node, kcs, seqlengths, k, threshold){
    if(!is.list(node) & length(attr(node, "sequences")) > 1){
      ## fork leaves only
      seqs <- kcs[attr(node, "sequences"), , drop = FALSE]
      lens <- seqlengths[attr(node, "sequences")]
      fths <- farthest2(seqs, k = k, seqlengths = lens)
      if(attr(fths, "distance") <= dthresh) return(node)
      km <- if(nrow(seqs) > 2){
        tryCatch(kmeans(seqs, centers = 2, nstart = 20),
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
    }
    return(node)
  }
  otur <- function(tree, kcs, seqlengths, k, threshold){ # kcs is the kcount matrix
    tree <- otun(tree, kcs, seqlengths, k, threshold)
    if(is.list(tree)) tree[] <- lapply(tree, otur, kcs, seqlengths, k, threshold)
    return(tree)
  }
  ##  build tree recursively
  tree <- otur(tree, kcs = kcounts, seqlengths = xlengths, k = k, threshold = threshold)
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
