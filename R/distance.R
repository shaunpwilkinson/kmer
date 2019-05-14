#' K-mer counting.
#'
#' Count all k-letter words in a sequence or set of sequences
#'   with a sliding window of length k.
#'
#' @param x a matrix of aligned sequences, a list of unaligned sequences,
#'   or a vector representing a single sequence.
#'   Accepted modes are "character" and "raw" (the latter being applicable
#'   for "DNAbin" and "AAbin" objects).
#' @param k integer representing the k-mer size. Defaults to 5.
#'   Note that high values of k may be slow to compute and use a lot of
#'   memory due to the large numbers of calculations required,
#'   particularly when the residue alphabet is also large.
#' @param residues either NULL (default; the residue alphabet is automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Specifying the residue alphabet is therefore
#'   recommended unless x is a "DNAbin" or "AAbin" object.
#' @param gap the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} and \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @param named logical. Should the k-mers be returned as column names in
#'   the returned matrix? Defaults to TRUE.
#' @param compress logical indicating whether to compress AAbin sequences
#'   using the Dayhoff(6) alphabet for k-mer sizes exceeding 4.
#'   Defaults to TRUE to avoid memory overflow and excessive computation time.
#' @param encode logical indicating if the resulting matrix should be encoded
#'   in raw bytes (output matrix can be decoded with \code{kmer:::.decodekc()}).
#'   Note that the output will be rounded and have maximum k-mer count of 15.
#' @return Returns a matrix of k-mer counts with one row for each sequence
#'   and \emph{n}^\emph{k} columns (where \emph{n} is the size of the
#'   residue alphabet and \emph{k} is the k-mer size)
#' @details
#'   This function computes a vector or matrix of k-mer counts
#'   from a sequence or set of sequences using a sliding a window of length k.
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
#' @seealso \code{\link{kdistance}} for k-mer distance matrix computation.
#'
#' @examples
#'   ## compute a matrix of k-mer counts for the woodmouse
#'   ## data (ape package) using a k-mer size of 3
#'   library(ape)
#'   data(woodmouse)
#'   x <- kcount(woodmouse, k = 3)
#'   x
#'   ## 64 columns for nucleotide 3-mers AAA, AAC, ... TTT
#'   ## convert to AAbin object and repeat the operation
#'   y <- kcount(ape::trans(woodmouse, 2), k = 2)
#'   y
#'   ## 400 columns for amino acid 2-mers AA, AB, ... , YY
################################################################################
kcount <- function(x, k = 5, residues = NULL, gap = "-", named = TRUE,
                   compress = TRUE, encode = FALSE){

  DNA <- .isDNA(x)
  AA <- .isAA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- .alphadetect(x, residues = residues, gap = gap)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  if(is.matrix(x)) x <- .unalign(x, gap = gap)
  if(!is.list(x)) x <- list(x)
  x <- lapply(x, function(s) s[s != gap])
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  ## divide job to prevent overflow
  if(length(x) > 10000){
    nmats <- length(x) %/% 10000 + 1
    f <- rep(seq_len(nmats), each = 10000)[seq_along(x)]
    y <- split(x, f)
    kcounts <- lapply(y, kcount, k, residues, gap, named, compress, encode)
    kcounts <- do.call("rbind", kcounts)
    return(kcounts)
  }
  nseq <- length(x)
  if(is.null(names(x))) names(x) <- paste0("S", 1:nseq)
  seqalongx <- seq_along(x)
  if(DNA){
    if(k > 12) stop("Maximum kmer size is 12 for DNA\n")
    arity <- 4
    x <- lapply(x, function(s) s[!(s %in% as.raw(c(2, 4, 240)))])
    seqlengths <- sapply(x, length)
    if(min(seqlengths) < k) stop("minimum sequence length is less than k")
    kcounts <- .kcountDNA(x, k = k)
  }else{
    tuplecount <- function(y, k, arity){
      tuplemat <- matrix(nrow = k, ncol = length(y) - k + 1)
      for(i in 1:k) tuplemat[i, ] <- y[i:(length(y) - (k - i))]
      res <- apply(tuplemat, 2, .decimal, from = arity) + 1
      res <- tabulate(res, nbins = arity^k)
      return(res)
    }
    if(AA){
      if(k > 3 & compress){
        message("Converting to Dayhoff(6) compressed alphabet for k > 3")
        message("Classes: AGPST, C, DENQ, FWY, HKR, ILMV\n")
        residues <- c("A", "C", "D", "F", "H", "I")
        arity <- 6
      }else{
        arity <- 20
      }
      x <- .encodeAA(x, arity = arity, na.rm = TRUE)
    }else{
      arity <- length(residues)
      if(k > 2 & arity >= 20) stop("Unable to calculate distance matrix for
                               large k and large alphabet size. If residues
                               are amino acids consider converting to AAbin
                               object for compression")
      modes <- lapply(x, mode)
      if(!(all(modes == "integer") | all(modes == "numeric"))){
        x <- .encodeCH(x, residues = residues, na.rm = TRUE)
      }
    }
    seqlengths <- sapply(x, length)
    if(min(seqlengths) < k) stop("minimum sequence length is less than k")
    kcounts <- t(sapply(x, tuplecount, k, arity))
  }
  # label columns with k-mer words
  indices <- matrix(1L, nrow = k, ncol = arity^k)
  ntimes <- 1
  counter <- k - 1
  for(i in 1:k){
    indices[i, ] <- rep(rep(1:arity, each = arity^counter), times = ntimes)
    ntimes <- ntimes * arity
    counter <- counter - 1
  }
  colnames(kcounts) <- if(named){
    apply(indices, 2, function(i) paste0(residues[i], collapse = ""))
  }else NULL
  if(encode) kcounts <- .encodekc(kcounts)
  return(kcounts)
}
################################################################################
#' K-mer distance matrix computation.
#'
#' Computes the matrix of k-mer distances between all pairwise comparisons
#'   of a set of sequences.
#'
#' @param x a matrix of aligned sequences or a list of unaligned sequences.
#'   Accepted modes are "character" and "raw" (the latter being applicable
#'   for "DNAbin" and "AAbin" objects).
#' @param k integer representing the k-mer size to be used for calculating
#'   the distance matrix. Defaults to 5. Note that high values of k
#'   may be slow to compute and use a lot of memory due to the large numbers
#'   of calculations required, particularly when the residue alphabet is
#'   also large.
#' @param method a character string giving the k-mer distance measure
#'   to be used. Currently the available options are \code{"edgar"} (default;
#'   see Edgar (2004) for details) and the standard methods available for
#'   the base function "dist" ("euclidean", "maximum", "manhattan", "canberra",
#'   "binary" and "minkowski").
#' @param residues either NULL (default; the residue alphabet is automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Specifying the residue alphabet is therefore
#'   recommended unless x is a "DNAbin" or "AAbin" object.
#' @param gap the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @param compress logical indicating whether to compress AAbin sequences
#'   using the Dayhoff(6) alphabet for k-mer sizes exceeding 4.
#'   Defaults to TRUE to avoid memory overflow and excessive computation time.
#' @param ... further arguments to be passed to \code{"as.dist"}.
#' @return an object of class \code{"dist"}.
#' @details
#'   This function computes the \emph{n} * \emph{n} k-mer distance matrix
#'   (where \emph{n} is the number of sequences), returning an object of class
#'   \code{"dist"}. DNA and amino acid sequences can be passed to the function
#'   either as a list of non-aligned sequences or as a matrix of aligned sequences,
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
#' @seealso \code{\link{kcount}} for k-mer counting, and
#'   \code{\link{mbed}} for leaner distance matrices
#'
#' @examples
#'   ## compute a k-mer distance matrix for the woodmouse
#'   ## dataset (ape package) using a k-mer size of 5
#'   library(ape)
#'   data(woodmouse)
#'   ### subset global alignment by removing gappy ends
#'   woodmouse <- woodmouse[, apply(woodmouse, 2, function(v) !any(v == 0xf0))]
#'   ### compute the distance matrix
#'   woodmouse.dist <- kdistance(woodmouse, k = 5)
#'   ### cluster and plot UPGMA tree
#'   woodmouse.tree <- as.dendrogram(hclust(woodmouse.dist, "average"))
#'   plot(woodmouse.tree)
################################################################################
kdistance <- function(x, k = 5, method = "edgar", residues = NULL,
                      gap = "-", compress = TRUE, ...){
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- .alphadetect(x, residues = residues, gap = gap)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  if(is.matrix(x)) x <- .unalign(x, gap = gap)
  nseq <- length(x)
  if(is.null(names(x))) names(x) <- paste0("S", 1:nseq)
  seqalongx <- seq_along(x)
  kcounts <- kcount(x, k = k, residues = residues, gap = gap, compress = compress)
  seqlengths <- apply(kcounts, 1, sum) + k - 1
  ## not sapply(x, length) in case of gaps, unknowns etc
  ## which will be picked up by kcount
  if(method == "edgar"){
    d <- .kdist(kcounts, from = seqalongx - 1, to = seqalongx - 1,
                        seqlengths = seqlengths, k = k)
  }else{
    freqs <- kcounts/(seqlengths - k + 1)
    d <- dist(freqs, method = method)
  }
  return(as.dist(d, ... = ...))
}
################################################################################
#' Convert sequences to vectors of distances to a subset of seed sequences.
#'
#' This function computes a matrix of
#'   distances from each sequence to a subset of 'seed' sequences using
#'   the method outlined in Blacksheilds et al (2010).
#'
#' @param x a matrix of aligned sequences or a list of unaligned sequences.
#'   Accepted modes are "character" and "raw" (the latter is for "DNAbin"
#'   and "AAbin" objects).
#' @param k integer representing the k-mer size to be used for calculating
#'   the distance matrix. Defaults to 5. Note that high values of k
#'   may be slow to compute and use a lot of memory due to the large numbers
#'   of calculations required, particularly when the residue alphabet is
#'   also large.
#' @param seeds optional integer vector indicating which sequences should
#'   be used as the seed sequences. If \code{seeds = NULL} a set of
#'   log(\emph{n}, 2)^2 non-identical sequences is randomly selected from the
#'   sequence set (where \emph{n} is the number of sequences; see Blacksheilds et al.
#'   2010). Alternatively, if \code{seeds = 'all'} a standard \emph{n} * \emph{n}
#'   distance matrix is computed.
#' @param residues either NULL (default; emitted residues are automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Specifying the residue alphabet is therefore
#'   recommended unless x is a "DNAbin" or "AAbin" object.
#' @param gap the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @param counts logical indicating whether the (usually large) matrix of
#'   k-mer counts should be returned as an attribute of the returned
#'   object. Defaults to FALSE.
#' @return Returns an object of class \code{"mbed"}, whose primary object is
#'   an \emph{n} * log(\emph{n}, 2)^2 matrix
#'   (where \emph{n} is the number of sequences). The returned
#'   object contains additional attributes including an
#'   integer vector of seed sequence indices ("seeds"), a logical vector
#'   identifying the duplicated sequences ("duplicates"), an integer vector
#'   giving the matching indices of the non-duplicated sequences ("pointers"),
#'   a character vector of MD5 digests of the sequences ("hashes"),
#'   an integer vector of sequence lengths ("seqlengths"), and if
#'   \code{counts = TRUE}, the matrix of k-mer counts ("kcounts";
#'   see \code{\link{kcount}} for details).
#' @details
#'   This function computes a \emph{n} * log(\emph{n}, 2)^2 k-mer distance matrix
#'   (where \emph{n} is the number of sequences), returning an object of class
#'   \code{"mbed"}. If the number of sequences is less than or equal to 19, the full
#'   \emph{n} * \emph{n} distance matrix is produced (since the rounded up value of
#'   log(\emph{19}, 2)^2 is 19). Currently the only distance measure supported is
#'   that of Edgar (2004).
#'
#'   For maximum information retention following the embedding process
#'   it is generally desirable to select the seed sequences based on their
#'   uniqueness, rather than simply selecting a random subset
#'   (Blackshields et al. 2010).
#'   Hence if 'seeds' is set to NULL (the default setting) the the `mbed`
#'   function selects the subset by clustering the sequence set into
#'   \emph{t} groups using the k-means algorithm (\emph{k} = \emph{t}),
#'   and choosing one representative from each group.
#'   Users can alternatively pass an integer vector (as in the above example)
#'   to specify the seeds manually. See Blackshields et al (2010) for other
#'   seed selection options.
#'
#'   DNA and amino acid sequences can be passed to the function
#'   either as a list of non-aligned sequences or as a matrix of aligned sequences,
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
#'   Note that agglomerative (bottom-up) tree-building methods
#'   such as neighbor-joining and UPGMA depend on a full
#'   \emph{n} * \emph{n} distance matrix.
#'   See the \code{\link{kdistance}} function for details on computing
#'   symmetrical distance matrices.
#'
#' @author Shaun Wilkinson
#'
#' @references
#'   Blackshields G, Sievers F, Shi W, Wilm A, Higgins DG (2010) Sequence embedding
#'   for fast construction of guide trees for multiple sequence alignment.
#'   \emph{Algorithms for Molecular Biology}, \strong{5}, 21.
#'
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
#' @seealso \code{\link{kdistance}} for full \emph{n} * \emph{n} distance
#'   matrix computation.
#'
#' @examples
#'   ## compute an embedded k-mer distance matrix for the woodmouse
#'   ## dataset (ape package) using a k-mer size of 5
#'   library(ape)
#'   data(woodmouse)
#'   ## randomly select three sequences as seeds
#'   suppressWarnings(RNGversion("3.5.0"))
#'   set.seed(999)
#'   seeds <- sample(1:15, size = 3)
#'   ## embed the woodmouse dataset in three dimensions
#'   woodmouse.mbed <- mbed(woodmouse, seeds = seeds, k = 5)
#'   ## print the distance matrix (without attributes)
#'   print(woodmouse.mbed[,], digits = 2)
################################################################################
mbed <- function(x, seeds = NULL, k = 5, residues = NULL, gap = "-",
                 counts = FALSE){
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- .alphadetect(x, residues = residues, gap = gap)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  if(is.matrix(x)) x <- .unalign(x, gap = gap)
  nseq <- length(x)
  if(is.null(names(x))) names(x) <- paste0("S", 1:nseq)
  if(!is.null(seeds)){
    if(identical(seeds, "all")) seeds <- seq_along(x)
    stopifnot(mode(seeds) %in% c("numeric", "integer"),
              max(seeds) <= nseq,
              min(seeds) > 0)
  }
  hashes <- .digest(x)
  duplicates <- duplicated(hashes)
  nuseq <- sum(!duplicates)
  pointers <- .point(hashes)
  catchnames <- names(x)
  x <- x[!duplicates]
  seqalongx <- seq_along(x)
  if(DNA){
    x <- lapply(x, function(s) s[!(s %in% as.raw(c(2, 4, 240)))])
    seqlengths <- sapply(x, length)
    if(min(seqlengths) < k) stop("minimum sequence length is less than k")
    kcounts <- .kcountDNA(x, k = k)
  }else{
    tuplecount <- function(y, k, arity){
      tuplemat <- matrix(nrow = k, ncol = length(y) - k + 1)
      for(i in 1:k) tuplemat[i, ] <- y[i:(length(y) - (k - i))]
      res <- apply(tuplemat, 2, .decimal, from = arity) + 1
      res <- tabulate(res, nbins = arity^k)
      return(res)
    }
    if(AA){
      arity <- if(k > 2) 6 else 20 # compress AA alphabet for high k values
      x <- .encodeAA(x, arity = arity, na.rm = TRUE)
    }else{
      arity <- length(residues)
      if(k > 2 & arity >= 20) stop("Unable to calculate distance matrix for
                               large k and large alphabet size. If residues
                               are amino acids consider converting to AAbin
                               object for automatic compression, else reduce k")
      modes <- lapply(x, mode)
      if(!(all(modes == "integer") | all(modes == "numeric"))){
        x <- .encodeCH(x, residues = residues, na.rm = TRUE)
      }
    }
    seqlengths <- sapply(x, length)
    if(min(seqlengths) < k) stop("minimum sequence length is less than k")
    kcounts <- t(sapply(x, tuplecount, k, arity))
  }
  if(is.null(seeds)){
    nseeds <- ceiling(log(nuseq, 2)^2)
    allseeds <- nseeds >= nuseq
    if(allseeds){
      nseeds <- nuseq
      seeds <- seqalongx
    }else{
      suppressWarnings(groups <- kmeans(kcounts, centers = nseeds)$cluster)
      seeds <- match(1:nseeds, groups)
    }
    ## LLR algorithm see Blacksheilds et al. 2010
  }else{
    seeds <- unique(pointers[seeds])
    nseeds <- length(seeds)
  }
  res <- .kdist(kcounts, from = seqalongx - 1, to = seeds - 1,
                seqlengths = seqlengths, k = k)
  if(any(duplicates)){
    tmp <- matrix(nrow = nseq, ncol = ncol(res))
    rownames(tmp) <- catchnames
    colnames(tmp) <- names(x)[seeds]
    tmp[!duplicates, ] <- res
    tmp[duplicates, ] <- res[pointers[duplicates], ]
    res <- tmp
    if(counts){
      tmpkc <- matrix(nrow = nseq, ncol = ncol(kcounts))
      rownames(tmpkc) <- catchnames
      tmpkc[!duplicates, ] <- kcounts
      tmpkc[duplicates, ] <- kcounts[pointers[duplicates], ]
      kcounts <- tmpkc
    }
  }
  attr(res, "seeds") <- seeds # integer vector
  if(counts) attr(res, "kcounts") <- kcounts
  rm(kcounts)
  gc()
  attr(res, "duplicates") <- duplicates
  attr(res, "pointers") <- pointers
  attr(res, "hashes") <- hashes
  attr(res, "seqlengths") <- seqlengths
  class(res) <- "mbed"
  return(res)
}
################################################################################
