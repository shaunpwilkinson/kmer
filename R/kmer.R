#' Fast K-mer Counting and Clustering for Biological Sequence Analysis.
#'
#' The kmer package contains tools for rapidly computing
#'   distance matrices, building large trees, and clustering
#'   operational taxonomic units using fast alignment-free
#'   k-mer counting and divisive clustering techniques.
#'
#' @section Functions:
#' A breif description of the primary \pkg{kmer} functions are
#'   provided with links to their help pages below.
#'
#'
#' @section K-mer counting:
#' \itemize{
#' \item \code{\link{kcount}} counts all k-letter words in a
#'   sequence or set of sequences using a sliding window of
#'   length k
#' }
#'
#' @section Distance matrix computation:
#' \itemize{
#' \item \code{\link{kdistance}} calculates pairwise
#'   distances between sequences based on k-mer counts
#' \item \code{\link{mbed}} embeds sequences as vectors of
#'   k-mer distances to a set of seed' sequences
#' }
#'
#' @section Alignment-free clustering:
#' \itemize{
#' \item \code{\link{cluster}} builds a phylogenetic tree by
#'   successively splitting a set of sequences
#'   (recursive partitioning) based on k-mer counts
#' \item \code{\link{otu}} heirarchically clusters a set of sequences
#'   until a predefined furthest neighbor dissimilarity threshold is reached.
#' }
#'
#'
#' @docType package
#' @name kmer
################################################################################
NULL
