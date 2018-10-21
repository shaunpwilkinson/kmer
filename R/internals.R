# Internal functions.

# Detect residue alphabet
.alphadetect <- function(sequences, residues = NULL, gap = "-",
                         endchar = "?"){
  if(identical(toupper(residues), "RNA")){
    residues <- c("A", "C", "G", "U")
  }else if(.isDNA(sequences) | identical(toupper(residues), "DNA")){
    residues <- c("A", "C", "G", "T")
  }else if(.isAA(sequences) | identical(residues, "AA") |
           identical (toupper(residues), "AMINO")){
    residues <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
  }
  else if(is.null(residues)){
    residues <- sort(unique(as.vector(unlist(sequences, use.names = FALSE))))
    if(!is.null(gap)) residues <- residues[residues != gap]
    if(!is.null(endchar)) residues <- residues[residues != endchar]
  }else{
    if(!is.null(gap)) residues <- residues[residues != gap]
    if(!is.null(endchar)) residues <- residues[residues != endchar]
  }
  if(!(length(residues) > 0)){# & mode(residues) == "character")){
    stop("invalid residues argument")
  }
  return(residues)
}

# Convert a vector in any arity to a decimal integer
.decimal <- function(x, from) sum(x * from^rev(seq_along(x) - 1))

# Check if object is DNA
.isDNA <- function(x){
  if(inherits(x, "DNAbin")){
    return(TRUE)
  }else if(inherits(x, "AAbin")){
    return(FALSE)
  }else if(mode(x) == "character"){
    return(FALSE)
  }else if(mode(x) == "raw"){
    return(all(x %in% as.raw(c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48,
                               224, 176, 208, 112, 240, 4, 2))))
  }else if(mode(x) == "list"){
    if(length(x) > 0){
      return(all(unlist(x, use.names = FALSE) %in%
                   as.raw(c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48,
                                         224, 176, 208, 112, 240, 4, 2))))
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

# Check if object is amino acid sequence
.isAA <- function(x){
  if(inherits(x, "AAbin")){
    return(TRUE)
  }else if(inherits(x, "DNAbin")){
    return(FALSE)
  }else if(mode(x) == "character"){
    return(FALSE)
  }else if(mode(x) == "raw"){
    return(all(x %in% as.raw(c(65:90, 42, 45, 63))))
  }else if(mode(x) == "list"){
    if(length(x) > 0){
      return(all(unlist(x, use.names = FALSE) %in% as.raw(c(65:90, 42, 45, 63))))
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

# Remove ambiguous residues from amino acid sequence

.disambiguateAA <- function(a, probs = rep(0.05, 20), random = TRUE){
  # a is a raw byte in AAbin format
  guide <- as.raw(c(65:90, 42, 45)) #length = 28
  nonambigs <- guide[1:25][-c(2, 10, 15, 21, 24)]
  #structure(guide, class = "AAbin")
  if(a == guide[24]){
    if(random){
      sample(nonambigs, size = 1, prob = probs)
    }else{
      return(nonambigs[which.max(probs)])
    }
  }else if(a == guide[2]){# B
    if(random){
      sample(nonambigs[c(3, 12)], size = 1, prob = probs[c(3, 12)]) # D or N
    }else{
      return(nonambigs[which.max(probs[c(3, 12)])])
    }
  }else if(a == guide[10]){# J
    if(random){
      sample(nonambigs[c(8, 10)], size = 1, prob = probs[c(8, 10)]) # I or L
    }else{
      return(nonambigs[which.max(probs[c(8, 10)])])
    }
  }else if(a == guide[26]){ # Z
    if(random){
      sample(nonambigs[c(4, 14)], size = 1, prob = probs[c(4, 14)]) # E or Q
    }else{
      return(nonambigs[which.max(probs[c(4, 14)])])
    }
  }else if(a == guide[15]){# O (Pyrrolysine)
    return(nonambigs[9]) # K (Lysine)
  }else if(a == guide[21]){# U (Selenocysteine)
    return(nonambigs[2]) #C )(Cysteine)
  }else if(a == guide[27] | a == guide[28]){
    return(as.raw(0)) # note NULL cant be used with sapply, list cant be simplified.
  }else stop("invalid byte for class 'AAbin'")
}

# Convert amino acid sequence to integers
.encodeAA <- function(x, arity = 20, probs = NULL, random = TRUE, na.rm = FALSE){
  # x is a vector in AAbin format, possibly containing ambiguties
  # arity is an integer, either 20, 22, 24, 26, 27, or 6 (Dayhoff6 compression)
  if(is.null(probs)) probs <- rep(0.05, 20)
  if(arity == 20){
    fun20 <- function(v, probs, random, na.rm){
      ambigs <- !(v %in% as.raw((65:89)[-c(2, 10, 15, 21, 24)]))
      if(any(ambigs)) v[ambigs] <- sapply(unclass(v[ambigs]), .disambiguateAA, probs, random)
      v <- v[v != as.raw(0)]
      bits <- (65:89)[-c(2, 10, 15, 21, 24)]
      ints <- 0:19
      res <- ints[match(as.numeric(v), bits)]
      attributes(res) <- attributes(unclass(v))
      if(na.rm) if(any(is.na(res))) res <- res[!is.na(res)]
      return(res)
    }
    if(is.list(x)){
      res <- lapply(x, fun20, probs, random, na.rm)
    }else{
      res <- fun20(x, probs, random, na.rm)
    }
    return(res)
  }else if(arity == 22){
    # for use with Gonnet matrix
    # return order "C" "S" "T" "P" "A" "G" "N" "D" "E" "Q" "H" "R"
    # "K" "M" "I" "L" "V" "F" "Y" "W" "X" "*"
    # Ambig codes B, J and Z, special codes O and Z, are returned as 20 (X),
    # and gaps are returned as NA
    fun22 <- function(v){
      bits <- c(67, 83, 84, 80, 65, 71, 78, 68, 69, 81, 72, 82, 75, 77,
                73, 76, 86, 70, 89, 87, 88, 42, 66, 74, 79, 85, 90)
      ints <- c(0:21, rep(20, 5))
      res <- ints[match(as.numeric(v), bits)]
      attributes(res) <- attributes(unclass(v))
      isnares <- is.na(res) #placeholder for ambig treatment
      if(na.rm) if(any(isnares)) res <- res[!isnares]
      return(res)
    }
    res <- if(is.list(x)) lapply(x, fun22) else fun22(x)
    return(res)
  }else if(arity == 24){
    # for use with PAM and BLOSUM matrices
    # return order "A" "R" "N" "D" "C" "Q" "E" "G" "H" "I" "L" "K" "M"
    # "F" "P" "S" "T" "W" "Y" "V" "B" "Z" "X" "*"
    # Ambig code J, and special codes O and U are returned as 22 (X).
    # Gaps are returned as NA or removed if na.rm = T
    fun24 <- function(v){
      bits <- c(65, 82, 78, 68, 67, 81, 69, 71, 72, 73, 76, 75,
                77, 70, 80, 83, 84, 87, 89, 86, 66, 90, 88, 42, 74, 79, 85)
      ints <- c(0:23, 22, 22, 22)
      res <- ints[match(as.numeric(v), bits)]
      attributes(res) <- attributes(unclass(v))
      isnares <- is.na(res)
      if(na.rm) if(any(isnares)) res <- res[!isnares]
      return(res)
    }
    res <- if(is.list(x)) lapply(x, fun24) else fun24(x)
    return(res)
  }else if(arity == 26){
    ### return order = LETTERS
    fun26 <- function(v){
      res <- as.integer(v) - 65
      res[res < 0 | res > 25] <- NA
      attributes(res) <- attributes(unclass(v))
      if(na.rm) if(any(is.na(res))) res <- res[!is.na(res)]
      return(res)
    }
    res <- if(is.list(x)) lapply(x, fun26) else fun26(x)
    return(res)
  }else if(arity == 27){
    ### return order ACDEFGHIKLMNPQRSTVWY, X, BJZ, OU, *
    ### for input into .probAA
    fun27 <- function(v){
      bits <- c(65, 67, 68, 69, 70, 71, 72, 73, 75, 76, 77, 78, 80,
                81, 82, 83, 84, 86, 87, 89, 88, 66, 74, 90,  79,85, 42)
      ints <- c(0:26)
      res <- ints[match(as.numeric(v), bits)]
      attributes(res) <- attributes(unclass(v))
      isnares <- is.na(res)
      if(na.rm) if(any(isnares)) res <- res[!isnares]
      return(res)
    }
    res <- if(is.list(x)) lapply(x, fun27) else fun27(x)
    return(res)
  }else if(arity == 6){
    fun6 <- function(v){
      v <- unclass(v)
      bits <- 65:90
      ints <- c(0, 2, 1, 2, 2, 3, 0, 4, 5, 5, 4, 5, 5, 2, 4, 0, 2, 4, 0, 0, 1, 5, 3, NA, 3, 2)
      res <- ints[match(as.numeric(v), bits)]
      attributes(res) <- attributes(v)
      if(na.rm) if(any(is.na(res))) res <- res[!is.na(res)]
      return(res)
    }
    res <- if(is.list(x)) lapply(x, fun6) else fun6(x)
    return(res)
  }else stop("invalid 'arity' argument")
}

# Convert character sequence to integers
.encodeCH <- function(x, residues, na.rm = FALSE){
  fun <- function(v, residues, na.rm = FALSE){
    ints <- seq(0, length(residues) - 1)
    res <- ints[match(v, residues)]
    attributes(res) <- attributes(v)
    if(na.rm) res <- res[!is.na(res)]
    return(res)
  }
  if(is.list(x)) lapply(x, fun, residues, na.rm) else fun(x, residues, na.rm)
}

# Decompose an alignment matrix to a list of non-aligned sequences
.unalign <- function(x, gap = "-"){
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  if(is.list(x)){
    if(length(x) == 1){
      tmpname <- names(x)
      x <- x[[1]]
      if(is.null(dim(x))){
        x <- matrix(x, nrow = 1)
        rownames(x) <- tmpname
      }
    }
  }
  res <- vector(mode = "list", length = nrow(x))
  for(i in 1:nrow(x)) res[[i]] <- x[i, x[i, ] != gap, drop = TRUE]
  if(AA){
    res <- lapply(res, unclass)
    class(res) <- "AAbin"
  }else if(DNA){
    res <- lapply(res, unclass)
    class(res) <- "DNAbin"
  }
  names(res) <- rownames(x)
  return(res)
}

# Find MD5 hash
.digest <- function(x){
  if(.isDNA(x)){
    x <- .dna2char(x)
  }else if(.isAA(x)){
    x <- .aa2char(x)
  }else if(is.list(x)){
    if(length(x) == 0) return(NULL)
    if(length(x[[1]] > 1)){
      x <- vapply(x, paste0, "", collapse = "")
    }else{
      x <- vapply(x, head, "", 1)
    }
  }
  if(mode(x) != "character") stop("sequences must be in raw or character format\n")
  tmpnames <- names(x)
  res <- openssl::md5(x)
  names(res) <- tmpnames
  return(res)
}

# Find re-replication indices
.point <- function(h){
  uh <- unique(h)
  pointers <- seq_along(uh)
  names(pointers) <- uh
  unname(pointers[h])
}

# y is a kcount (or other) matrix
# returns two element integer vector giving
.farthest2 <- function(y, k, seqlengths){
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

## y is a kmer count matrix (not row-normalized)
## returns integer giving row of most central seq
.central1 <- function(y, k, seqlengths, maxdist_central = FALSE){
  y <- as.matrix(y) # in case y is a df
  stopifnot(nrow(y) > 0)
  if(nrow(y) == 2){
    res <- 1
    dis <- .kdist(y, from = 0, to = 1, seqlengths = seqlengths, k = k)[1, ]
    attr(res, "maxdist_centroid") <- dis/2
    if(maxdist_central) attr(res, "maxdist_central") <- dis
    return(res)
  }
  means <- apply(y, 2, mean)
  y2 <- rbind(means, y)
  dists <- .kdist(y2, from = 0, to = 1:nrow(y),
                  seqlengths = c(mean(seqlengths), seqlengths),
                  k = k)[1, ]
  res <- which.min(dists)
  attr(res, "maxdist_centroid") <- max(dists)
  if(maxdist_central){
    dists2 <- .kdist(y, from = res - 1, to = seq(0, nrow(y) - 1),
                     seqlengths = seqlengths,
                     k = k)[1, ]
    attr(res, "maxdist_central") <- max(dists2)
  }
  return(res)
}


.dna2char <- function(x){
  cbytes <- as.raw(c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77, 66, 86, 72, 68, 78, 45, 63))
  indices <- c(136, 24, 72, 40, 96, 144, 192, 48, 80 ,160, 112, 224, 176, 208, 240, 4, 2)
  vec <- raw(240)
  vec[indices] <- cbytes
  if(is.list(x)){
    fun <- function(s, vec) rawToChar(vec[as.integer(s)])
    res <- vapply(x, fun, "", vec)
  }else{
    res <- rawToChar(vec[as.integer(x)])
  }
  return(res)
}

.aa2char <- function(x) if(is.list(x)) vapply(x, rawToChar, "") else rawToChar(x)


#' @noRd
.encodekc <- function(xx){  ## kmer count matrix
  if(mode(xx) == "raw") return(xx)
  xx <- round(xx)
  if(ncol(xx) %% 2 == 1) xx <- cbind(xx, 0L)
  xx[xx > 15] <- 15
  dimnames(xx) <- NULL
  dims <- dim(xx)
  xx <- as.raw(xx)
  dim(xx) <- dims
  fun <- function(xxx){ # 2-col int mat
    m1 <- matrix(rawToBits(xxx[, 1]), ncol = 8, byrow = TRUE)[, 1:4]
    m2 <- matrix(rawToBits(xxx[, 2]), ncol = 8, byrow = TRUE)[, 1:4]
    return(packBits(t(cbind(m2, m1))))
  }
  res <- matrix(as.raw(0), nrow = dims[1], ncol = dims[2]/2)
  guide <- split(seq(1, dims[2]), f = rep(seq(1, dims[2]/2), each = 2))
  for(i in seq_along(guide)){
    res[, i] <- fun(xx[, guide[[i]]])
  }
  return(res)
}

#' @noRd
.decodekc <- function(zz){ ## kmer count matrix (max count 15)
  if(mode(zz) != "raw") return(zz)
  dims <- dim(zz)
  fun <- function(zzz){
    # takes a raw vector
    # returns a 2 col matrix
    mymat <- matrix(rawToBits(zzz), ncol = 8L, byrow = TRUE)
    m1 <- m2 <- matrix(as.raw(0L), ncol = 8L, nrow = nrow(mymat))
    m1[, 1:4] <- mymat[, 5:8]
    m2[, 1:4] <- mymat[, 1:4]
    m1 <- as.integer(packBits(t(m1)))
    m2 <- as.integer(packBits(t(m2)))
    return(cbind(m1, m2))
  }
  guide <- split(seq(1, dims[2] * 2), f = rep(seq(1, dims[2]), each = 2))
  out <- matrix(0L, nrow = dims[1], ncol = dims[2] * 2)
  for(i in seq_along(guide)) out[, guide[[i]]] <- fun(zz[, i])
  return(out)
}
