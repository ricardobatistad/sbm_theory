# SBM Theory
#
# This package provides supporting functions for experiements based on
# Peixoto's non-parametric MCMC method for community recovery.


#' Generate entries of prefernce matrix
#'
#' Creates the full parameter list from a given n, d, p, lambda, and sign (i.e., +/-).
#' Notice that some values are illegal entries (i.e., values lead to probabilities that
#' are negative or greater than 1). This issue is taken care of later.
#' @param n Number of vertices (i.e., nodes) (integer)
#' @param d Average vertex degree (integer)
#' @param p Proportion of vertices belonging to the 2nd group (double)
#' @param lambda Relevant eigenvalue of matrix R (see Lelarge) (double)
#' @param sign Sign (i.e., +/-) associated with b value calculation (character)
#' @return A one-row data.frame with the full parameter set
#' @export
paramGenerator <- function(n, d, p, lambda, sign){

  ## lambda = d*(1 - b)^2  so  b = 1 +- sqrt(lambda/d)   # Lelarge
  b <- numeric()
  if (sign == "pos"){
    b <- 1 + sqrt(lambda/d)
  } else {
    b <- 1 - sqrt(lambda/d)
  }

  ## pa + (1 - p)b = pb + (1 - p)c = 1                   # Lelarge
  a <- (1 - (1 - p)*b)/p
  c <- (1 - p*b)/(1 - p)

  return(data.frame(n, d, p, lambda, a, b, c, sign))

}



#' Generate M matrix (see Lelarge)
#'
#' Notice that some entries may be illegal. This issue is taken care of in script below.
#' @param params Full set of parameters
#' @return Matrix M
#' @export
matrixGenerator <- function(params){

  n <- as.integer(params["n"])
  d <- as.integer(params["d"])
  a <- as.double(params["a"])
  b <- as.double(params["b"])
  c <- as.double(params["c"])

  M <- matrix(c(a, b, b, c), nrow = 2)
  M <- (d/n)*M

  return(M)

}



#' Determine whether a prefernce matrix is viable
#'
#' @param M preference matrix
#' @return whether the graph is viable
#' @export
nullFinder <- function(M){

  illegalProbs <- sum(M < 0 | 1 < M)
  return(illegalProbs == 0)

}



#' Generate graph from SBM
#'
#' @param graphNum graph number
#' @return graph edge list
#' @export
graph_generator <- function(graphNum){

  params <- parameters_100[graphNum, ]

  n      <- as.integer(params["n"])
  d      <- as.integer(params["d"])
  p      <- params["p"]
  lambda <- params["lambda"]
  sign   <- as.character(params["sign"][1,1])

  blockSize <- as.integer(c(n*p, n*(1 - p)))
  M         <- matrixGenerator(params)

  E <- NULL

  if (nullFinder(M)) {

    set.seed(graphNum)
    G <- igraph::sample_sbm(n, pref.matrix = M, block.sizes = blockSize)
    E <- as.matrix(igraph::get.edgelist(G))

  }

  empty    <- if(is.null(E)) "_NULL" else ""

  fileName <- sprintf("graphs/graph_%d_%d_%d_%.2f_%.1f_%s%s.csv",
                      graphNum, n, d, p, lambda, sign, empty)
  utils::write.csv(E, fileName)

  # cat(fileName)

  return(E)

}



#' Produce stats for a given graphs
#'
#' @param params row of parameter matrix
#' @return graph summary df
#' @export
graph_stats <- function(params){

  n <- params$n
  d <- params$d
  p <- params$p

  M <- matrixGenerator(params)

  # Illegal probs
  validIdx <- nullFinder(M)


  # M matrix
  a_scaled = M[1,1]
  b_scaled = M[1,2]
  c_scaled = M[2,2]


  # Average degree per block
  block_0 = n*p
  block_1 = n*(1 - p)

  d_0   = (block_0 - 1)*a_scaled
  d_1   = (block_1 - 1)*c_scaled
  d_out = (n - 1)*b_scaled


  # SBM stats
  threshold <- 2*sqrt(d)  # B = 2

  c_in  <- M[1,1]*p + M[2,2]*(1 - p)
  c_out <- M[1,2]
  e     <- n*(c_in - c_out)

  recov <- NULL
  if (e < 0) {
    recov <- "Unknown"
  } else {
    recov <- if (threshold < e) "Yes" else "No"
  }

  return(data.frame(validIdx,
                    a_scaled, b_scaled, c_scaled,
                    d_0, d_1, d_out,
                    c_in, c_out, e, threshold, recov))

}


##################################################################################
################################ Overlap #########################################
##################################################################################



#' Produce separation (i.e., c) params
#'
#' @param c sepration
#' @param d avg degree of graph
#' @param p cluster sizes
#' @return entries of M matrix
#' @export
separation_params <- function(c, d, p){

  c_max <- 2*d
  c_in  <- p*(c_max + c)
  c_out <- p*(c_max - c)

  return(c("c_in" = c_in, "c_out" = c_out, "c_max" = c_max))

}



#' Produce block sizes
#'
#' @param n num of vertices
#' @param p cluster sizes
#' @param s overlap
#' @return whether the graph is viable
#' @export
blockSizes_gen <- function(n, p, s) {

  blockSize <- NULL

  if (s == 0) {

    block_0 <- n*p
    block_1 <- n*(1-p)

    blockSize <- c("block_0" = block_0, "block_1" = block_1)

  } else {

    block_0 <- as.integer(n*(1 - s)*p)
    block_1 <- as.integer(n*(1 - s)*(1 - p))
    block_2 <- n - block_0 - block_1

    blockSize <- c("block_0" = block_0, "block_1" = block_1, "block_2" = block_2)

  }

  return(blockSize)

}



#' Produce block sizes (normalized)
#'
#' @param blockSize vector with block sizes
#' @return normalized block sizes
#' @export
blockSizes_gen_norm <- function(blockSize) {

  case <- if(length(blockSize) == 3) "default" else "2 groups"

  blockSize_norm <- switch(case, "2 groups" = c(blockSize, "block_2" = 0),
                                 "default"  = blockSize)

  return(blockSize_norm)

}


#' Generate M matrix for overlap
#'
#' @param params Full set of parameters
#' @return Matrix M
#' @export
matrixGenerator_overlap <- function(params){

  n <- params$n
  d <- params$d
  c <- params$c
  p <- params$p
  s <- params$s

  C <- separation_params(c, d, p)
  c_in  <- C["c_in"]
  c_out <- C["c_out"]
  c_max <- C["c_max"]

  M <- NULL
  if (s == 0) {
    M <- matrix(c(c_in, c_out, c_out, c_in), nrow = 2)
  } else {
    M <- matrix(c(c_in, c_out, c_max, c_out, c_in, c_max, c_max, c_max, 2*c_max), nrow = 3)
  }

  return(M/n)

}



#' Generate graph from SBM (overlap)
#'
#' @param graphNum graph number
#' @return graph edge list
#' @export
graph_generator_overlap <- function(graphNum) {

  params <- parameters_overlap_100[graphNum, ]

  M <- matrixGenerator_overlap(params)

  n <- params$n
  d <- params$d
  p <- params$p
  s <- params$s
  c <- params$c
  ks <- params$ks

  blockSize <- blockSizes_gen(n, p, s)

  set.seed(graphNum)
  G <- sample_sbm(n, pref.matrix = M, block.sizes = blockSize)
  E <- as.matrix(get.edgelist(G))

  fileName <- sprintf("graphs_overlap/graph_%d_%d_%.2f_%.1f_%.2f_%.2f.csv",
                      graphNum, n, d, s, c, ks)

  utils::write.csv(E, fileName)

  # cat(fileName)

  return(E)

}



#' Produce stats for a given graphs
#'
#' @param params row of parameter matrix
#' @return graph summary df
#' @export
graph_stats_overlap <- function(params){

  n <- params$n
  p <- params$p
  d <- params$d
  s <- params$s
  c <- params$c

  # block sizes
  blockSize      <- blockSizes_gen(n, p, s)
  B              <- length(blockSize)
  blockSize_norm <- blockSizes_gen_norm(blockSize)

  block_0 <- blockSize_norm["block_0"]
  block_1 <- blockSize_norm["block_1"]
  block_2 <- blockSize_norm["block_2"]


  # Number of edges
  M <- matrixGenerator_overlap(params)

  num_edges = 0
  if (B == 2){

    e_00 <- M[1,1]*(block_0*(block_0 - 1)/2)
    e_11 <- M[2,2]*(block_1*(block_1 - 1)/2)
    e_01 <- M[1,2]*(block_0*block_1)

    num_edges = e_00 + e_11 + e_01

  } else {

    e_00 <- M[1,1]*(block_0*(block_0 - 1)/2)
    e_11 <- M[2,2]*(block_1*(block_1 - 1)/2)
    e_22 <- M[3,3]*(block_2*(block_2 - 1)/2)
    e_01 <- M[1,2]*(block_0*block_1)
    e_02 <- M[1,3]*(block_0*block_2)
    e_12 <- M[2,3]*(block_1*block_2)

    num_edges = e_00 + e_11 + e_22 + e_01 + e_02 + e_12

  }


  # c_in, c_out, c_max
  C <- separation_params(c, d, p)
  c_in  <- C["c_in"]
  c_out <- C["c_out"]
  c_max <- C["c_max"]


  # c_in, c_out, c_max scaled
  c_in_scaled  <- c_in/n
  c_out_scaled <- c_out/n
  c_max_scaled <- c_max/n


  # recoverability
  threshold <- B*sqrt(d)
  e <- n*(c_in_scaled - c_out_scaled)

  recov <- NULL
  if (e < 0) {
    recov <- "Unknown"
  } else {
    recov <- if (threshold < e) "Yes" else "No"
  }

  return(data.frame(B, block_0, block_1, block_2, num_edges,
                    c_in_scaled, c_out_scaled, c_max_scaled,
                    c_in, c_out, c_max,
                    e, threshold, recov))

}



#' #' Produce stats for a given graph
#' #'
#' #' @param graphNum graph index
#' #' @return summary stats (mostly reverse-engineering the parameters)
#' #' @export
#' stats_overlap_samp <- function(graphNum){
#'
#'   # Produce graph
#'   E <- graph_generator_overlap(graphNum)
#'   G <- graph_from_edgelist(E, directed = FALSE)
#'   A <- as.matrix(as_adjacency_matrix(G))
#'
#'   params <- parameters_overlap_100[graphNum, ]
#'   n <- params$n
#'   d <- params$d
#'   p <- params$p
#'   s <- params$s
#'
#'   # edges per block
#'   case <- if(s == 1) "1 group" else if(s == 0) "2 groups" else "default"
#'
#'   blockSize      <- blockSizes_gen(n, p, s)
#'   blockSize_norm <- blockSizes_gen_norm(blockSize)
#'
#'   size_0 <- blockSize_norm["block0"]
#'   size_1 <- blockSize_norm["block1"]
#'   size_2 <- blockSize_norm["block2"]
#'
#'
#'
#'   # block 0
#'   block_0 <- switch(case, "1 group"  = NULL,
#'                           "2 groups" = A[1:size_0, 1:size_0],
#'                           "default"  = A[1:size_0, 1:size_0])
#'
#'   block_0_numEdges     <- sum(upper.tri(block_0, diag = FALSE))
#'   block_0_numEdges_avg <- mean(rowSums(block_0))
#'
#'
#'   # block 1
#'   block_1 <- switch(case, "1 group"  = NULL,
#'                           "2 groups" = A[(size_0 + 1):n, (size_0 + 1):n],  # could've also been calculated as in the 'default' case
#'                           "default"  = A[(size_0 + 1):(size_0 + size_1), (size_0 + 1):(size_0 + size_1)])
#'
#'   block_1_numEdges     <- sum(upper.tri(block_1, diag = FALSE))
#'   block_1_numEdges_avg <- mean(rowSums(block_1))
#'
#'
#'   # block 2
#'   block_2 <- switch(case, "1 group"  = A,
#'                           "2 groups" = NULL,
#'                           "default"  = A[(size_0 + size_1 + 1):n, (size_0 + size_1 + 1):n])
#'
#'   block_2_numEdges     <- sum(upper.tri(block_2, diag = FALSE))
#'   block_2_numEdges_avg <- mean(rowSums(block_2))
#'
#'
#'   # block out 01-10
#'   block_out_01 <- switch(case, "1 group"  = NULL,
#'                                "2 groups" = A[1:size_0, (size_0 + 1):n],   # could've also been calculated as in the 'default' case
#'                                "default"  = A[1:size_0, (size_0 + 1):(size_0 + size_1)])
#'
#'   block_out_01_numEdges     <- sum(block_out_01)   # recall we're recreating upper triangular of A!
#'   block_out_01_numEdges_avg <- mean(c(rowSums(block_out_01), colSums(block_out_01)))
#'
#'
#'   # block out 02-20
#'   block_out_02 <- switch(case, "1 group"  = NULL,
#'                                "2 groups" = NULL,
#'                                "default"  = A[1:size_0, (size_0 + size_1 + 1):n])
#'
#'   block_out_02_numEdges     <- sum(block_out_02)   # recall we're recreating upper triangular of A!
#'   block_out_02_numEdges_avg <- mean(c(rowSums(block_out_02), colSums(block_out_02)))
#'
#'
#'   # block out 12-21
#'   block_out_12 <- switch(case, "1 group"  = NULL,
#'                                "2 groups" = NULL,
#'                                "default"  = A[(size_0 + 1):(size_0 + size_1), (size_0 + size_1 + 1):n])
#'
#'   block_out_12_numEdges     <- sum(block_out_12)   # recall we're recreating upper triangular of A!
#'   block_out_12_numEdges_avg <- mean(c(rowSums(block_out_12), colSums(block_out_12)))
#'
#'
#'   c_in         <- block_0_edges_avg + block_1_edges_avg + block_2_edges_avg
#'   c_out        <- block_out_01_numEdges_avg + block_out_02_numEdges_avg + block_out_12_numEdges_avg
#'   c_in_scaled  <- c_in/n
#'   c_out_scaled <- c_out/n
#'
#'
#'   # average degree
#'
#'   ## Easy
#'   numEdges    <- dim(E)[1]
#'   d_samp_easy <- (2*numEdges)/n
#'
#'   ## Medium
#'   numEdges_med <- sum(upper.tri(A, diag = FALSE))
#'   d_samp_med   <- mean(rowSums(A))
#'
#'   ## Hard
#'   numEdges_hard <- block_0_numEdges + block_1_numEdges + block_2_numEdges +
#'                    block_out_01_numEdges + block_out_02_numEdges + block_out_12_numEdges
#'   d_samp_hard   <- block_0_edges_avg*size_0 + block_1_edges_avg*size_1 + block_2_edges_avg*size_2 +
#'                    block_out_01_numEdges_avg*size_0*size_1 + block_out_02_numEdges_avg*size_0*size_2 +
#'                    block_out_12_numEdges_avg*size_1*size_2
#'
#'   ## Check
#'   print(c(d, d_samp_easy, d_samp_med, d_samp_hard))
#'
#'   # c_max
#'   c_max        <- 2*d_samp_easy
#'   c_max_scaled <- c_max/n
#'
#'   # e
#'   e <- n*(c_in_scaled - c_out_scaled)
#'
#'   sampStats <- data.frame(numEdges, d_samp = d_samp_easy,
#'                           c_in_samp = c_in, c_out_samp = c_out, c_max_samp = c_max,
#'                           e_samp = e,
#'                           c_in_scaled_samp = c_in_scaled, c_out_scaled_samp = c_out_scaled, c_max_scaled_samp = c_max_scaled,
#'                           B0_deg = block_0_numEdges_avg, B1_deg = block_1_numEdges_avg, B2_deg = block_2_numEdges_avg,
#'                           B01_deg = block_out_01_numEdges_avg, B02_deg = block_out_02_numEdges_avg,
#'                           B12_deg = block_out_12_numEdges_avg)
#'
#'   return(cbind(params, sampStats))
#'
#' }


##################################################################################
################################   3x3   #########################################
##################################################################################


#' Generate entries of prefernce matrix
#'
#' Creates the full parameter list from a given n, d, p, lambda, and sign (i.e., +/-).
#' Notice that some values are illegal entries (i.e., values lead to probabilities that
#' are negative or greater than 1). This issue is taken care of later.
#' @param n Number of vertices (i.e., nodes) (integer)
#' @param d Average vertex degree (integer)
#' @param p Proportion of vertices belonging to each group (vector)
#' @param lambda Relevant eigenvalue of matrix R (see Lelarge) (double)
#' @param sign Sign (i.e., +/-) associated with b value calculation (character)
#' @return A one-row data.frame with the full parameter set
#' @export
paramGenerator_3x3 <- function(n, d, p, lambda){

  ## lambda = d*(1 - b)^2  so  b = 1 +- sqrt(lambda/d)   # Lelarge
  b <- numeric()
  if (sign == "pos"){
    b <- 1 + sqrt(lambda/d)
  } else {
    b <- 1 - sqrt(lambda/d)
  }

  ## pa + (1 - p)b = pb + (1 - p)c = 1                   # Lelarge
  a <- (1 - (1 - p)*b)/p
  c <- (1 - p*b)/(1 - p)

  return(data.frame(n, d, p, lambda, a, b, c, sign))

}



#' Generate graph from SBM (overlap)
#'
#' @param graphNum graph number
#' @return graph edge list
#' @export
graph_generator_3x3 <- function(graphNum) {



}


