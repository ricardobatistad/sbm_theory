#' Generate true labels
#'
#' @param n number of nodes
#' @param p proportion of nodes in 1st cluster
#' @param s degree of overlap
#' @export
Z_gen <- function(n, p, s) {

  s1 <- as.integer(n*(1 - s)*p)
  s2 <- as.integer(n*(1 - s)*(1 - p))
  s3 <- n - s1 - s2

  Z <- rep(c(0, 1, 2), c(s1, s2, s3))

  return(Z)

}




#' Shannon's entropy
#'
#' Sources:
#'  https://en.wikipedia.org/wiki/Entropy
#' @param X sample
#' @return entropy
#' @export
entropy <- function(X) {

  H <- 0
  for(x in 0:1){

    p_x <- sum(X == x)/length(X)
    H   <- H + if(p_x == 0) 0 else p_x*log2(p_x)

  }

  return(-H)

}



#' Normalized mutual information
#'
#' Sources:
#'  https://www.researchgate.net/post/How_can_i_calculate_Mutual_Information_theory_from_a_simple_dataset
#'  https://en.wikipedia.org/wiki/Mutual_information
#' @param X True labels
#' @param Y Predicted labels
#' @return NMI score
#' @export
NMI <- function(X, Y){

  I <- 0
  n <- length(X)

  for (y in 0:1) {

    for (x in 0:1) {

      p_xy <- sum(X == x & Y == y)/n
      p_x  <- sum(X == x)/n
      p_y  <- sum(Y == y)/n

      I <- I + if(p_xy == 0) 0 else p_xy*log2(p_xy/(p_x*p_y))

    }

  }

  H_X <- entropy(X)
  H_Y <- entropy(Y)

  score <- if (I == 0) 0 else I/sqrt(H_X*H_Y)

  return(score)

}



#' Mean squared error
#' NOTE: Modify when using the latest python notebook (i.e., when pv for s = 0 actually has only 2 cols).
#'
#' @param pv sufficient stats
#' @param Z true labels
#' @return MSE
#' @export
MSE <- function(pv, Z) {

  Z <- apply(cbind((Z == 0), (Z == 1), (Z == 2)), 2, as.integer)

  samps <- sum(pv[1, ])
  pv    <- pv/samps

  return(mean(sqrt(rowSums((pv - Z)^2))^2))

}
