#' @useDynLib PCOR
#' @importFrom Rcpp sourceCpp
NULL

# get the angle between any two vectors in a matrix
sub_mat <- function(mat, r, n, p){

  # get the difference between one and the other vectors
  mat_ctr <- mat - matrix(1,n,1)%*%matrix(mat[r,],1,p)
  mat_ctr <- mat_ctr[-r,]

  # standardize the difference matrix
  mat_std <- mat_ctr/(matrix(sqrt(rowSums(mat_ctr^2)),n-1,1)%*%matrix(1,1,p))

  # get the angle from there difference
  A <- suppressWarnings(acos(eigenMapMatMult(mat_std, t(mat_std))))
  A[is.nan(A)] <- 0

  return(A)
}


# get a group of pcov value from sample
one_pcov <- function(r, X, Y, n, p, q, estimation.method){

  # initialize the result
  out <- c(0,0)

  # get the matrix multiplication between X and Y
  mx <- sub_mat(X, r, n, p)
  my <- sub_mat(Y, r, n, q)
  xy_mm <- eigenMapMatMult(mx, my)
  xy_ele <- mx*my

  # calculate the u-statistics result and v statistics result
  if(estimation.method=="u"){

    
    mmx <- matrix(diag(mx),n-1,1)%*%matrix(1, 1, n-1)*my
    mmy <- matrix(diag(my),n-1,1)%*%matrix(1, 1, n-1)*mx
    
    a <- sum(diag(xy_mm))
    b <- sum(diag(mx)*diag(my))
    sr1 <- a-b
    sr2 <- sum(xy_mm)-a+2*b-sum(mmx)-sum(mmy)
    sr3 <- (sum(mx)-sum(diag(mx)))*(sum(my)-sum(diag(my)))-2*sr1-4*sr2
  
    # calculate the S2 and pcov
    out[1] <- sr3/((n-1)*(n-2)*(n-3)*(n-4))
    out[2] <- sr1/((n-1)*(n-2)) - 2*sr2/((n-1)*(n-2)*(n-3)) + sr3/((n-1)*(n-2)*(n-3)*(n-4))
  
  }else if(estimation.method=="v"){

    sr <- sum(xy_mm)/n^3

    # calculate the S2 and pcov
    out[1] <- sum(mx)*sum(my)/n^4
    out[2] <- sum(xy_ele)/n^2+out[1]-2*sr

  }else{
    stop("The parameter \"estimation.method\" should be \"u\" or \"v\". ")
  }

  return(out)
}


# get all pcov from all the sample
pcov_va <- function(X, Y, estimation.method){

  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  r <- 1:n

  return(apply(sapply(r, FUN = function(i) one_pcov(i, X=X, Y=Y, n=n, p=p, q=q, estimation.method)), 1, mean))
}


# package the three functions above for parallel computing
pcov_para <- function(cl, X, Y, estimation.method){

  sub_mat <- function(mat, r, n, p){
    mat_ctr <- mat - matrix(1,n,1)%*%matrix(mat[r,],1,p)
    mat_ctr <- mat_ctr[-r,]
    mat_std <- mat_ctr/(matrix(sqrt(rowSums(mat_ctr^2)),n-1,1)%*%matrix(1,1,p))
    A <- suppressWarnings(acos(eigenMapMatMult(mat_std, t(mat_std))))
    A[is.nan(A)] <- 0
    return(A)
  }

  one_pcov <- function(r, X, Y, n, p, q, estimation.method){
    out <- c(0,0)
    mx <- sub_mat(X, r, n, p)
    my <- sub_mat(Y, r, n, q)
    xy_mm <- eigenMapMatMult(mx, my)
    xy_ele <- mx*my
    if(estimation.method=="u"){    
	    mmx <- matrix(diag(mx),n-1,1)%*%matrix(1, 1, n-1)*my
        mmy <- matrix(diag(my),n-1,1)%*%matrix(1, 1, n-1)*mx    
        a <- sum(diag(xy_mm))
        b <- sum(diag(mx)*diag(my))
        sr1 <- a-b
        sr2 <- sum(xy_mm)-a+2*b-sum(mmx)-sum(mmy)
        sr3 <- (sum(mx)-sum(diag(mx)))*(sum(my)-sum(diag(my)))-2*sr1-4*sr2
        out[1] <- sr3/((n-1)*(n-2)*(n-3)*(n-4))
        out[2] <- sr1/((n-1)*(n-2)) - 2*sr2/((n-1)*(n-2)*(n-3)) + sr3/((n-1)*(n-2)*(n-3)*(n-4))
    }else if(estimation.method=="v"){
      sr <- sum(xy_mm)/n^3
      out[1] <- sum(mx)*sum(my)/n^4
      out[2] <- sum(xy_ele)/n^2+out[1]-2*sr
    }else{
      stop("The parameter \"estimation.method\" should be \"u\" or \"v\". ")
    }
    return(out)
  }

  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  registerDoParallel(cl)
  return(foreach(r=1:n) %dopar% one_pcov(r=r, X=X, Y=Y, n=n, p=p, q=q, estimation.method=estimation.method))
}

# get the parallel
pcov_va_para <- function(X, Y, estimation.method){
  cl <- makeCluster(detectCores(), type = "SOCK")
  return(apply(matrix(unlist(pcov_para(cl, X, Y, estimation.method)), ncol = 2, byrow = TRUE), 2, mean))
}


#' Projection Covariance
#'
#' Calculate the projection covariance of two random vectors. Two vectors can be with different dimensions, but with equal sample sizes.
#' @param X A numeric matrix, n*p, each row of the matrix is i.i.d generated from one vector
#' @param Y A numeric matrix, n*q, each row of the matrix is i.i.d generated from the other vector
#' @param estimation.method A character, "u" or "v". When "u", it will use U statistics to estimate projection covariance and other statistics. When "v", it will use V statistics instead.
#' @param parallel A logical, TRUE or FALSE. When TRUE, it will compute parallel, which is faster when sample is larger than around 100. Default is FALSE.
#' @return The projection covariance of X and Y
#' @examples X = matrix(rnorm(100*34,1),100,34)
#' @examples Y = matrix(rnorm(100*62,1),100,62)
#' @examples pcov(X, Y, estimation.method="u", parallel=FALSE)
#' @seealso \code{\link{pcor}} \code{\link{pcor.test}}
#' @references L.Zhu, K.Xu, R.Li, W.Zhong(2017). Projection correlation between two random vectors. Biometrika, Volume 104, Pages 829-843. https://doi.org/10.1093/biomet/asx043
#' @export
pcov <- function(X, Y, estimation.method="u", parallel=FALSE){
  if(nrow(X)>150) message("If the number of your sample is larger than 150, we suggest you setting the parameter \"parallel\" TURE, which is faster.")

  # whether the numbers of two sample are equal
  if(nrow(X)!=nrow(Y)) {
    stop("The numbers of row in two matrix should be equal.")
  }

  # whether n,p,q is big enough.
  if(nrow(X)<=3) {
    stop("The dimension of X should larger than 3.")
  }
  if(ncol(X)<=1) {
    stop("The dimension of X should larger than 1.")
  }
  if(ncol(Y)<=1) {
    stop("The dimension of Y should larger than 1.")
  }

  # calculate the pcov value
  if(parallel == FALSE){
    t <- pcov_va(X, Y, estimation.method)[2]
  }else if(parallel == TRUE){
    t <- pcov_va_para(X, Y, estimation.method)[2]
  }else{
    stop("The parameter \"parallel\" should be TRUE or FALSE.")
  }

  if(t<0){
    warning("The square of projection correlation from U statistics estimation is negative, please try v statistics instead.")
    return(0.000)
  }else{
    return(sqrt(t))
  }
}



#' Projection Correlation
#'
#' Calculate the projection correlation of two random vectors. Two vectors can be with different dimensions, but with equal sample sizes.
#' @param X A numeric matrix, n*p, each row of the matrix is i.i.d generated from one vector
#' @param Y A numeric matrix, n*q, each row of the matrix is i.i.d generated from the other vector
#' @param estimation.method A character, "u" or "v". When "u", it will use U statistics to estimate projection covariance and other statistics. When "v", it will use V statistics instead.
#' @param parallel A logical, TRUE or FALSE. When TRUE, it will compute parallel, which is faster when sample is larger than around 100. Default is FALSE.
#' @return The projection correlation of X and Y
#' @examples X = matrix(rnorm(100*34,1),100,34)
#' @examples Y = matrix(rnorm(100*62,1),100,62)
#' @examples pcor(X, Y, estimation.method="u", parallel=FALSE)
#' @seealso \code{\link{pcov}} \code{\link{pcor.test}}
#' @references L.Zhu, K.Xu, R.Li, W.Zhong(2017). Projection correlation between two random vectors. Biometrika, Volume 104, Pages 829-843. https://doi.org/10.1093/biomet/asx043
#' @export
pcor <- function(X, Y, estimation.method="u", parallel=FALSE){
  if(nrow(X)>150) message("If the number of your sample is larger than 150, we suggest you setting the parameter \"parallel\" TURE, which is faster.")

  t1 <- suppressMessages(pcov(X, Y, estimation.method, parallel))
  t2 <- suppressMessages(pcov(X, X, estimation.method, parallel))
  t3 <- suppressMessages(pcov(Y, Y, estimation.method, parallel))

  if(t1<0 | t2<=0 | t3<=0){
    warning("The square of projection correlation from U statistics estimation is negative, please try v statistics instead.")
    return(0)
  }else{
    return(t1/sqrt(t2)/sqrt(t3))
  }
}



# Projection Correlation Test Statistics
#
# Calculate the projection correlation test statistics of two random vectors. Two vectors can be with different dimensions, but with equal sample sizes.
chisq_va <- function(X, Y, estimation.method){

  # whether the numbers of two sample are equal
  if(nrow(X)!=nrow(Y)) {
    stop("The numbers of row in two matrix should be equal.")
  }

  # whether n,p,q is big enough.
  if(nrow(X)<=3) {
    stop("The dimension of X should larger than 3.")
  }
  if(ncol(X)<=1) {
    stop("The dimension of X should larger than 1.")
  }
  if(ncol(Y)<=1) {
    stop("The dimension of Y should larger than 1.")
  }

  # calculate the chisq_va value
  t <- pcov_va(X, Y, estimation.method)

  return(t[2]/(pi^2-t[1]))
}



#' Projection Correlatoin Permutation Test
#'
#' Return the test result of projection correlation test. Test whether two vectors are independnet. Two vectors can be with different dimensions, but with equal sample sizes.
#' @param X A numeric matrix, n*p, each row of the matrix is i.i.d generated from one vector
#' @param Y A numeric matrix, n*q, each row of the matrix is i.i.d generated from the other vector
#' @param estimation.method A character, "u" or "v". When "u", it will use U statistics to estimate projection covariance and other statistics. When "v", it will use V statistics instead.
#' @param times An int, permutation times, Default is 999.
#' @return \code{pcor.value} projection correlation of X and Y
#' @return \code{p.value} the p-value under the null hypothesis two vectors are independent
#' @examples X = matrix(rnorm(10*7,1),10,7)
#' @examples Y = matrix(rnorm(10*6,1),10,6)
#' @examples pcor.test(X, Y, estimation.method="v", times=999)
#' @seealso \code{\link{pcov}} \code{\link{pcor}}
#' @references L.Zhu, K.Xu, R.Li, W.Zhong(2017). Projection correlation between two random vectors. Biometrika, Volume 104, Pages 829-843. https://doi.org/10.1093/biomet/asx043
#' @export
pcor.test <- function(X, Y, estimation.method="u", times=999){
  message("The permutation test may be slow. Please wait patiently.")
  message("Reset the permutation times if want faster.")

  n <- nrow(X)
  value <- chisq_va(X, Y, estimation.method)
  stat.value <- n*value

  pcor.permu <- replicate(times, chisq_va(X[sample(1:n,n),], Y, estimation.method), simplify = TRUE)
  p.value <- 1-mean(value>pcor.permu)

  dname <- paste(deparse(substitute(X)), "and", deparse(substitute(Y)))
  method <- "Projection Correlation Permutation Test of Independence"
  rval <- list(method = method, data.name = dname, stat.value = stat.value, p.value = p.value)

  return(rval)
}


