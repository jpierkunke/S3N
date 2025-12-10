# let q = number of prediction points, 
# n = number of observation points,
# m = number of obs neighbors for each pred point
# the neighbor variable must be a list of three objects:
# - nnIndx (q x m): ith row corresponds to prediction point i and consists of 
#   the indices of point i's m nearest obs neighbors
#   for prediction locations i that are both obs and preds points, this consists
#   of the m nearest neighbors including point i itself, regardless of how
#   many neighbors this point was conditioned on in the estimation stage
# - d (q x m): ith row corresponds to prediction point i and consists of 
#   the distances between point i and each of its m nearest obs neighbors
# - nnWght (q x m): ith row corresponds to prediction point i and consists of 
#   the spatial weights for the tail-up covariance between  point i and each of 
#   its m nearest obs neighbors
#
# prediction also requires obsobs_dist and obsobs_wt in order to look up distances
# and weights for each neighbor-neighbor pair among the neighbors for each
# prediction point

# intercept must be included in X.0 along with covariates
BRISC_prediction_stream <- function(BRISC_Out, coords.0, neighbor, 
                                    obsobs_dist, obsobs_wt,
                                    X.0 = NULL, n_omp = 1, verbose = TRUE){

  if(missing(BRISC_Out)){stop("error: BRISC_prediction expects BRISC_Out\n")}
  
  # number of prediction point locations
  n.0 <- nrow(coords.0)
  # if no X.0 provided, include an intercept but no covariates
  if(is.null(X.0)){
      X.0 <- matrix(1, nrow = n.0, ncol = 1)
  }

  # observed covariates, response
  X <- BRISC_Out$X
  y <- BRISC_Out$y
  n.neighbors <- BRISC_Out$n.neighbors
  # observed locations
  coords <- BRISC_Out$coords
  
  n.omp.threads <- as.integer(n_omp)
  n <- nrow(X)
  p <- ncol(X)
  Theta <- as.matrix(BRISC_Out$Theta)
  dim(Theta) <- c(length(BRISC_Out$Theta),1)
  Beta <- as.matrix(BRISC_Out$Beta)
  dim(Beta) <- c(length(BRISC_Out$Beta),1)
  
  cov.model <- BRISC_Out$cov.model
  cov.model.names <- c("exponential","spherical","matern","gaussian")
  cov.model.indx <- which(cov.model == cov.model.names) - 1
  
  ##check X.0 and coords.0
  if(missing(X.0)){stop("error: X.0 must be specified\n")}
  if(!any(is.data.frame(X.0), is.matrix(X.0))){stop("error: X.0 must be a data.frame or matrix\n")}
  if(ncol(X.0) != ncol(X)){ stop(paste("error: X.0 must have ",p," columns\n"))}

  if(missing(coords.0)){stop("error: coords.0 must be specified\n")}
  if(!any(is.data.frame(coords.0), is.matrix(coords.0))){stop("error: coords.0 must be a data.frame or matrix\n")}
  # if(!ncol(coords.0) == 2){stop("error: coords.0 must have two columns\n")}

  # if(nrow(nn.indx) != nrow(nn.dists) | nrow(nn.indx) != nrow(X.0)){
  #   stop("error: nn.indx, nn.dists, and X.0 must have the same number of rows, equal to the number of prediction locations")
  # }
  
  q <- nrow(X.0)
  
  nnIndx0 = neighbor$nnIndx
  nnDist0 = neighbor$nnDist
  nnWght0 = neighbor$nnWght
  
  storage.mode(X) <- "double"
  storage.mode(y) <- "double"
  storage.mode(coords) <- "double"
  storage.mode(n) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(n.neighbors) <- "integer"
  storage.mode(X.0) <- "double"
  storage.mode(coords.0) <- "double"
  storage.mode(q) <- "integer"
  storage.mode(Beta) <- "double"
  storage.mode(Theta) <- "double"
  storage.mode(cov.model.indx) <- "integer"
  storage.mode(n.omp.threads) <- "integer"
  storage.mode(verbose) <- "integer"
  
  storage.mode(nnIndx0) <- "integer"
  storage.mode(nnDist0) <- "double"
  storage.mode(nnWght0) <- "double"
  storage.mode(obsobs_dist) <- "double"
  storage.mode(obsobs_wt) <- "double"
  
  p5 <- proc.time()
  
  out <- .Call("BRISC_predictioncpp_s", X, y,
               coords, n, p, n.neighbors, X.0, 
               coords.0, q, 
               #nn.indx, nn.dists, pwdist_matrix,
               nnIndx0, nnDist0, nnWght0,
               obsobs_dist, obsobs_wt,
               Beta, Theta, cov.model.indx,
               n.omp.threads, verbose, PACKAGE = "S3N")

  p6 <- proc.time()
  
  output <- list()
  output$prediction <- out$p.y.0[,1]

  result_new <- matrix(0,length(output$prediction),2)
  result_new[,1] <- out$p.y.0[,1] - 1.96 * out$var.y.0[,1]
  result_new[,2] <- out$p.y.0[,1] + 1.96 * out$var.y.0[,1]

  output$prediction.ci <- result_new
  output$prediction.time <- p6 - p5
  output
}

BRISC_prediction <- function(BRISC_Out, coords.0, X.0 = NULL, 
                             n_omp = 1, verbose = TRUE, tol = 12){
  
  
  if(missing(BRISC_Out)){stop("error: BRISC_prediction expects BRISC_Out\n")}
  
  n.0 <- nrow(coords.0)
  if(is.null(X.0)){
    X.0 <- matrix(1, nrow = n.0, ncol = 1)
  }
  
  coords.0 <- round(coords.0, tol)
  X.0 <- round(X.0, tol)
  
  X <- BRISC_Out$X
  y <- BRISC_Out$y
  n.neighbors <- BRISC_Out$n.neighbors
  coords <- BRISC_Out$coords
  
  n.omp.threads <- as.integer(n_omp)
  n <- nrow(X)
  p <- ncol(X)
  Theta <- as.matrix(BRISC_Out$Theta)
  dim(Theta) <- c(length(BRISC_Out$Theta),1)
  Beta <- as.matrix(BRISC_Out$Beta)
  dim(Beta) <- c(length(BRISC_Out$Beta),1)
  
  cov.model <- BRISC_Out$cov.model
  cov.model.names <- c("exponential","spherical","matern","gaussian")
  cov.model.indx <- which(cov.model == cov.model.names) - 1
  
  
  ##check X.0 and coords.0
  if(missing(X.0)){stop("error: X.0 must be specified\n")}
  if(!any(is.data.frame(X.0), is.matrix(X.0))){
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  if(ncol(X.0) != ncol(X)){ stop(paste("error: X.0 must have ",p," columns\n"))}
  
  if(missing(coords.0)){stop("error: coords.0 must be specified\n")}
  if(!any(is.data.frame(coords.0), is.matrix(coords.0))){
    stop("error: coords.0 must be a data.frame or matrix\n")
    }
  if(!ncol(coords.0) == 2){stop("error: coords.0 must have two columns\n")}
  
  q <- nrow(X.0)
  
  ##get nn indx
  if(verbose == TRUE){
    cat(paste(("----------------------------------------"), collapse="   "), "\n")
    cat(paste(("\tGetting nearest neighbors\n\tfor prediction locations"), 
              collapse="   "), "\n")
  }
  # use nn2 from RANN package to find k nearest neighbors from obs (coords) 
  # for each preds point (coords.0)
  # - doesn't work for stream distance
  nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1
  # from ?nn2:
  # returns a list of length 2
  # nn.idx: N x k integer matrix returning the near neighbour indices
  # nn.dists: N x k matrix returning the near neighbour Euclidean distances
  
  # class(nn.indx.0)
  # [1] "matrix" "array" 
  #
  # dim(nn.indx.0) # no. of pred points x no. of neighbors (q x m)
  # [1] 168657     10
  #
  # str(nn.indx.0)
  # num [1:168657, 1:10] 257 306 257 306 306 306 306 306 306 306 ...
  #
  # head(nn.indx.0)
  #      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
  # [1,]  257  306  303  290  236  238  320  286  336   221
  # [2,]  306  257  303  290  236  320  238  286  336   354
  # [3,]  257  306  303  290  236  320  238  286  336   354
  # [4,]  306  257  303  290  320  236  238  286  336   354
  # [5,]  306  257  303  290  320  238  236  336  286   354
  # [6,]  306  257  290  303  320  238  336  354  346   236
  
  storage.mode(X) <- "double"
  storage.mode(y) <- "double"
  storage.mode(coords) <- "double"
  storage.mode(n) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(n.neighbors) <- "integer"
  storage.mode(X.0) <- "double"
  storage.mode(coords.0) <- "double"
  storage.mode(q) <- "integer"
  storage.mode(Beta) <- "double"
  storage.mode(Theta) <- "double"
  storage.mode(cov.model.indx) <- "integer"
  storage.mode(nn.indx.0) <- "integer"
  storage.mode(n.omp.threads) <- "integer"
  storage.mode(verbose) <- "integer"
  
  p5 <- proc.time()
  
  out <- .Call("BRISC_predictioncpp", X, y, coords, n, p, n.neighbors, X.0, 
               coords.0, q, nn.indx.0, Beta, Theta, cov.model.indx, 
               n.omp.threads, verbose, PACKAGE = "S3N")
  p6 <- proc.time()
  
  output <- list()
  output$prediction <- out$p.y.0[,1]
  
  result_new <- matrix(0,length(output$prediction),2)
  result_new[,1] <- out$p.y.0[,1] - 1.96 * out$var.y.0[,1]
  result_new[,2] <- out$p.y.0[,1] + 1.96 * out$var.y.0[,1]
  
  output$prediction.ci <- result_new
  output$prediction.time <- p6 - p5
  output
}

