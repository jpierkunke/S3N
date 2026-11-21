# inputs:
# - BRISC_Out : an object of class BRISC_Out, obtained as an output of BRISC_estimation
# - n_boot: number of bootstrap samples (Default value is 100)
# - h: number of cores to be used in parallel computing setup for bootstrap samples
#     - If h = 1, there is no parallelization (Default value is 1)
# - n_omp: number of threads to be used
#     - value can be more than 1 if source code is compiled with OpenMP support
#     - Default is 1.
# - init: keyword that specifies initialization scheme to be used
#     - "Initial" for initialization of parameter values for bootstrap samples 
#       with initial values used in BRISC_estimate
#     - "Estimate" for initialization of parameter values for bootstrap samples 
#       with estimated values of parameters in BRISC_estimate respectively.
# - verbose:
#     - if TRUE, model specifications along with information regarding OpenMP support
#       and progress of the algorithm is printed to the screen. Otherwise, nothing is
#       printed to the screen. Default value is TRUE.
# - nugget_status: 
#     - if nugget_status = 0, tau.sq is fixed to 0, 
#     - if nugget_status = 1 tau.sq is estimated
#     - Default value is 1.
#
# outputs:
# - boot.Theta:
#     - estimates of spatial covariance parameters corresponding to bootstrap samples
# - boot.Beta: 
#     - estimates of beta corresponding to bootstrap samples
# - confidence.interval: 
#     - confidence intervals corresponding to the parameters
# - boot.time: 
#     - time (in seconds) required to perform the bootstrapping after 
#       preprocessing data in R, reported using proc.time()
BRISC_bootstrap <- function(BRISC_Out, n_boot = 100, h = 1, n_omp = 1, 
                            init = "Initial", verbose = TRUE, nugget_status = 1){

  if(missing(BRISC_Out)){stop("error: BRISC_bootstrap expects BRISC_Out\n")}

  if(nugget_status == 0){fix_nugget = 0}
  if(nugget_status == 1){fix_nugget = 1}
  X <- BRISC_Out$X
  n.omp.threads <- as.integer(n_omp)
  n.neighbors <- BRISC_Out$n.neighbors
  eps <- BRISC_Out$eps
  cov.model <- BRISC_Out$cov.model
  p  <- ncol(X)
  n <- nrow(X)


  storage.mode(X) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"
  storage.mode(n.neighbors) <- "integer"
  storage.mode(n.omp.threads) <- "integer"
  storage.mode(eps) <- "double"

  cov.model.names <- c("exponential","spherical","matern","gaussian")
  cov.model.indx <- which(cov.model == cov.model.names) - 1
  storage.mode(cov.model.indx) <- "integer"

  cov.model <- BRISC_Out$cov.model
  B =  BRISC_Out$BRISC_Object$B
  efF = BRISC_Out$BRISC_Object$F
  norm.residual = BRISC_Out$BRISC_Object$norm.residual
  Xbeta = BRISC_Out$BRISC_Object$Xbeta
  d = BRISC_Out$BRISC_Object$d
  nnWt = BRISC_Out$BRISC_Object$nnWt
  nnIndx = BRISC_Out$BRISC_Object$nnIndx
  nnIndxLU = BRISC_Out$BRISC_Object$nnIndxLU
  D = BRISC_Out$BRISC_Object$D
  Dwt = BRISC_Out$BRISC_Object$Dwt
  CIndx = BRISC_Out$BRISC_Object$CIndx
  Length.D = BRISC_Out$BRISC_Object$Length.D

  # these three lines are new once I changed neighbor.R to handle stream distances
  storage.mode(nnIndx) <- "integer"
  storage.mode(nnIndxLU) <- "integer"
  storage.mode(CIndx) <- "integer"

  if(init == "Initial"){
    if(cov.model == "matern") {
      theta_boot_init <- c(
        BRISC_Out$init[2]/BRISC_Out$init[1], 
        BRISC_Out$init[3], 
        BRISC_Out$init[4])
    } else {
      theta_boot_init <- c(
        BRISC_Out$init[2]/BRISC_Out$init[1], 
        BRISC_Out$init[3])
    }
  }
  if(init == "Estimate"){
    if(cov.model == "matern") {
      theta_boot_init <- c(
        BRISC_Out$Theta[2]/BRISC_Out$Theta[1], 
        BRISC_Out$Theta[3], 
        BRISC_Out$Theta[4])
    } else {
      theta_boot_init <- c(
        BRISC_Out$Theta[2]/BRISC_Out$Theta[1], 
        BRISC_Out$Theta[3])}
  }

  theta_boot_init <- sqrt(theta_boot_init)

  p3 <- proc.time()

  # parallelization with h cores and n.omp.threads threads
  if(h > 1){
    cl <- makeCluster(h)
    tryCatch({
      clusterCall(cl, .libPaths, .libPaths())
      clusterEvalQ(cl, loadNamespace("S3N"))
      clusterExport(cl=cl, 
                    varlist=c("norm.residual", 
                              "X", "B", "efF", 
                              "Xbeta", "D", "d", 
                              "nnIndx", "nnIndxLU",
                              "CIndx", "nnWt", "Dwt",
                              "n", "p", 
                              "n.neighbors", "theta_boot_init", 
                              "cov.model.indx", "Length.D",
                              "n.omp.threads", "bootstrap_brisc", 
                              "eps", "fix_nugget"),
                    envir=environment())
      if(verbose == TRUE){
        cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tBootstrap Progress"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
        pboptions(type = "txt", char = "=")
        result <- pblapply(
          1:n_boot,
          bootstrap_brisc,
          norm.residual, 
          X, B, efF, Xbeta, D, d, 
          nnIndx, nnIndxLU, CIndx, nnWt, Dwt,
          n, p, n.neighbors, 
          theta_boot_init, cov.model.indx, 
          Length.D, n.omp.threads, eps, fix_nugget,
          cl = cl)
        } else { # not verbose
        result <- parLapply(
          cl,
          1:n_boot,
          bootstrap_brisc,
          norm.residual, 
          X, B, efF, Xbeta, D, d,
          nnIndx, nnIndxLU, CIndx, nnWt, Dwt,
          n, p, n.neighbors, 
          theta_boot_init, cov.model.indx, 
          Length.D, n.omp.threads, eps, fix_nugget)
        }
    }, finally = stopCluster(cl))
  }
  if(h == 1){ # no parallelization (1 core)
    if(verbose == TRUE){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tBootstrap Progress"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
      pboptions(type = "txt", char = "=")
      # print("Right before bootstrap_brisc\n")
      # print("All fine till here 1") # this is what it's calling; fine till here, so issue is in C code
      result <- pblapply( # lapply with progress bar (pb)
        1:n_boot, 
        bootstrap_brisc, 
        norm.residual, 
        X, B, efF, Xbeta, D, d, 
        nnIndx, nnIndxLU, CIndx, nnWt, Dwt, 
        n, p, n.neighbors, 
        theta_boot_init, cov.model.indx, 
        Length.D, n.omp.threads, eps, fix_nugget)
    } else {
      # print("All fine till here 2")
      result <- lapply(
        1:n_boot,
        bootstrap_brisc,
        norm.residual, 
        X, B, efF, Xbeta, D, d, 
        nnIndx, nnIndxLU, CIndx, nnWt, Dwt, 
        n, p, n.neighbors, 
        theta_boot_init, cov.model.indx, 
        Length.D, n.omp.threads, eps, fix_nugget)
    }
  }

  p4 <- proc.time()

  result_table = arrange(result)

  estimate <- c(BRISC_Out$Beta, BRISC_Out$Theta)
  result_CI <- matrix(0,2,length(estimate))

  for(i in 1:length(estimate)){
    result_CI[,i] <- 2*estimate[i] - quantile(result_table[,i], c(.975,.025))
  }

  result_list <- list()


  result_list$boot.Theta <- result_table[,(length(BRISC_Out$Beta) + 1):dim(result_table)[2]]
  if (cov.model != "matern") {colnames(result_list$boot.Theta) <- c("sigma.sq", "tau.sq", "phi")}
  if (cov.model == "matern") {colnames(result_list$boot.Theta) <- c("sigma.sq", "tau.sq", "phi", "nu")}
  result_list$boot.Beta <- as.matrix(result_table[,1:length(BRISC_Out$Beta)])
  colnames(result_list$boot.Beta) <- rep(0, length(BRISC_Out$Beta))
  for(i in 1:length(BRISC_Out$Beta)){
    name_beta <- paste0("beta_",i)
    colnames(result_list$boot.Beta)[i] <- name_beta
  }
  result_list$confidence.interval <- cbind(result_CI[,1:length(BRISC_Out$Beta)],
                                           pmax(
                                             result_CI[,(length(BRISC_Out$Beta) + 1):dim(result_table)[2]],
                                             0*result_CI[,(length(BRISC_Out$Beta) + 1):dim(result_table)[2]]))
  if (cov.model != "matern")  {colnames(result_list$confidence.interval)[(length(BRISC_Out$Beta) + 1):dim(result_table)[2]] <-
    c("sigma.sq", "tau.sq", "phi")}
  if (cov.model == "matern")  {colnames(result_list$confidence.interval)[(length(BRISC_Out$Beta) + 1):dim(result_table)[2]] <-
    c("sigma.sq", "tau.sq", "phi", "nu")}
  for(i in 1:length(BRISC_Out$Beta)){
    name_beta <- paste0("beta_",i)
    colnames(result_list$confidence.interval)[i] <- name_beta
  }
  result_list$boot.time = p4 - p3
  result_list
}
