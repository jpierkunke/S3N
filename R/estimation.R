BRISC_estimation_stream <- function(coords, y, x, neighbor,
                                    sigma.sq = 1, tau.sq = 0.1, phi = 1, nu = 1.5, 
                                    n_omp = 1, cov.model = "exponential", 
                                    search.type = "tree", stabilization = NULL,
                                    pred.stabilization = 1e-8, verbose = TRUE, 
                                    eps = 2e-05, nugget_status = 1
){
  if(nrow(y) != nrow(x) | nrow(y) != nrow(coords)){
    stop("Error: x, y, and coords must all have the same number of rows, equal to the number of observation locations")
  }
  n <- nrow(y)
  p <- ncol(x)
  if(nugget_status == 0){fix_nugget = 0}
  if(nugget_status == 1){fix_nugget = 1}
  if(is.matrix(y) & ncol(as.matrix(y)) > 1){
    stop("error: BRISC handles only univariate response. You have used a y (response) with more than one column.
          The ordering of inputs x (covariates) and y (response) in BRISC_estimation has been changed BRISC 0.2.0 onwards.
          Did you mean x (covariates) instead of y (response) by the multivariate input?
          Please check the new documentation with ?BRISC_estimation.")
  }
  
  if(tau.sq < 0 ){stop("error: tau.sq must be non-negative")}
  if(sigma.sq < 0 ){stop("error: sigma.sq must be non-negative")}
  if(phi < 0 ){stop("error: phi must be non-negative")}
  if(nu < 0 ){stop("error: nu must be non-negative")}
  
  ##Covariance model
  cov.model.names <- c("exponential","spherical","matern","gaussian")
  cov.model.indx <- which(cov.model == cov.model.names) - 1
  storage.mode(cov.model.indx) <- "integer"
  
  ##Stabilization
  if(is.null(stabilization)){
    if(cov.model == "exponential"){
      stabilization = FALSE
    }
    if(cov.model != "exponential"){
      stabilization = TRUE
    }
  }
  
  if(!isTRUE(stabilization)){
    if(cov.model != "exponential"){
      warning('We recommend using stabilization for all covariance models other than exponential')
    }
  }
  
  if(isTRUE(stabilization)){
    taus <- 10^-6 * sigma.sq
    artificial_noise <- sqrt(taus) * rnorm(n)
    y <- y + artificial_noise
    tau.sq <- tau.sq + taus
    fix_nugget <- 1
    if(verbose == TRUE) cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tUsing Stabilization"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
  }
  
  ##Initial values
  if(cov.model!="matern"){
    initiate <- c(sigma.sq, tau.sq, phi)
    names(initiate) <- c("sigma.sq", "tau.sq", "phi")
  }else{
    initiate <- c(sigma.sq, tau.sq, phi, nu)
    names(initiate) <- c("sigma.sq", "tau.sq", "phi", "nu")}
  
  alpha.sq.starting <- sqrt(tau.sq/sigma.sq)
  phi.starting <- sqrt(phi)
  nu.starting <- sqrt(nu)
  
  storage.mode(alpha.sq.starting) <- "double"
  storage.mode(phi.starting) <- "double"
  storage.mode(nu.starting) <- "double"
  
  ##Search type
  search.type.names <- c("brute", "tree", "cb")
  if(!search.type %in% search.type.names){
    stop("error: specified search.type '",search.type,"' is not a valid option; choose from ", paste(search.type.names, collapse=", ", sep="") ,".")
  }
  search.type.indx <- which(search.type == search.type.names)-1
  storage.mode(search.type.indx) <- "integer"
  
  ##Option for Multithreading if compiled with OpenMp support
  n.omp.threads <- as.integer(n_omp)
  storage.mode(n.omp.threads) <- "integer"
  
  ##type conversion
  storage.mode(y) <- "double"
  storage.mode(x) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"
  storage.mode(coords) <- "double"
  storage.mode(verbose) <- "integer"
  storage.mode(eps) <- "double"
  
  p1<- proc.time()
  
  if(nugget_status == 0){alpha.sq.starting = 0}
  
  n.neighbors = as.integer(max(neighbor$nnIndxLU[(n+1):(2*n)]))
  storage.mode(n.neighbors) <- "integer"
  if(verbose){
    message(paste("Using", n.neighbors, "neighbors"), fill = TRUE)
  }
  nnIndxLU = as.integer(neighbor$nnIndxLU)
  storage.mode(nnIndxLU) <- "integer"
  CIndx = as.integer(neighbor$CIndx)
  storage.mode(CIndx) <- "integer"
  D <- neighbor$D
  storage.mode(D) <- "double"
  d <- neighbor$d
  storage.mode(d) <- "double"
  nnIndx = as.integer(neighbor$nnIndx)
  storage.mode(nnIndx) <- "integer"
  j = as.integer(neighbor$Length.D)
  storage.mode(j) <- "integer"
  nnWt <- as.matrix(neighbor$nnWt)
  storage.mode(nnWt) <- "double"
  Dwt <- as.matrix(neighbor$Dwt)
  storage.mode(Dwt) <- "double"
  
  # cat("Using BRISC_estimateneighborcpp...", fill=TRUE)
  
  result_short <- .Call("BRISC_estimateneighborcpp",
                        y, x,
                        p, n, n.neighbors, coords, 
                        cov.model.indx, alpha.sq.starting, phi.starting, nu.starting, 
                        search.type.indx, n.omp.threads, verbose, eps, fix_nugget,
                        nnIndxLU, CIndx, D, d, nnIndx, j, nnWt, Dwt,
                        PACKAGE = "S3N")
  result <- append(result_short, neighbor)

  p2 <- proc.time()

  Theta <- result$theta
  if(!is.null(pred.stabilization)){
    Theta[2] <- max(result$theta[1] * pred.stabilization, result$theta[2])
  }
  if(cov.model!="matern"){
    names(Theta) <- c("sigma.sq", "tau.sq", "phi")
  }
  else{names(Theta) <- c("sigma.sq", "tau.sq", "phi", "nu")}

  Beta <- result$Beta
  for(i in 1:length(result$Beta)){
    name_beta <- paste0("beta_",i)
    names(Beta)[i] <- name_beta
  }

  llk <- -(result$log_likelihood + n + n * log(2*pi))/2

  result_list <- list ()
  # result_list$ord <- ord
  result_list$coords <- coords
  result_list$y <- y
  result_list$X <- x
  result_list$n.neighbors <- n.neighbors
  result_list$cov.model <- cov.model
  result_list$eps <- eps
  result_list$init <- initiate
  result_list$Beta <- Beta
  result_list$Theta <- Theta
  result_list$log_likelihood <- llk
  result_list$estimation.time <- p2 - p1
  result_list$BRISC_Object <- result
  class(result_list) <- "BRISC_Out"
  result_list
}

# Before 2024 Dec 1, to work with pwdist matrix and neighbor vars from neighbor.R
# BRISC_estimation_stream <- function(coords, pwdist_mat, y, x = NULL,
#                                     sigma.sq = 1, tau.sq = 0.1, phi = 1, nu = 1.5, 
#                                     n.neighbors = 15, n_omp = 1,
#                                     order = "Unif_comp", cov.model = "exponential", 
#                                     search.type = "tree", stabilization = NULL,
#                                     pred.stabilization = 1e-8, verbose = TRUE, 
#                                     eps = 2e-05, nugget_status = 1, ordering = NULL, 
#                                     neighbor = NULL
# ){
#   n <- nrow(coords)
#   if(is.null(x)){
#     x <- matrix(1, nrow = n, ncol = 1)
#   }
#   p <- ncol(x)
#   if(ncol(x) == 1) warning('The ordering of inputs x (covariates) and y (response) in BRISC_estimation has been changed BRISC 1.0.0 onwards.
# Please check the new documentation with ?BRISC_estimation.')
#   ##Coords and ordering
#   if(nugget_status == 0){fix_nugget = 0}
#   if(nugget_status == 1){fix_nugget = 1}
#   # if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
#   # if(order == "AMMD" && length(y) < 65){stop("error: Number of data points must be atleast 65 to use AMMD")}
#   # if(ncol(coords) != 2 || nrow(coords) != n){
#   #   stop("error: either the coords have more than two columns or the number of rows is different than
#   #        data used in the model formula")
#   # }
#   if(is.matrix(y) & ncol(as.matrix(y)) > 1)
#   {stop("error: BRISC handles only univariate response. You have used a y (response) with more than one column.
#           The ordering of inputs x (covariates) and y (response) in BRISC_estimation has been changed BRISC 0.2.0 onwards.
#           Did you mean x (covaiates) instead of y (response) by the multivariatte input?
#           Please check the new documentation with ?BRISC_estimation.")
#     
#   }
#   
#   # coords <- round(coords, tol)
#   # x <- round(x, tol)
#   # y <- round(y, tol)
#   
#   if(tau.sq < 0 ){stop("error: tau.sq must be non-negative")}
#   if(sigma.sq < 0 ){stop("error: sigma.sq must be non-negative")}
#   if(phi < 0 ){stop("error: phi must be non-negative")}
#   if(nu < 0 ){stop("error: nu must be non-negative")}
#   
#   if(is.null(ordering)){
#     if(order != "Unif_comp" && order != "Sum_coords"){
#       stop("error: Please insert a valid ordering scheme choice. Currently the only option is Unif_comp.")
#     }
#     # if(order != "AMMD" && order != "Sum_coords"){
#     #   stop("error: Please insert a valid ordering scheme choice given by 1 or 2.")
#     # }
#     if(verbose == TRUE){
#       cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tOrdering Coordinates"), collapse="   "), "\n")
#     }
#     # if(order == "AMMD"){ord <- orderMaxMinLocal(coords)}
#     if(order == "Sum_coords"){ord <- order(coords[,1] + coords[,2])}
#     if(order == "Unif_comp"){ord <- 1:n} #orderUnifComp(compIDs)}
#   }
#   if(!is.null(ordering)){
#     ord <- ordering
#   }
#   # need drop=FALSE so that if coords has only one column, its class is still
#   # matrix after reordering
#   coords <- coords[ord, ,drop=FALSE]
#   
#   ##Input data and ordering
#   X <- x[ord, ,drop=FALSE]
#   y <- y[ord, ,drop=FALSE]
#   
#   
#   
#   ##Covariance model
#   cov.model.names <- c("exponential","spherical","matern","gaussian")
#   cov.model.indx <- which(cov.model == cov.model.names) - 1
#   storage.mode(cov.model.indx) <- "integer"
#   
#   
#   ##Stabilization
#   if(is.null(stabilization)){
#     if(cov.model == "exponential"){
#       stabilization = FALSE
#     }
#     if(cov.model != "exponential"){
#       stabilization = TRUE
#     }
#   }
#   
#   if(!isTRUE(stabilization)){
#     if(cov.model != "exponential"){
#       warning('We recommend using stabilization for all covariance models other than exponential')
#     }
#   }
#   
#   if(isTRUE(stabilization)){
#     taus <- 10^-6 * sigma.sq
#     artificial_noise <- sqrt(taus) * rnorm(n)
#     y <- y + artificial_noise
#     tau.sq <- tau.sq + taus
#     fix_nugget <- 1
#     if(verbose == TRUE) cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tUsing Stabilization"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
#   }
#   
#   
#   ##Initial values
#   if(cov.model!="matern"){
#     initiate <- c(sigma.sq, tau.sq, phi)
#     names(initiate) <- c("sigma.sq", "tau.sq", "phi")
#   }
#   else{
#     initiate <- c(sigma.sq, tau.sq, phi, nu)
#     names(initiate) <- c("sigma.sq", "tau.sq", "phi", "nu")}
#   
#   alpha.sq.starting <- sqrt(tau.sq/sigma.sq)
#   phi.starting <- sqrt(phi)
#   nu.starting <- sqrt(nu)
#   
#   storage.mode(alpha.sq.starting) <- "double"
#   storage.mode(phi.starting) <- "double"
#   storage.mode(nu.starting) <- "double"
#   
#   
#   ##Search type
#   search.type.names <- c("brute", "tree", "cb")
#   if(!search.type %in% search.type.names){
#     stop("error: specified search.type '",search.type,"' is not a valid option; choose from ", paste(search.type.names, collapse=", ", sep="") ,".")
#   }
#   search.type.indx <- which(search.type == search.type.names)-1
#   storage.mode(search.type.indx) <- "integer"
#   
#   
#   ##Option for Multithreading if compiled with OpenMp support
#   n.omp.threads <- as.integer(n_omp)
#   storage.mode(n.omp.threads) <- "integer"
#   
#   
#   ##type conversion
#   storage.mode(y) <- "double"
#   storage.mode(X) <- "double"
#   storage.mode(p) <- "integer"
#   storage.mode(n) <- "integer"
#   storage.mode(coords) <- "double"
#   storage.mode(n.neighbors) <- "integer"
#   storage.mode(verbose) <- "integer"
#   storage.mode(eps) <- "double"
#   
#   
#   
#   p1<- proc.time()
#   
#   if(nugget_status == 0){alpha.sq.starting = 0}
#   
#   
#   ##note: this isn't going to work for streams
#   if(is.null(neighbor)){
#     cat("Getting neighbor variables...", fill=TRUE)
#     neighbors = BRISC_neighbor_stream(
#       coords = as.matrix(1:length(obs$COMID)),
#       pwdist_matrix = pwdist_matrix,
#       n.neighbors = n_neighbors,
#       order = "Unif_comp", 
#       n_omp = 1, search.type = "tree",
#       verbose = TRUE, ordering = NULL, tol = 12
#     )
#   }
#   
#   if(!is.null(neighbor)){
#     nnIndxLU = neighbor$nnIndxLU
#     # nnIndxLU = as.integer(neighbor$nnIndxLU)
#     storage.mode(nnIndxLU) <- "integer"
#     CIndx = neighbor$CIndx
#     # CIndx = as.integer(neighbor$CIndx)
#     storage.mode(CIndx) <- "integer"
#     D <- neighbor$D
#     storage.mode(D) <- "double"
#     d <- neighbor$d
#     storage.mode(d) <- "double"
#     nnIndx = neighbor$nnIndx
#     # nnIndx = as.integer(neighbor$nnIndx)
#     storage.mode(nnIndx) <- "integer"
#     j = neighbor$Length.D
#     # j = as.integer(neighbor$j)
#     storage.mode(j) <- "integer"
#     
#     # convert pwdist_matrix into vector (by column)
#     pwdist_mat = as.vector(pwdist_mat)
#     storage.mode(pwdist_mat) <- "double"
#     
#     print("Using BRISC_estimateneighborcpp_s...")
#     
#     result_short <- .Call("BRISC_estimateneighborcpp_s",
#                           y, X, pwdist_mat,
#                           p, n, n.neighbors, coords, 
#                           cov.model.indx, alpha.sq.starting, phi.starting, nu.starting, 
#                           search.type.indx, n.omp.threads, verbose, eps, fix_nugget,
#                           nnIndxLU, CIndx, D, d, nnIndx, j,
#                           PACKAGE = "BRISC")
#     result <- append(result_short, neighbor)
#     
#   }
#   
#   # result_short
#   
#   
#   p2 <- proc.time()
# 
#   Theta <- result$theta
#   if(!is.null(pred.stabilization)){
#     Theta[2] <- max(result$theta[1] * pred.stabilization, result$theta[2])
#   }
#   if(cov.model!="matern"){
#     names(Theta) <- c("sigma.sq", "tau.sq", "phi")
#   }
#   else{names(Theta) <- c("sigma.sq", "tau.sq", "phi", "nu")}
# 
#   Beta <- result$Beta
#   for(i in 1:length(result$Beta)){
#     name_beta <- paste0("beta_",i)
#     names(Beta)[i] <- name_beta
#   }
# 
#   llk <- -(result$log_likelihood + n + n * log(2*pi))/2
#   
#   result_list <- list ()
#   result_list$ord <- ord
#   result_list$coords <- coords
#   result_list$y <- y
#   result_list$X <- X
#   result_list$n.neighbors <- n.neighbors
#   result_list$cov.model <- cov.model
#   result_list$eps <- eps
#   result_list$init <- initiate
#   result_list$Beta <- Beta
#   result_list$Theta <- Theta
#   result_list$log_likelihood <- llk
#   result_list$estimation.time <- p2 - p1
#   result_list$BRISC_Object <- result
#   class(result_list) <- "BRISC_Out"
#   result_list
# }




BRISC_estimation <- function(coords, y, x = NULL, sigma.sq = 1, tau.sq = 0.1, phi = 1, nu = 1.5, n.neighbors = 15, n_omp = 1,
                             order = "Sum_coords", cov.model = "exponential", search.type = "tree", stabilization = NULL,
                             pred.stabilization = 1e-8, verbose = TRUE, eps = 2e-05, nugget_status = 1, ordering = NULL, neighbor = NULL, tol = 12
){
  n <- nrow(coords)
  if(is.null(x)){
    x <- matrix(1, nrow = n, ncol = 1)
  }
  p <- ncol(x)
  if(ncol(x) == 1) warning('The ordering of inputs x (covariates) and y (response) in BRISC_estimation has been changed BRISC 1.0.0 onwards.
Please check the new documentation with ?BRISC_estimation.')
  ##Coords and ordering
  if(nugget_status == 0){fix_nugget = 0}
  if(nugget_status == 1){fix_nugget = 1}
  if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
  if(order == "AMMD" && length(y) < 65){stop("error: Number of data points must be atleast 65 to use AMMD")}
  if(ncol(coords) != 2 || nrow(coords) != n){
    stop("error: either the coords have more than two columns or the number of rows is different than
         data used in the model formula")
  }
  if(is.matrix(y) & ncol(as.matrix(y)) > 1)
    {stop("error: BRISC handles only univariate response. You have used a y (response) with more than one column.
          The ordering of inputs x (covariates) and y (response) in BRISC_estimation has been changed BRISC 0.2.0 onwards.
          Did you mean x (covaiates) instead of y (response) by the multivariatte input?
          Please check the new documentation with ?BRISC_estimation.")

  }

  coords <- round(coords, tol)
  x <- round(x, tol)
  y <- round(y, tol)

  if(tau.sq < 0 ){stop("error: tau.sq must be non-negative")}
  if(sigma.sq < 0 ){stop("error: sigma.sq must be non-negative")}
  if(phi < 0 ){stop("error: phi must be non-negative")}
  if(nu < 0 ){stop("error: nu must be non-negative")}

  if(is.null(ordering)){
    if(order != "AMMD" && order != "Sum_coords"){
      stop("error: Please insert a valid ordering scheme choice given by 1 or 2.")
    }
    if(verbose == TRUE){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tOrdering Coordinates"), collapse="   "), "\n")
    }
    if(order == "AMMD"){ord <- orderMaxMinLocal(coords)}
    if(order == "Sum_coords"){ord <- order(coords[,1] + coords[,2])}
  }
  if(!is.null(ordering)){
    ord <- ordering
  }
  coords <- coords[ord,]

  ##Input data and ordering
  X <- x[ord,,drop=FALSE]
  y <- y[ord]



  ##Covariance model
  cov.model.names <- c("exponential","spherical","matern","gaussian")
  cov.model.indx <- which(cov.model == cov.model.names) - 1
  storage.mode(cov.model.indx) <- "integer"


  ##Stabilization
  if(is.null(stabilization)){
    if(cov.model == "exponential"){
      stabilization = FALSE
    }
    if(cov.model != "exponential"){
      stabilization = TRUE
    }
  }

  if(!isTRUE(stabilization)){
    if(cov.model != "exponential"){
      warning('We recommend using stabilization for spherical, Matern and Gaussian covariance model')
    }
  }

  if(isTRUE(stabilization)){
    taus <- 10^-6 * sigma.sq
    artificial_noise <- sqrt(taus) * rnorm(n)
    y <- y + artificial_noise
    tau.sq <- tau.sq + taus
    fix_nugget <- 1
    if(verbose == TRUE) cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tUsing Stabilization"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
  }


  ##Initial values
  if(cov.model!="matern"){
    initiate <- c(sigma.sq, tau.sq, phi)
    names(initiate) <- c("sigma.sq", "tau.sq", "phi")
  }
  else{
    initiate <- c(sigma.sq, tau.sq, phi, nu)
    names(initiate) <- c("sigma.sq", "tau.sq", "phi", "nu")}

  alpha.sq.starting <- sqrt(tau.sq/sigma.sq)
  phi.starting <- sqrt(phi)
  nu.starting <- sqrt(nu)

  storage.mode(alpha.sq.starting) <- "double"
  storage.mode(phi.starting) <- "double"
  storage.mode(nu.starting) <- "double"


  ##Search type
  search.type.names <- c("brute", "tree", "cb")
  if(!search.type %in% search.type.names){
    stop("error: specified search.type '",search.type,"' is not a valid option; choose from ", paste(search.type.names, collapse=", ", sep="") ,".")
  }
  search.type.indx <- which(search.type == search.type.names)-1
  storage.mode(search.type.indx) <- "integer"


  ##Option for Multithreading if compiled with OpenMp support
  n.omp.threads <- as.integer(n_omp)
  storage.mode(n.omp.threads) <- "integer"


  ##type conversion
  storage.mode(y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"
  storage.mode(coords) <- "double"
  storage.mode(n.neighbors) <- "integer"
  storage.mode(verbose) <- "integer"
  storage.mode(eps) <- "double"



  p1<- proc.time()

  if(nugget_status == 0){alpha.sq.starting = 0}


  ##estimtion
  if(is.null(neighbor)){
    print("Using BRISC_estimatecpp...")
    result <- .Call("BRISC_estimatecpp", y, X, p, n, n.neighbors, coords, cov.model.indx, alpha.sq.starting, phi.starting, nu.starting, search.type.indx, n.omp.threads, verbose, eps, fix_nugget, PACKAGE = "S3N")
  }

  if(!is.null(neighbor)){
    nnIndxLU = neighbor$nnIndxLU
    # nnIndxLU = as.integer(neighbor$nnIndxLU)
    storage.mode(nnIndxLU) <- "integer"
    CIndx = neighbor$CIndx
    # CIndx = as.integer(neighbor$CIndx)
    storage.mode(CIndx) <- "integer"
    D <- neighbor$D
    storage.mode(D) <- "double"
    d <- neighbor$d
    storage.mode(d) <- "double"
    nnIndx = neighbor$nnIndx
    # nnIndx = as.integer(neighbor$nnIndx)
    storage.mode(nnIndx) <- "integer"
    j = neighbor$j
    # j = as.integer(neighbor$j)
    storage.mode(j) <- "integer"
    
    print("Using BRISC_estimateneighborcpp...")

    result_short <- .Call("BRISC_estimateneighborcpp", y, X, p, n, n.neighbors, coords, cov.model.indx, alpha.sq.starting, phi.starting, nu.starting, search.type.indx, n.omp.threads, verbose, eps, fix_nugget,
                    nnIndxLU, CIndx, D, d, nnIndx, j, PACKAGE = "S3N")
    result <- append(result_short, neighbor)

  }


  p2 <- proc.time()

  Theta <- result$theta
  if(!is.null(pred.stabilization)){
    Theta[2] <- max(result$theta[1] * pred.stabilization, result$theta[2])
  }
  if(cov.model!="matern"){
    names(Theta) <- c("sigma.sq", "tau.sq", "phi")
  }
  else{names(Theta) <- c("sigma.sq", "tau.sq", "phi", "nu")}

  Beta <- result$Beta
  for(i in 1:length(result$Beta)){
    name_beta <- paste0("beta_",i)
    names(Beta)[i] <- name_beta
  }

  llk <- -(result$log_likelihood + n + n * log(2*pi))/2

  result_list <- list ()
  result_list$ord <- ord
  result_list$coords <- coords
  result_list$y <- y
  result_list$X <- X
  result_list$n.neighbors <- n.neighbors
  result_list$cov.model <- cov.model
  result_list$eps <- eps
  result_list$init <- initiate
  result_list$Beta <- Beta
  result_list$Theta <- Theta
  result_list$log_likelihood <- llk
  result_list$estimation.time <- p2 - p1
  result_list$BRISC_Object <- result
  class(result_list) <- "BRISC_Out"
  result_list
}
