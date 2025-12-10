bootstrap_brisc <- function(i, # 1:n_boot (seed)
                            e, # norm.residual
                            X_data, B, efF, Xbeta, D, d, 
                            nnIndx, nnIndxLU, CIndx, nnWt, Dwt,
                            n, p, n.neighbors, 
                            theta_not, cov.model.indx, 
                            Length.D, n.omp.threads, eps, fix_nugget){
  # print("Starting bootstrap_brisc\n")
  set.seed(i)
  boot_indices = sample(1:n, n, replace = TRUE)
  e_boot <- e[boot_indices]
  # print("Fine till BRISC_bootstrapcpp") #<= issue is in BRISC_bootstrapcpp
  boot_statistics <- .Call(
    "BRISC_bootstrapcpp", X_data, B, efF, Xbeta, e_boot,
    D, d, nnIndx, nnIndxLU, CIndx, nnWt, Dwt,
    n, p, n.neighbors, theta_not, cov.model.indx, Length.D, 
    n.omp.threads, eps, fix_nugget, PACKAGE = "S3N")
  # print("Fine after BRISC_bootstrapcpp")
  if(cov.model.indx == 2){
    result = c(
      Beta = boot_statistics$Beta, 
      sigma.sq = boot_statistics$theta[1], 
      tau.sq = boot_statistics$theta[2], 
      phi = boot_statistics$theta[3], 
      nu = boot_statistics$theta[4])
  }
  if(cov.model.indx != 2){
    result = c(
      Beta = boot_statistics$Beta, 
      sigma.sq = boot_statistics$theta[1], 
      tau.sq = boot_statistics$theta[2], 
      phi = boot_statistics$theta[3])
  }
  result
}
