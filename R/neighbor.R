# Neighbor variables:
# nnIndxLU: length 2n
#   - essentially two rows concatenated into one:
#       1) index of nnIndx that corresponds to the first neighbor of point i
#           e.g. nnIndx[2] equals the index of the first neighbor of the third point (zero-indexed)
#       2) number of neighbors for each point
#           e.g. nnIndx[n+2] equals the number of neighbors the third point has (zero-indexed)
#
# CIndx: length 2n
#   - same as nnIndxLU except for covariance and distance matrices:
#       1) index of C (cov matrix vector) that corresponds to the first element of the cov matrix for point i
#       2) number of cov matrix elements for each point
#   - number of elements in the covariance/distance matrices of the neighbors 
#     for point i = (# neighbors of point i)^2 = (nnIndxLU[n+i])^2
#
# D: length = sum over i of (# neighbors of point i)^2; this equals CIndx[n]+CIndx[2n]
#
# d: length = total # neighbors = m*(1+m)/2 + m*(n-m-1)
#
# nnIndx: length = total # neighbors
#   - vector of indices of the neighbor points
#
# Length.D = length of D


stream_dist = function(i,j){
  # if(pwdist_matrix[i,j] == 0){
  #   print(paste0("For i=", i, " j=", j, ": d = 0"))
  #   break()
  # }
  return(pwdist_matrix[i,j])
}

#' orderUnifComp
#' 
#' @param compIDs connected component IDs
#' 
#' @return The ordering
#' @importFrom dplyr n_distinct
#' @export
orderUnifComp = function(compIDs){
  n = length(compIDs)
  inds = 1:n
  if(n_distinct(compIDs) == 1){
    ord = sample(inds, n)
    # draw a permutation of the connected components uniformly at random
    comps_perm = sample(unique(compIDs), n_distinct(compIDs))
    ord = unlist(lapply(comps_perm, function(x) sample(inds[compIDs==x], sum(compIDs==x))))
  }
  return(ord)
}



# for stream distance, coords is still two columns but they're COMID and connected component ID, not X and Y
# also the only ordering scheme we allow currently is uniform by connected component (Unif_comp)

BRISC_neighbor_stream <- function(coords, pwdist_matrix, n.neighbors = 15, n_omp = 1,
                                  order = "Unif_comp", search.type = "tree", 
                                  verbose = TRUE, ordering = NULL, tol = 12
){
  print("Without pwdist_matrix argument")
  n <- nrow(coords)
  print(paste0("ordering method: ", order))
  
  # do I want to have an analog of this part of BRISC_neighbor in BRISC_neighbor_stream?
  #   coords <- round(coords, tol)
  
  # if(!is.matrix(coords)){stop("error: coords must be n-by-2 matrix of xy-coordinate locations")}
  # if(order == "AMMD" && nrow(coords) < 65){stop("error: Number of data points must be atleast 65 to use AMMD")}
  # if(ncol(coords) != 2 || nrow(coords) != n){
  #   stop("error: either the coords have more than two columns or the number of rows is different than
  #        data used in the model formula")
  # }

  # coords <- round(coords, tol)


  if(is.null(ordering)){
    if(order != "Unif_comp"){
      stop("error: Please insert a valid ordering scheme choice. Currently the only option is Unif_comp.")
    }
    # if(order != "AMMD" && order != "Sum_coords"){
    #   stop("error: Please insert a valid ordering scheme choice given by 1 or 2.")
    # }
    if(verbose == TRUE){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tOrdering Coordinates"), collapse="   "), "\n")
    }
    # if(order == "AMMD"){ord <- orderMaxMinLocal(coords)}
    # if(order == "Sum_coords"){ord <- order(coords[,1] + coords[,2])}
    if(order == "Unif_comp"){ord <- 1:n} #orderUnifComp(compIDs)}
  }
  if(!is.null(ordering)){
    ord <- ordering
  }
  # need drop=FALSE so that if coords has only one column, its class is still
  # matrix after reordering
  coords <- coords[ord, ,drop=FALSE]

  ##Search type
  search.type.names <- c("brute", "tree", "cb")
  if(!search.type %in% search.type.names){
    stop("error: specified search.type '",search.type,"' is not a valid option; choose from ", paste(search.type.names, collapse=", ", sep="") ,".")
  }
  search.type.indx <- which(search.type == search.type.names)-1
  # storage.mode(search.type.indx) <- "integer"

  ##Option for Multithreading if compiled with OpenMp support
  n.omp.threads <- as.integer(n_omp)
  # storage.mode(n.omp.threads) <- "integer"
  # 
  # 
  # ##type conversion
  # storage.mode(n) <- "integer"
  # storage.mode(coords) <- "double"
  # storage.mode(n.neighbors) <- "integer"
  # storage.mode(verbose) <- "integer"



  p1 <- proc.time()
  
  # computing the six neighbor objects in R for streams ----------------------
  
  # this first part is analogous to getNNIndx and mkNNIndx -------------------
  # note: R is 1-indexed, C is 0-indexed
  # this code needs to be 1-indexed so that these six objects are created correctly in R
  # but the BRISC C code that uses them later is 0-indexed (doesn't affect this R code)
  m = n.neighbors
  
  # first we build nIndx, d, nnIndxLU, and nnIndx
  nIndx = m*(1+m)/2 + m*(n-m-1) # total number of neighbors
  d = rep(Inf, nIndx)
  nnIndxLU = rep(Inf, 2*n)
  nnIndx = rep(Inf, nIndx)
  
  for(i in 0:(n-1)){
    Ri = i+1
    # getNNIndx(i, m, iNNIndx, iNN)
    # here, iNNIndx is the index of the first neighbor and
    # iNN is the number of neighbors for this point i (indexed Ri=i+1 in R)
    if(i == 0){
      iNNIndx = 0
      iNN = 0
    }else if(i < m){
      iNNIndx = i/2*(i-1)
      iNN = i
    }else{
      iNNIndx = m/2*(m-1)+(i-m)*m
      iNN = m
    }
    
    nnIndxLU[Ri] = iNNIndx
    nnIndxLU[n+Ri] = iNN
    
    if(i != 0){
      for(j in 0:(i-1)){
        # print(paste0("i=", i, " j=", j))
        Rj = j+1
        # print(paste0("Ri=", Ri, " Rj=", Rj))
        thisd = stream_dist(coords[Ri,1], coords[Rj,1])
        d_ind = iNNIndx+iNN
        # print(paste0("d_ind: ", d_ind))
        if(thisd < d[d_ind]){
          d[d_ind] = thisd
          nnIndx[d_ind] = j
          # sort in increasing order just within the neighbors for point Ri
          ind_to_sort = (iNNIndx+1):d_ind
          # print(ind_to_sort)
          # d[ind_to_sort] = sort(d[ind_to_sort])
          # nnIndx[ind_to_sort] = sort(nnIndx[ind_to_sort])
          nnIndx[ind_to_sort] = nnIndx[ind_to_sort][order(d[ind_to_sort])]
          d[ind_to_sort] = d[ind_to_sort][order(d[ind_to_sort])]
        }
      }
    }
    # iNNIndx; iNN; nnIndx; d; nnIndxLU
  }
  
  # [1] 5
  # [1] 2
  # [1] 0 1 0 1 2 1 2
  # [1] 2.150628 3.026919 5.177546 4.469291 4.623351 3.876049 4.030109
  # [1] 0 0 1 3 5 0 1 2 2 2
  # these dists corresp. to (2,1), (3,2), (3,1), (4,2), (4,3), (5,2), (5,3)
  # now let's verify these results:
  # pwdists_5 = pwdists %>%
  #   dplyr::filter(COMID1 %in% coords[,1], COMID2 %in% coords[,1])
  # pwdists_5$C1 = c(1,1,1,1,2,2,2,3,3,4)
  # pwdists_5$C2 = c(2,3,4,5,3,4,5,4,5,5)
  # using pwdists_5 to look up the pairwise distances between the points above,
  # i.e. (2,1), (3,2), (3,1), (4,2), (4,3), (5,2), (5,3), gives us the same d:
  # 2.150628 3.026919 5.177546 4.469291 4.623351 3.876049 4.030109
  
  # this next part is analogous to mkNNIndx -------------------
  
  # now let's make CIndx
  CIndx = rep(NA, 2*n)
  
  for(i in 1:n){
    if(i == 1){
      CIndx[n+i] = 0
      CIndx[i] = 0
    }else{
      CIndx[n+i] = nnIndxLU[n+i]*nnIndxLU[n+i]
      CIndx[i] = CIndx[n+i-1] + CIndx[i-1]
    }
  }
  # for five points and two neighbors, looks like this:
  # [1] 0 0 1 5 9 0 1 4 4 4
  
  # now let's make D
  D = rep(0, CIndx[n]+CIndx[2*n])
  
  # ind = 1
  for(i in 1:(n-1)){ # for each point i...
    # print(paste0("loop ", ind))
    Ri = i+1 # R is 1-indexed
    # print(paste0("i=", i))
    for(k in 0:(nnIndxLU[n+Ri]-1)){ # for each of its neighbors, of which there are nnIndxLU[n+i] many...
      # print(paste0("k=", k))
      for(l in 0:k){ # we'll make a k by k matrix as a vector
        # print(paste0("i=", i, ", k=", k, ", l=", l))
        # Ri = i+1; Rk = k+1; Rl = l+1
        # print(paste0("CIndx[Ri]=", CIndx[Ri]))
        # print(paste0("nnIndxLU[Ri] = ", nnIndxLU[Ri]))
        # print(paste0("nnIndxLU[n+Ri]=", nnIndxLU[n+Ri]))
        # print(paste0("nnIndx[nnIndxLU[Ri]+k+1]+1 = ", nnIndx[nnIndxLU[Ri]+k+1]+1))
        # print(paste0("nnIndx[nnIndxLU[Ri]+l+1]+1 = ", nnIndx[nnIndxLU[Ri]+l+1]+1))
        # print(paste0("coords[nnIndx[nnIndxLU[Ri]+k+1]+1, 1] = ", coords[nnIndx[nnIndxLU[Ri]+k+1]+1, 1]))
        # print(paste0("coords[nnIndx[nnIndxLU[Ri]+l+1]+1, 1] = ", coords[nnIndx[nnIndxLU[Ri]+l+1]+1, 1]))
        # print(paste0("k=", k))
        # print(paste0("res = ", CIndx[Ri]+l*nnIndxLU[n+Ri]+k))
        ind_D = CIndx[i]+l*nnIndxLU[n+i]+k
        D[ind_D+1] = stream_dist(coords[nnIndx[nnIndxLU[Ri]+k+1]+1, 1],
                                 coords[nnIndx[nnIndxLU[Ri]+l+1]+1, 1])
      }
    }
    # ind = ind+1
  }
  
  Length.D = length(D)
  
  nnIndxLU = as.integer(nnIndxLU)
  CIndx = as.integer(CIndx)
  nnIndx = as.integer(nnIndx)
  
  result = list(nnIndxLU = nnIndxLU, CIndx = CIndx, D = D, d = d, nnIndx = nnIndx, Length.D = Length.D)

  p2 <- proc.time()

  result
}





BRISC_neighbor <- function(coords, n.neighbors = 15, n_omp = 1,
                           order = "Sum_coords", search.type = "tree", verbose = TRUE, ordering = NULL, tol = 12
){
  n <- nrow(coords)
  if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
  if(order == "AMMD" && nrow(coords) < 65){stop("error: Number of data points must be atleast 65 to use AMMD")}
  if(ncol(coords) != 2 || nrow(coords) != n){
    stop("error: either the coords have more than two columns or the number of rows is different than
         data used in the model formula")
  }
  
  coords <- round(coords, tol)
  
  
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
  storage.mode(n) <- "integer"
  storage.mode(coords) <- "double"
  storage.mode(n.neighbors) <- "integer"
  storage.mode(verbose) <- "integer"
  
  
  
  p1<- proc.time()
  
  ##estimtion
  result <- .Call("BRISC_neighborcpp", n, n.neighbors, coords, search.type.indx, n.omp.threads, verbose, PACKAGE = "S3N")
  
  p2 <- proc.time()
  
  result
}
