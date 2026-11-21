#define USE_FC_LEN_T
#include <string>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rconfig.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#ifndef FCONE
#define FCONE
#endif

extern "C" {

//   SEXP BRISC_predictioncpp(SEXP X_r, SEXP y_r, SEXP afv_r, SEXP coords_r, SEXP n_r, SEXP p_r, SEXP m_r, SEXP X0_r, SEXP coords0_r, SEXP q_r, SEXP nnIndx0_r,
// 		    SEXP betaSamples_r, SEXP thetaSamples_r, SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r){
// 
//     int i, k, l, info, nProtect=0;
//     const int inc = 1;
//     const double one = 1.0;
//     const double zero = 0.0;
//     char const *lower = "L";
// 
// 
//     //get args
//     double *X = REAL(X_r);
//     double *y = REAL(y_r);
//     double *afv = REAL(afv_r);
//     double *coords = REAL(coords_r);
//     int n = INTEGER(n_r)[0];
//     int p = INTEGER(p_r)[0];
//     int m = INTEGER(m_r)[0];
//     int mm = m*m;
// 
//     double *X0 = REAL(X0_r);
//     double *coords0 = REAL(coords0_r);
//     int q = INTEGER(q_r)[0];
// 
//     int *nnIndx0 = INTEGER(nnIndx0_r);
//     double *beta = REAL(betaSamples_r);
//     double *theta = REAL(thetaSamples_r);
//     int covModel = INTEGER(covModel_r)[0];
//     std::string corName = getCorName(covModel);
//     int nThreads = INTEGER(nThreads_r)[0];
//     int verbose = INTEGER(verbose_r)[0];
// 
// #ifdef _OPENMP
//     omp_set_num_threads(nThreads);
// #else
//     if(nThreads > 1){
//       warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
//       nThreads = 1;
//     }
// #endif
// 
//     if(verbose){
//       Rprintf("----------------------------------------\n");
//       Rprintf("\tPrediction description\n");
//       Rprintf("----------------------------------------\n");
//       Rprintf("BRISC model fit with %i observations.\n\n", n);
//       Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
//       Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
//       Rprintf("Using %i nearest neighbors.\n\n", m);
//       Rprintf("Predicting at %i locations.\n\n", q);
// #ifdef _OPENMP
//       Rprintf("\nSource compiled with OpenMP support and model fit using %i threads.\n", nThreads);
// #else
//       Rprintf("\n\nSource not compiled with OpenMP support.\n");
// #endif
//     }
// 
//     //parameters
//     int nTheta, sigmaSqIndx, tauSqIndx, phiIndx, nuIndx;
// 
//     if(corName != "matern"){
//       nTheta = 3;//sigma^2, tau^2, phi
//       sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2;
//     }else{
//       nTheta = 4;//sigma^2, tau^2, phi, nu
//       sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2; nuIndx = 3;
//     }
// 
//     //get max nu
//     double nuMax = 0;
//     int nb = 0;
// 
//     if(corName == "matern"){
// 	if(theta[nuIndx] > nb){
// 	  nb = theta[nuIndx];
// 	}
// 
//       nb = 1+static_cast<int>(floor(nuMax));
//     }
// 
//     double *bk = (double *) R_alloc(nThreads*nb, sizeof(double));
// 
//     double *C = (double *) R_alloc(nThreads*mm, sizeof(double)); zeros(C, nThreads*mm);
//     double *c = (double *) R_alloc(nThreads*m, sizeof(double)); zeros(c, nThreads*m);
//     double *tmp_m  = (double *) R_alloc(nThreads*m, sizeof(double));
//     double phi = 0, nu = 0, sigmaSq = 0, tauSq = 0, d;
//     int threadID = 0;
// 
//     SEXP y0_r;
//     SEXP vary0_r;
//     PROTECT(vary0_r = allocMatrix(REALSXP, q, 1)); nProtect++;
//     double *vary0 = REAL(vary0_r);
//     PROTECT(y0_r = allocMatrix(REALSXP, q, 1)); nProtect++;
//     double *y0 = REAL(y0_r);
// 
//     if(verbose){
//       Rprintf("-------------------------------------------------\n");
//       Rprintf("\t\tPredicting\n");
//       Rprintf("-------------------------------------------------\n");
//       #ifdef Win32
//         R_FlushConsole();
//       #endif
//     }
// 
//     double *z = (double *) R_alloc(q, sizeof(double));
//     GetRNGstate();
//     for(i = 0; i < q; i++){
//       z[i] = rnorm(0.0,1.0);
//     }
//     PutRNGstate();
// 
// 
// 
// 
// #ifdef _OPENMP
// #pragma omp parallel for private(threadID, phi, nu, sigmaSq, tauSq, k, l, d, info)
// #endif
//     for(i = 0; i < q; i++){
// #ifdef _OPENMP
// 	threadID = omp_get_thread_num();
// #endif
// 	phi = theta[phiIndx];
// 	if(corName == "matern"){
// 	  nu = theta[nuIndx];
// 	}
// 	sigmaSq = theta[sigmaSqIndx];
// 	tauSq = theta[tauSqIndx];
// 
// 	for(k = 0; k < m; k++){
// 	  d = dist2(coords[nnIndx0[i+q*k]], coords[n+nnIndx0[i+q*k]], coords0[i], coords0[q+i]);
// 	  c[threadID*m+k] = sigmaSq*spCor(d, phi, nu, afv[i], afv[nnIndxLU[i]+k], covModel, &bk[threadID*nb]);
// 	  for(l = 0; l < m; l++){
// 	    d = dist2(coords[nnIndx0[i+q*k]], coords[n+nnIndx0[i+q*k]], coords[nnIndx0[i+q*l]], coords[n+nnIndx0[i+q*l]]);
// 	    C[threadID*mm+l*m+k] = sigmaSq*spCor(d, phi, nu, afv[nnIndxLU[i]+l], afv[nnIndxLU[i]+k], covModel, &bk[threadID*nb]);
// 	    if(k == l){
// 	      C[threadID*mm+l*m+k] += tauSq;
// 	    }
// 	  }
// 	}
// 
// 	F77_NAME(dpotrf)(lower, &m, &C[threadID*mm], &m, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
// 	F77_NAME(dpotri)(lower, &m, &C[threadID*mm], &m, &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
// 
// 	F77_NAME(dsymv)(lower, &m, &one, &C[threadID*mm], &m, &c[threadID*m], &inc, &zero, &tmp_m[threadID*m], &inc FCONE);
// 
// 	d = 0;
// 	for(k = 0; k < m; k++){
// 	  d += tmp_m[threadID*m+k]*(y[nnIndx0[i+q*k]] - F77_NAME(ddot)(&p, &X[nnIndx0[i+q*k]], &n, &beta[0], &inc));
// 	}
// 
// 	y0[i] = F77_NAME(ddot)(&p, &X0[i], &q, &beta[0], &inc) + d;
//     vary0[i] = std::max(0.0, sqrt(sigmaSq + tauSq - F77_NAME(ddot)(&m, &tmp_m[threadID*m], &inc, &c[threadID*m], &inc)));
//     }
// 
// 
//     //make return object
//     SEXP result_r, resultName_r;
//     int nResultListObjs = 1 + 1;
// 
//     PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
//     PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
// 
//     SET_VECTOR_ELT(result_r, 0, y0_r);
//     SET_VECTOR_ELT(resultName_r, 0, mkChar("p.y.0"));
// 
//     SET_VECTOR_ELT(result_r, 1, vary0_r);
//     SET_VECTOR_ELT(resultName_r, 1, mkChar("var.y.0"));
// 
//     namesgets(result_r, resultName_r);
// 
//     //unprotect
//     UNPROTECT(nProtect);
// 
//     return(result_r);
// 
//   }

    
SEXP BRISC_predictioncpp_s(SEXP X_r, SEXP y_r, SEXP coords_r, SEXP n_r, SEXP p_r, 
                           SEXP m_r, SEXP X0_r, SEXP coords0_r, SEXP q_r, 
                           SEXP nnIndx_r, SEXP nnDist_r, SEXP nnWght_r,
                           SEXP obsobs_d_r, SEXP obsobs_wt_r,
                           SEXP betaSamples_r, SEXP thetaSamples_r,
                           SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r){

  int i, k, l, info, nProtect=0;
  const int inc = 1;
  const double one = 1.0;
  const double zero = 0.0;
  char const *lower = "L";

  //get args
  double *X = REAL(X_r);
  double *y = REAL(y_r);
  double *coords = REAL(coords_r);
  int n = INTEGER(n_r)[0];
  int p = INTEGER(p_r)[0];
  int m = INTEGER(m_r)[0];
  int mm = m*m;

  double *X0 = REAL(X0_r);
  double *coords0 = REAL(coords0_r);
  int q = INTEGER(q_r)[0];
  // Rprintf("q: %i", q);
  
  // pred-obs nearest neighbor indices and distances
  int *nnIndx0 = INTEGER(nnIndx_r);
  double *nnDist0 = REAL(nnDist_r);
  double *nnWght0 = REAL(nnWght_r);
  // obs-obs nearest neighbor distances (pwdist_matrix)
  double *obsobs_d = REAL(obsobs_d_r);
  double *obsobs_wt = REAL(obsobs_wt_r);
  
  // the first four values go down the first column
  // obsobs_d[n-1] is the bottom of the first column
  // obsobs_d[n] is the top of the second column
  // obsobs_d[3*n + 5] is the distance between the 3rd and 5th points
  // obsobs_d[5*n + 3] returns zero because it's above the diagonal
  // Rprintf("obsobs_d: %f %f %f %f %f %f %f %f %f", obsobs_d[0], obsobs_d[1], obsobs_d[2], obsobs_d[3], obsobs_d[n-1], obsobs_d[n], obsobs_d[n+1], obsobs_d[5*n + 3], obsobs_d[3*n + 5]);
  
  double *beta = REAL(betaSamples_r);
  double *theta = REAL(thetaSamples_r);
  int covModel = INTEGER(covModel_r)[0];
  std::string corName = getCorName(covModel);
  int nThreads = INTEGER(nThreads_r)[0];
  int verbose = INTEGER(verbose_r)[0];
  
  // Rprintf("nnIndx0: %i %i %i %i %i", nnIndx0[0], nnIndx0[1], nnIndx0[2], nnIndx0[3], nnIndx0[q]);
  // Rprintf("nnDist0: %f %f %f %f %f %f %f %f %f", nnDist0[0], nnDist0[1], nnDist0[2], nnDist0[3], nnDist0[8906], nnDist0[8907], nnDist0[q-1], nnDist0[q], nnDist0[2*q]);
  // should match this in R: nn.dists[c(0:3, 8906, 8907, q-1, q, 2*q) + 1]
  // return(0);
  // so the indexing goes down columns first, not across rows
  // nn.indx[1:4,1:4]
  // [,1] [,2] [,3] [,4]
  // [1,] 3114 3126 3120 5403
  // [2,] 3114 3126 3120 5403
  // [3,] 3114 3126 3120 5403
  // [4,] 3114 3126 3120 3116
  
  // Rprintf("nnDist[3113, 3113]: %f \n", nnDist[3113, 3113]);

  if(nThreads < 1){
    nThreads = 1;
  }
#ifndef _OPENMP
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif

  if(verbose){
    Rprintf("----------------------------------------\n");
    Rprintf("\tPrediction description\n");
    Rprintf("----------------------------------------\n");
    Rprintf("BRISC model fit with %i observations.\n\n", n);
    Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
    Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
    Rprintf("Using %i nearest neighbors.\n\n", m);
    Rprintf("Predicting at %i locations.\n\n", q);
#ifdef _OPENMP
    Rprintf("\nSource compiled with OpenMP support and model fit using %i threads.\n", nThreads);
#else
    Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
  }

  //parameters
  int nTheta, sigmaSqIndx, tauSqIndx, phiIndx, nuIndx;

  if(corName != "matern"){
    nTheta = 3;//sigma^2, tau^2, phi
    sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2;
  }else{
    nTheta = 4;//sigma^2, tau^2, phi, nu
    sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2; nuIndx = 3;
  }

  //get max nu
  double nuMax = 0;
  int nb = 0;

  // i = 3, k = 2, that's the 2nd row, 3rd col, so pred 2, obs 3
  // printf("%i", nnIndx0[3 + q*2]);

  // return(0);

  if(corName == "matern"){
    if(theta[nuIndx] > nuMax){
      nuMax = theta[nuIndx];
    }

    nb = 1+static_cast<int>(floor(nuMax));
  }

  double *bk = (double *) R_alloc(nThreads*nb, sizeof(double));

  double *C = (double *) R_alloc(nThreads*mm, sizeof(double)); zeros(C, nThreads*mm);
  double *c = (double *) R_alloc(nThreads*m, sizeof(double)); zeros(c, nThreads*m);
  double *tmp_m  = (double *) R_alloc(nThreads*m, sizeof(double));
  double phi = 0, nu = 0, sigmaSq = 0, tauSq = 0, d;
  int threadID = 0;
  int err_flag = 0;

  SEXP y0_r;
  SEXP vary0_r;
  PROTECT(vary0_r = allocMatrix(REALSXP, q, 1)); nProtect++;
  double *vary0 = REAL(vary0_r);
  PROTECT(y0_r = allocMatrix(REALSXP, q, 1)); nProtect++;
  double *y0 = REAL(y0_r);

  if(verbose){
    Rprintf("-------------------------------------------------\n");
    Rprintf("\t\tPredicting\n");
    Rprintf("-------------------------------------------------\n");
#ifdef Win32
    R_FlushConsole();
#endif
  }

  double *z = (double *) R_alloc(q, sizeof(double));
  GetRNGstate();
  for(i = 0; i < q; i++){
    z[i] = rnorm(0.0,1.0);
  }
  PutRNGstate();




#ifdef _OPENMP
#pragma omp parallel for num_threads(nThreads) private(threadID, phi, nu, sigmaSq, tauSq, k, l, d, info)
#endif

  // for each prediction point i from 0 to q-1
  for(i = 0; i < q; i++){ //3; i++){ //q; i++){
    // Rprintf("\n\n-------------------------------------------------\n");
    // Rprintf("-------------------------------------------------\n");
    
#ifdef _OPENMP
    threadID = omp_get_thread_num();
    // Rprintf("\n threadID: %i \n", threadID);
#endif
    phi = theta[phiIndx];
    if(corName == "matern"){
      nu = theta[nuIndx];
    }
    sigmaSq = theta[sigmaSqIndx];
    tauSq = theta[tauSqIndx];

    for(k = 0; k < m; k++){
      // Rprintf("\n i: %i, k: %i \n", i, k);
      // distance between pred point i and nearest obs neighbor k
      // d = dist2(coords[nnIndx0[i+q*k]], coords[n+nnIndx0[i+q*k]], coords0[i], coords0[q+i]);
      // if(i < 20){
      //   Rprintf("\n i: %i, i+q*k: %i, nnIndx0[i+q*k] %i \n", i, i+q*k, nnIndx0[i+q*k]);
      // }
      // k_ind = i+q*k; //this was when I had a matrix format
      int k_ind = m*i + k;
      // Rprintf("d(i, k): %f \n", d);
      // afv[nnIndx0[i+q*l]-1]
      // Rprintf("i %i, k %i, k_ind %i, nnWght0[%i] %f, nnDist0[%i] %f\n", i, k, k_ind, k_ind, nnWght0[k_ind], k_ind, nnDist0[k_ind]);
      c[threadID*m+k] = nnWght0[k_ind]*sigmaSq*spCor(nnDist0[k_ind], phi, nu, covModel, &bk[threadID*nb]);
      // if(i==7064){
      //   Rprintf("i %i, k %i, k_ind %i, nnWght0[%i] %f, nnDist0[%i] %f\n", i, k, k_ind, k_ind, nnWght0[k_ind], k_ind, nnDist0[k_ind]);
      // }
      for(l = 0; l < m; l++){
        // Rprintf("%i ", threadID*mm+l*m+k);
        // Rprintf("\n k,l: %i, %i", k, l);
        // distance between pred point i's nearest obs neighbors k and l
        // d = dist2(coords[nnIndx0[i+q*k]], coords[n+nnIndx0[i+q*k]],
        //           coords[nnIndx0[i+q*l]], coords[n+nnIndx0[i+q*l]]);
        
        // nnIndx0 used to be 1-indexed, but now it's already 0-indexed
        // Rprintf("\n nnIndx0: %i, %i", nnIndx0[i+q*k]-1, nnIndx0[i+q*l]-1);
        // k_and_l = n*(nnIndx0[k_ind]-1) + (nnIndx0[m*i+l]-1);
        
        // Rprintf("\n nnIndx0: %i, %i", nnIndx0[i+q*k], nnIndx0[i+q*l]);
        int k_and_l = n*(nnIndx0[k_ind]) + (nnIndx0[m*i+l]);
        
        // if(i < 2 && l < 5 && k < 5){
        //   Rprintf("%f ", d);
        // }
        // afv[nnIndx0[i+q*l]-1], afv[nnIndx0[i+q*k]-1]
        // Rprintf("\t l %i, nnIndx0: %i, %i, k_and_l %i\n", l, nnIndx0[k_ind]-1, nnIndx0[m*i+l]-1, k_and_l);
        // Rprintf("\t l %i, obsobs_wt[%i] %f, obsobs_d[%i] %f\n", l, k_and_l, obsobs_wt[k_ind], k_and_l, obsobs_d[k_ind]);
        C[threadID*mm+l*m+k] = obsobs_wt[k_and_l]*sigmaSq*spCor(obsobs_d[k_and_l], phi, nu, covModel, &bk[threadID*nb]);
        if(k == l){
          C[threadID*mm+l*m+k] += tauSq;
        }
        // if(i < 2 && l < 5 && k < 5){
          // Rprintf("%f ", C[threadID*mm+l*m+k]);
        // }
        // if(i==7064){ //n*(nnIndx0[k_ind]) + (nnIndx0[m*i+l])
        //   Rprintf("\t l %i, nnIndx0[%i] = %i, nnIndx0[%i] = %i, k_and_l = %i, obsobs_wt[%i] %f, obsobs_d[%i] %f\n", l, k_ind, nnIndx0[k_ind], m*i+l, nnIndx0[m*i+l], k_and_l, k_and_l, obsobs_wt[k_and_l], k_and_l, obsobs_d[k_and_l]);
        //   // Rprintf("C[%i] = %f ", threadID*mm+l*m+k, C[threadID*mm+l*m+k]);
        // }
      }

      // if(i < 2 && k < 5){
      //   Rprintf("\n");
      // }
    }
    // if(i <2){
    // Rprintf("\n");
    // }

    // Rprintf("\n %d", C[threadID*mm]);
    // Rprintf("\n %f %i", C[threadID*mm], m);

    // if(i == 7064){
    //   Rprintf("C:\n");
    //   for(k = 0; k < m; k++){
    //     for(l = 0; l < m; l++){
    //       Rprintf("%f ", C[threadID*mm + l*m + k]);
    //       // Rprintf("%i ", threadID*mm + l*m + k);
    //     }
    //     Rprintf("\n");
    //   }
    //   Rprintf("\n");
    // }
    // 
    // if(i == 7064){
    //   Rprintf("c:\n");
    //   for(k = 0; k < m; k++){
    //       Rprintf("%f ", c[threadID*m+k]);
    //   }
    //   Rprintf("\n");
    // }

    // Cholesky factor of C to find the inverse
    F77_NAME(dpotrf)(lower, &m, &C[threadID*mm], &m, &info FCONE);
    // Rprintf("\n %d", info);
    if(info != 0){
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        if(err_flag == 0){ err_flag = 1; }
      }
      continue;
    }

    // This line replaces C with inverse of C
    F77_NAME(dpotri)(lower, &m, &C[threadID*mm], &m, &info FCONE);
    if(info != 0){
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        if(err_flag == 0){ err_flag = 2; }
      }
      continue;
    }

    // if(i == 7064){
    //   Rprintf("C inverse:\n");
    //   for(k = 0; k < m; k++){
    //     for(l = 0; l < m; l++){
    //       Rprintf("%f ", C[threadID*mm + l*m + k]);
    //       // Rprintf("%i ", threadID*mm + l*m + k);
    //     }
    //     Rprintf("\n");
    //   }
    //   Rprintf("\n");
    // }
    
    // This line stores C^{-1} dot c as tmp_m (i.e. this is KK^{-1} as in Rasmussen and Williams)
    F77_NAME(dsymv)(lower, &m, &one, &C[threadID*mm], &m, &c[threadID*m], &inc, &zero, &tmp_m[threadID*m], &inc FCONE);
    
    // if(i == 7064){
    //   Rprintf("tmp_m:\n");
    //   for(k = 0; k < m; k++){
    //     Rprintf("%f ", tmp_m[threadID*m+k]);
    //   }
    //   Rprintf("\n");
    // }

    // By the end of this loop, d is the predictive mean, KK^{-1}f, for prediction point i
    // where f = y - X*beta (observed)
    d = 0;
    for(k = 0; k < m; k++){
      int k_ind = m*i + k;
      d += tmp_m[threadID*m+k]*(y[nnIndx0[k_ind]] - F77_NAME(ddot)(&p, &X[nnIndx0[k_ind]], &n, &beta[0], &inc));
      // if(i == 7064){
      //   Rprintf("d=%f\n", d);
      // }
    }
    // if(i < 10){
    //   Rprintf("\n %f %f \n", d, y[i] - F77_NAME(ddot)(&p, &X[nnIndx0[i]-1], &n, &beta[0], &inc));
    // }
    // Rprintf("%f ", d);
    // Rprintf("\n");
    // if(i > 9000 && i < 9010){

    // if(i < 10){
    //   for(k = 0; k < m; k++){
    //     Rprintf("%f ", tmp_m[threadID*m+k]);
    //   }
    //   Rprintf("\n");
    // }



    // This line adds X*beta back in to go from f (zero-mean, called d here) to y (predicted fish density)
    y0[i] = F77_NAME(ddot)(&p, &X0[i], &q, &beta[0], &inc) + d;
    // if(i == 7064){
    //   Rprintf("y0[i]=%f\n", y0[i]);
    // }

    // if(i < 20){
    //   Rprintf("\n %i %f %f \n", i, F77_NAME(ddot)(&p, &X[i], &n, &beta[0], &inc), F77_NAME(ddot)(&p, &X0[i], &q, &beta[0], &inc));
    // }

    // if(i < 10){
    //   Rprintf("\n %f %f %f \n", d, y[i] - F77_NAME(ddot)(&p, &X0[i], &q, &beta[0], &inc), y0[i] - F77_NAME(ddot)(&p, &X0[i], &q, &beta[0], &inc));
    // }
    // if(i == 9000){
    //   Rprintf("\n\n\n %f %f %f \n", d, y[i] - F77_NAME(ddot)(&p, &X0[i], &q, &beta[0], &inc), y0[i] - F77_NAME(ddot)(&p, &X0[i], &q, &beta[0], &inc));
    // }
    // if(i > 9001 && i < 9010){
    //   Rprintf("\n %f %f %f \n", d, y[i] - F77_NAME(ddot)(&p, &X0[i], &q, &beta[0], &inc), y0[i] - F77_NAME(ddot)(&p, &X0[i], &q, &beta[0], &inc));
    // }
  // This line computes the predictive standard deviation;
  // tauSq is included here because BRISC models the response, not the latent process, as an NNGP
  vary0[i] = std::max(0.0, sqrt(sigmaSq + tauSq - F77_NAME(ddot)(&m, &tmp_m[threadID*m], &inc, &c[threadID*m], &inc)));
  }

  if(err_flag == 1){error("c++ error: dpotrf failed\n");}
  if(err_flag == 2){error("c++ error: dpotri failed\n");}

  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 1 + 1;

  PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

  SET_VECTOR_ELT(result_r, 0, y0_r);
  SET_VECTOR_ELT(resultName_r, 0, mkChar("p.y.0"));

  SET_VECTOR_ELT(result_r, 1, vary0_r);
  SET_VECTOR_ELT(resultName_r, 1, mkChar("var.y.0"));

  namesgets(result_r, resultName_r);

  //unprotect
  UNPROTECT(nProtect);

  return(result_r);
  
  // SEXP res = PROTECT(allocVector(REALSXP, 1)); nProtect++; res = 0;
  // UNPROTECT(nProtect);
  // return(res);

  }
  
  
//   // for hash tables in C:
//   // https://stackoverflow.com/questions/3578083/what-is-the-best-way-to-use-a-hashmap-in-c
//   
//   // here's how I would write this C function to take inputs from and send output to R
//   // later I might prefer to have it take inputs from/send output to C/C++
//   // inputs:
//   // pred_binID_r:      a single binary ID (string) whose nearest observed neighbors we want to find
//   // obs_binIDs_r:      vector of all binary IDs in the observed points data frame
//   // preds_binIDs_r:    vector of all binary IDs in the prediction points data frame
//   // streams_binIDs_r:  vector of all binary IDs in the stream segments data frame
//   // preds_updists_r:   vector of upstream distances for the prediction points
//   // streams_sums_r:    vector of upstream distance to the upstream node (i.e. upstream distance +
//   //                       segment length) of each stream segment
//   SEXP breadth_first_single_obs(SEXP pred_binID_r, SEXP obs_binIDs_r, SEXP preds_binIDs_r, SEXP preds_updists_r, 
//                                 SEXP streams_binIDs_r, SEXP streams_sums_r, SEXP m_r){
//     
//     // auto start_time = high_resolution_clock::now();
//     
//     // test code:
//     // SEXP out = PROTECT(allocVector(STRSXP, 3));
//     // 
//     // SET_STRING_ELT(out, 0, mkChar("1010"));
//     // SET_STRING_ELT(out, 1, mkChar("01010"));
//     // SET_STRING_ELT(out, 2, mkChar("111"));
//     // 
//     // UNPROTECT(1);
//     // 
//     // return out;
//     
//     // initialize local variables
//     // std::string pred_binID = CHAR(STRING_ELT(pred_binID_r, 0));
//     const char* pred_binID = CHAR(asChar(pred_binID_r));
//     
//     const char** obs_binIDs = CHAR(asChar(obs_binIDs_r));
//     // std::string preds_binIDs = CHAR(preds_binIDs_r);
//     // std::string streams_binIDs = CHAR(streams_binIDs_r);
//     // 
//     // double *preds_updists = REAL(preds_updists_r);
//     // double *streams_sums = REAL(streams_sums_r);
//     
//     // int m = INTEGER(m_r)[0];
//     
//     Rprintf("pred_binID %s \n", pred_binID);
//     Rprintf("obs_binIDs %s %s \n", obs_binIDs[0], obs_binIDs[1]);
//     // Rprintf("preds_binIDs %s \n", preds_binIDs);
//     // Rprintf("streams_binIDs %s \n", streams_binIDs);
//     // Rprintf("preds_updists %f \n", preds_updists);
//     // Rprintf("streams_sums %f \n", streams_sums);
//     // Rprintf("m %i \n", m);
//     
//     return ScalarReal(0);
//   }

}
