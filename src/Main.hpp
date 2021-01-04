
#ifndef MAIN_HPP
#define MAIN_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

arma::fmat getSymmMat(arma::fmat& xMat1,    // n x p
                      arma::fmat& xMat2,    // n x p
                      int p);

arma::fmat getfmatMulti(arma::fmat& xMat1,  // n x p
                        arma::fmat& xMat2); // n x p

Rcpp::List MAIN_REGION(std::vector<std::string> t_MarkerReqstd,
                       double t_NonZero_cutoff,
                       double t_StdStat_cutoff,
                       int t_maxMarkers,
                       std::string t_outputFile,
                       double t_missingRate_cutoff,
                       double t_maxMAF_cutoff,
                       std::string t_kernel,
                       arma::vec t_wBeta);

Rcpp::List MAIN_MARKER(std::vector<std::string> t_MarkerReqstd,
                       double t_StdStat_cutoff,
                       double t_missingRate_cutoff,
                       double t_minMAF_cutoff,
                       int t_minMAC_cutoff,
                       double t_varRatio);

#endif
