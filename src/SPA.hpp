
#ifndef SPA_HPP
#define SPA_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

double K0(double t_x,
          arma::mat t_muMat,     // N x (J-1)
          arma::mat t_cMat,      // N x (J-1)
          double t_m1);          // sum(muMat * cMat)

arma::vec K12(double t_x,
              arma::mat t_muMat,
              arma::mat t_cMat,
              double t_m1);

Rcpp::List fastgetroot_K1(double t_Stat,
                          double t_initX,
                          double t_Ratio0,
                          arma::mat t_muMat,
                          arma::mat t_cMat,
                          double t_m1);

double fastGet_Saddle_Prob(double t_Stat,
                           double t_zeta,
                           double t_K2,
                           double t_Ratio0,
                           arma::mat t_muMat,
                           arma::mat t_cMat,
                           double t_m1,          // sum(muMat * cMat)
                           bool t_lowerTail);

// add partial normal approximation to speed up the SPA
Rcpp::List fastSaddle_Prob(double t_Stat,
                           double t_VarP,
                           double t_VarW,
                           double t_Ratio0,      // Ratio of variance (G==0)
                           arma::vec t_K1roots,
                           arma::vec t_adjGVec1, // N1 x 1, where N1 is length(G!=0)
                           arma::mat t_muMat1,   // N1 x (J-1)
                           arma::mat t_iRMat1);  // N1 x (J-1)

#endif
