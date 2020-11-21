
#ifndef UTIL_HPP
#define UTIL_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

double getWeights(std::string t_kernel, 
                  double t_freq, 
                  arma::vec t_wBeta);

void imputeGeno(arma::vec& GVec, 
                double freq, 
                std::vector<uint32_t> posMissingGeno);

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat);

// duplicate each element for (J-1) times: n x 1 -> n(J-1) x 1 
arma::vec Vec2LongVec(arma::vec t_xVec, int n, int J);

// sum up each (J-1) elements: n(J-1) x 1 -> n x 1
arma::vec LongVec2Vec(arma::vec t_xVec, int n, int J);

// convert: n(J-1) x 1 -> n x (J-1) 
arma::mat Vec2Mat(arma::vec xVec, int n, int J);

// convert: n x (J-1) -> n(J-1) x 1
arma::vec Mat2Vec(arma::mat xMat, int n, int J);

// duplicate each row for (J-1) times: n x p -> n(J-1) x p
arma::mat getCovaMat(arma::mat Cova, int n, int J, int p);

arma::mat sumCols(arma::mat t_xMat, int J);

arma::vec getRPsiR(arma::mat t_muMat, arma::mat t_iRMat, int t_n, int t_J, int t_p); 

#endif

