
#ifndef ERSPA_H
#define ERSPA_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::umat makeSeqMat(int t_n,         // number of subjects with G != 0 
                      int t_J);         // number of levels

// [[Rcpp::export]]
arma::umat updateSeqMat(arma::umat t_SeqMat, // n x J^n matrix
                        int t_n1,            // number of subjects, should be < n
                        int t_J);             // number of levels

// [[Rcpp::export]]
arma::vec getStatVec(arma::umat t_SeqMat,   // n x J^n matrix
                     arma::vec t_GVec,      // n x 1 vector, where n is number of subjects with Geno != 0
                     arma::mat t_muMat,     // n x J matrix, where n is number of subjects with Geno != 0
                     arma::mat t_iRMat);     // n x (J-1) matrix

double getProbOne(arma::uvec t_SeqVec,  // n x 1
                  arma::mat t_muMat);    // n x J

// [[Rcpp::export]]
double getProb(arma::umat t_SeqMat,          // n x m matrix, where m \leq J^n is the number of resampling with abs(stat) > stat_obs
               arma::mat t_muMat);            // n x J matrix

double getPvalER(arma::uvec t_yVec,     // n x 1 vector, from 0 to J-1
                 arma::vec t_GVec,      // n x 1 vector,
                 arma::mat t_muMat,     // n x J matrix,
                 arma::mat t_iRMat);     // n x (J-1) matrix


#endif
