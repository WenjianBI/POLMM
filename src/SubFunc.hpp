
#ifndef SUBFUNC_H
#define SUBFUNC_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;

// count number of lines of a text file
int countLine(string t_file);

// record computation time
arma::vec getTime();
void printTime(arma::vec t1, arma::vec t2, std::string message);

// freq --> invStd
double getinvStd(double t_freq);

// set random seed
void set_seed(unsigned int seed);

arma::vec nb(int n);

double calCV(arma::vec t_xVec);

// C++ version of which(). Note: start from 0, not 1 
std::vector<int> whichCPP(Rcpp::StringVector strVec, 
                          std::string strValue);


arma::mat getCovaMat(arma::mat Cova,   // matrix: n x p
                     int n, int J, int p); 

arma::mat convert2(arma::vec xVec, // n(J-1) x 1 
                   int n, int J);

arma::vec convert1(arma::mat xMat, // matrix: n x (J-1)
                   int n, int J);

arma::vec getRowSums(arma::mat t_xMat);

Rcpp::List getobjP(arma::mat t_Cova,     // matrix: n x p
                   arma::mat t_yMat,
                   arma::mat t_muMat,    // matrix: n x J
                   arma::mat t_iRMat);   // matrix: n x (J-1)

arma::mat getPsixMat(arma::mat t_xMat,    // matrix: n x (J-1)
                     arma::mat t_muMat);   // matrix: n x J

arma::mat sumCols(arma::mat t_xMat,
                  int J);

arma::vec getRPsiR(arma::mat t_muMat,
                   arma::mat t_iRMat,
                   int t_n, int t_J, int t_p);

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat);

Rcpp::List outputadjGFast(arma::vec GVec,
                          Rcpp::List objP);

arma::vec getadjGFast(arma::vec GVec,
                      arma::mat XXR_Psi_RX_new,   // XXR_Psi_RX_new ( n x p )
                      arma::mat XR_Psi_R_new,     // XR_Psi_R_new ( p x n ), sum up XR_Psi_R ( p x n(J-1) ) for each subject 
                      int n, int p);

double getStatFast(arma::vec GVec,         // n x 1
                   arma::vec RymuVec);     // n x 1: row sum of the n x (J-1) matrix R %*% (yMat - muMat)

double getVarWFast(arma::vec adjGVec,     // n x 1
                   arma::vec RPsiRVec);   // n x 1
              

arma::mat getyMatR(arma::mat yVec, int n, int J);
  
  
#endif
