
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

#include "ER_SPA.hpp"

using namespace Rcpp;


arma::umat makeSeqMat(int t_n,         // number of subjects with G != 0 
                      int t_J)         // number of levels
{
  int nER = pow(t_J, t_n);       // J^n
  arma::uvec y = arma::linspace<arma::uvec>(0, nER-1, nER);  // nER x 1 matrix: seq(0,nER-1,1)
  arma::umat SeqMat(t_n, nER);
  int powJ = nER / t_J;         // J^(n-1)
  for(int i = 0; i < t_n; i++){
    int pos_row = t_n - 1 - i;
    arma::uvec SeqVec = y / powJ;  // nER x 1
    y = y - SeqVec * powJ;
    powJ = powJ / t_J;           // J^(n-2), ..., J^0
    SeqMat.row(pos_row) = SeqVec.t();
  }
  
  return SeqMat;
}


arma::umat updateSeqMat(arma::umat t_SeqMat, // n x J^n matrix
                        int t_n1,            // number of subjects, should be < n
                        int t_J)             // number of levels
{
  int nER = pow(t_J, t_n1);
  arma::umat PartSeqMat = t_SeqMat.submat(0, 0, t_n1-1, nER-1);
  return PartSeqMat;
}


arma::vec getStatVec(arma::umat t_SeqMat,   // n x J^n matrix
                     arma::vec t_GVec,      // n x 1 vector, where n is number of subjects with Geno != 0
                     arma::mat t_muMat,     // n x J matrix, where n is number of subjects with Geno != 0
                     arma::mat t_iRMat)     // n x (J-1) matrix
{
  int n = t_muMat.n_rows;
  int J = t_muMat.n_cols;
  int nER = t_SeqMat.n_cols;
  
  arma::vec StatVec(nER);
  
  arma::mat A(n, J-1);
  for(int i = 0; i < J-1; i++){
    A.col(i) = t_GVec / t_iRMat.col(i);
  }
  
  double a1 = arma::accu(A % t_muMat.cols(0, J-2));
  
  for(int i = 0; i < nER; i++){
    double a2 = 0;
    for(int j = 0; j < n; j++){
      int idxL = t_SeqMat(j, i); // from 0 to J-1, level index
      if(idxL != J-1){
        a2 += A(j,idxL);
      }
    }
    
    StatVec(i) = a2 - a1;
  }
  
  return StatVec;
}

double getProbOne(arma::uvec t_SeqVec,  // n x 1
                  arma::mat t_muMat)    // n x J
{
  int n = t_muMat.n_rows;
  double tempProb = 1;
  for(int j = 0; j < n; j++){
    tempProb *= t_muMat(j, t_SeqVec(j));
  }
  return tempProb;
}


double getProb(arma::umat t_SeqMat,          // n x m matrix, where m \leq J^n is the number of resampling with abs(stat) > stat_obs
               arma::mat t_muMat)            // n x J matrix
{
  int nER = t_SeqMat.n_cols;
  
  double prob = 0;
  
  for(int i = 0; i < nER; i++){
    arma::uvec t_SeqVec = t_SeqMat.col(i);
    double tempProb = getProbOne(t_SeqVec, t_muMat);
    prob += tempProb;
  }
  
  return prob;
}

// Main function: note that n is the number of subjects with Geno != 0
double getPvalER(arma::uvec t_yVec,     // n x 1 vector, from 1 to J
                 arma::vec t_GVec,      // n x 1 vector,
                 arma::mat t_muMat,     // n x J matrix,
                 arma::mat t_iRMat)     // n x (J-1) matrix
{
  int n = t_muMat.n_rows;
  int J = t_muMat.n_cols;
  
  arma::uvec yVec(n);  // yVec: from 0 to J-1; t_yVec: from 1 to J
  for(int i = 0; i < n; i++){
    yVec(i) = t_yVec(i) - 1;
  }
  
  arma::umat SeqMat = makeSeqMat(n, J);
  int nER = SeqMat.n_cols;
  arma::vec StatVec = getStatVec(SeqMat, t_GVec, t_muMat, t_iRMat);
  double StatObs = arma::as_scalar(getStatVec(yVec, t_GVec, t_muMat, t_iRMat));
  
  double eps = 1e-10;
  
  double pvalER = 0;
  double absStatObs = std::abs(StatObs);
  for(int i = 0; i < nER; i++){
    double absStatTmp = std::abs(StatVec(i));
    if(absStatObs < absStatTmp - eps){
      pvalER += getProb(SeqMat.col(i), t_muMat);
    }else if(absStatObs < absStatTmp + eps){
      pvalER += 0.5 * getProb(SeqMat.col(i), t_muMat);
    }
  }
  
  // arma::uvec idxVec1 = arma::find(abs(StatVec) > std::abs(StatObs) + eps); // abs(double) return int;
  // double pvalER1 = getProb(SeqMat.cols(idxVec1), t_muMat);
  // 
  // arma::uvec idxVec2 = arma::find(abs(StatVec) > std::abs(StatObs) - eps && abs(StatVec) < std::abs(StatObs) + eps); // abs(double) return int;
  // double pvalER2 = getProb(SeqMat.cols(idxVec2), t_muMat);
  // 
  // double pvalER = pvalER1 + 0.5 * pvalER2;
  // double pvalER = 1 - pvalER
  
  // Rcpp::List outList = Rcpp::List::create(Named("StatVec") = StatVec,
  //                                         Named("StatObs") = StatObs,
  //                                         Named("idxVec1") = idxVec1,
  //                                         Named("idxVec2") = idxVec2,
  //                                         Named("pvalER") = pvalER);
  // return outList;
  return pvalER;
}
