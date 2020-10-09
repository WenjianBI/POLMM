#ifndef POLMMGENE_H
#define POLMMGENE_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// int global_test = 1987;
// 
// // [[Rcpp::export]]
// int get_global_test(){
//   return global_test;
// }
// 
// // [[Rcpp::export]]
// void set_global_test(int t_global_test){
//   global_test = t_global_test;
// }

namespace POLMMGENE{

using namespace Rcpp;
using namespace std;

class POLMMGENEClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  // dimension: sample size, number of categories, number of covariates
  int m_n, m_J, m_p;
  
  // parameters
  double m_tau;
  bool m_showInfo;
  
  // control paramters
  int m_maxiterPCG; 
  double m_tolPCG;
  
  // SparseGRM 
  bool m_flagSparseGRM;
  Rcpp::List m_SparseGRM;
  arma::sp_mat m_SparseGRM_all;
  
  // working vectors/matrix
  arma::vec m_eta;
  arma::mat m_muMat, m_iRMat, m_iSigmaX_XSigmaX;
  
  // input data
  arma::uvec m_yVec; // should start from 1, not 0
  arma::mat m_Cova;
  
  // almost input data
  arma::mat m_yMat, m_CovaMat;
  
  //
  Rcpp::List m_objP, m_LOCOList;
  string m_excludechr;
  arma::vec m_RPsiRVec;
    
  // for Efficient Resampling (ER)
  arma::umat m_SeqMat;

  ////////////////////// -------------------- functions ---------------------------------- //////////////////////
  
  // yMat: matrix with dim of n x J
  arma::mat getyMat();
  arma::vec get_ZPZ_adjGVec(arma::vec t_adjGVec, string t_excludechr);
  void getPCGofSigmaAndVector(arma::vec t_y1Vec, arma::vec& t_xVec, string t_excludechr);
  void getPCGofSigmaAndCovaMat(arma::mat t_xMat, arma::mat& t_iSigma_xMat, string t_excludechr);
  arma::mat solverBlockDiagSigma(arma::cube& InvBlockDiagSigma,   // (J-1) x (J-1) x n
                                 arma::mat& xMat);
  arma::cube getInvBlockDiagSigma();
  arma::mat getSigmaxMat(arma::mat t_xMat,   // matrix: n x (J-1) 
                         string t_excludechr);
  arma::vec getKinbVecPOLMM(arma::vec t_bVec, string t_excludeChr);
  arma::vec tZMat(arma::vec t_xVec);
  arma::vec ZMat(arma::vec t_xVec);
  arma::mat getiPsixMat(arma::mat t_xMat);
  
public:
  
  POLMMGENEClass(int t_maxiterPCG,
                 double t_tolPCG,
                 arma::mat t_Cova,
                 arma::uvec t_yVec,     // should be from 1 to J
                 double t_tau,
                 Rcpp::List t_SparseGRM,    // results of function getKinMatList()
                 Rcpp::List t_LOCOList,
                 arma::vec t_eta,
                 int t_nMaxNonZero);
  
  void setPOLMMGENEchr(Rcpp::List t_LOCOList, string t_excludechr);
  
  Rcpp::List getStatVarS(arma::mat t_GMat, 
                         double t_NonZero_cutoff,
                         double t_StdStat_cutoff);
  
  // get single marker p values from ER or SPA
  double getPvalERinClass(arma::vec t_GVec);
  
  void check_ZPZ_adjGVec(arma::vec t_adjGVec);
  
  double checkError();
  
};

}


#endif
