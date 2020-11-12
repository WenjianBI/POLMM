
#ifndef POLMM_HPP
#define POLMM_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace POLMM{

class POLMMClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  // dimensions: sample size, number of categories, number of covariates
  int m_n, m_J, m_p;
  
  arma::mat m_XXR_Psi_RX;  // XXR_Psi_RX ( n x p )
  arma::mat m_XR_Psi_R;    // XR_Psi_R ( p x n ), sum up XR_Psi_R ( p x n(J-1) ) for each subject
  arma::vec m_RymuVec;     // n x 1: row sum of the n x (J-1) matrix R %*% (yMat - muMat)
  arma::vec m_RPsiR;
  
  arma::mat m_iSigmaX_XSigmaX; // n(J-1) x p
  arma::mat m_CovaMat;     // n(J-1) x p
  
  arma::mat m_Cova, m_yMat;
  
  arma::sp_mat m_SparseGRM;
  arma::mat m_muMat;  // n x J
  arma::mat m_iRMat;
  
  arma::cube m_InvBlockDiagSigma;
  
  double m_tolPCG;
  int m_maxiterPCG;
  
  double m_tau;
  bool m_printPCGInfo;
  
  // for Efficient Resampling (ER)
  arma::umat m_SeqMat;

  ////////////////////// -------------------- functions ---------------------------------- //////////////////////
  

public:
  
  POLMMClass(arma::mat t_muMat,
             arma::mat t_iRMat,
             arma::mat t_Cova,
             arma::vec t_yVec,
             arma::sp_mat t_SparseGRM,
             double t_tau,
             bool t_printPCGInfo,
             double t_tolPCG,
             int t_maxiterPCG);
  
  arma::vec getadjGFast(arma::vec t_GVec);
  double getStatFast(arma::vec t_adjGVec);
  arma::vec get_ZPZ_adjGVec(arma::vec t_adjGVec);
  void getPCGofSigmaAndCovaMat(arma::mat t_xMat,             // matrix with dim of n(J-1) x p
                                           arma::mat& t_iSigma_xMat);    // matrix with dim of n(J-1) x p
  void getPCGofSigmaAndVector(arma::vec t_y1Vec,    // vector with length of n(J-1)
                              arma::vec& t_xVec);    // vector with length of n(J-1)
  arma::mat getSigmaxMat(arma::mat& t_xMat);
  
  arma::mat getiPsixMat(arma::mat t_xMat);
  arma::mat getPsixMat(arma::mat t_xMat);
  
  arma::cube getInvBlockDiagSigma();
  arma::mat solverBlockDiagSigma(arma::mat& t_xMat);

  void setRPsiR();
  arma::vec getVarWVec(arma::vec adjGVec);
  
  Rcpp::List MAIN_SPA(double t_Stat,
                      arma::vec t_adjGVec,
                      arma::vec t_K1roots,
                      double t_VarP,
                      double t_VarW,
                      double t_Ratio0,
                      arma::uvec t_posG1);
    
};

double K0(double t_x,
          arma::mat t_muMat,     // N x (J-1)
          arma::mat t_cMat,      // N x (J-1)
          double t_m1);           // sum(muMat * cMat)

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

}

#endif
