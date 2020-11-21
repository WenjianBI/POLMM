
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
  
  arma::uvec m_yVec;
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
  arma::Mat<uint8_t> m_SeqMat;
  
  ////////////////////// -------------------- functions ---------------------------------- //////////////////////
  
  
public:
  
  POLMMClass(arma::mat t_muMat,
             arma::mat t_iRMat,
             arma::mat t_Cova,
             arma::uvec t_yVec,
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
  
  double MAIN_ER(arma::vec t_GVec,
                 arma::uvec t_posG1);
  
  void setSeqMat(int t_NonZero_cutoff);
  
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

double getPvalER(arma::uvec t_yVec,     // N1 x 1 vector, from 0 to J-1
                 arma::vec t_GVec,      // N1 x 1 vector,
                 arma::mat t_muMat,     // N1 x J matrix,
                 arma::mat t_iRMat,     // N1 x (J-1) matrix
                 arma::Mat<uint8_t> t_SeqMat);   // N1 x nER

arma::vec getStatVec(arma::Mat<uint8_t> t_SeqMat,   // n x J^n matrix
                     arma::vec t_GVec,         // n x 1 vector, where n is number of subjects with Geno != 0
                     arma::mat t_muMat,        // n x J matrix, where n is number of subjects with Geno != 0
                     arma::mat t_iRMat);       // n x (J-1) matrix

double getProbOne(arma::Col<uint8_t> t_SeqVec,  // n x 1
                  arma::mat t_muMat);           // n x J

double getProb(arma::Mat<uint8_t> t_SeqMat,  // n x m matrix, where m \leq J^n is the number of resampling with abs(stat) > stat_obs
               arma::mat t_muMat);           // n x J matrix

}

#endif
