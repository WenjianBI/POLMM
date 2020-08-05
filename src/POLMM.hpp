#ifndef POLMM_H
#define POLMM_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "DenseGRM.hpp"
#include "Plink.hpp"
#include "SubFunc.hpp"
#include "POLMM.hpp"

namespace POLMM {

using namespace Rcpp;
using namespace std;
using namespace Plink;
using namespace DenseGRM;

class POLMMClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  // flags to indicate sparse/dense GRM, using GMatRatio or Plink to estimate variance ratio
  bool m_flagSparseGRM, m_flagGMatRatio;
  
  // dimension: sample size, number of categories, number of covariates, markers in Plink files
  int m_n, m_J, m_p, m_M;
  
  // input data
  arma::Col<int> m_yVec; // should start from 1, not 0
  arma::mat m_Cova, m_GMatRatio;
  
  // almost input data
  arma::mat m_yMat, m_CovaMat;
  
  // parameters
  arma::vec m_beta, m_bVec, m_eps;
  double m_tau;
  int m_iter;
  
  // control paramters
  int m_maxiter, m_maxiterPCG, m_maxiterEps; 
  int m_tracenrun, m_seed, m_nSNPsVarRatio, m_grainSize;
  double m_tolBeta, m_tolTau, m_tolPCG, m_tolEps; 
  double m_CVcutoff, m_minMafVarRatio, m_maxMissingVarRatio;
  double m_memoryChunk, m_minMafGRM, m_maxMissingGRM;
  bool m_LOCO, m_showInfo;
  Rcpp::List m_LOCOList;
  
  // working vectors/matrix
  arma::mat m_WMat, m_muMat, m_mMat, m_nuMat, m_iRMat, m_YMat; 
  arma::mat m_iSigma_CovaMat, m_iSigmaX_XSigmaX;
  arma::vec m_eta, m_iSigma_YVec, m_iSigma_VPYVec;
  
  arma::mat m_TraceRandMat, m_V_TRM, m_iSigma_V_TRM;
  
  PlinkClass* m_ptrPlinkObj;
  DenseGRMClass* m_ptrDenseGRMObj;
  
  Rcpp::List m_initParList, m_controlList;
  
  // SparseGRM 
  Rcpp::List m_SparseGRM;
  arma::sp_mat m_SparseGRM_all;
  arma::sp_mat m_ZMat_sp;
  arma::sp_mat m_SigmaMat_sp;
  
  ////////////////////// -------------------- functions ---------------------------------- //////////////////////
  
  // functions
  void setSigmaMat_sp();
  
  // functions in setNullModel()
  arma::mat getyMat();
  
  // set up m_TraceRandMat (TRM) and m_V_TRM, only used once at setPOLMMObj()
  void getTraceRandMat();
  
  // sum up each (J-1) elements: n(J-1) x 1 -> n x 1
  arma::vec ZMat(arma::vec t_xVec); 
  
  // duplicate each element for (J-1) times: n x 1 -> n(J-1) x 1 
  arma::vec tZMat(arma::vec t_xVec);
  
  // update (m_eta, m_WMat, m_muMat, m_mMat, m_nuMat, m_iRMat, m_YMat) based on (m_beta, m_bVec, m_eps)
  void updateMats();
  
  // outMat = iPsiMat %*% xMat, iPsiMat is determined by muMat
  arma::mat getiPsixMat(arma::mat t_xMat);
  
  // use PCG to calculate iSigma_xMat = Sigma^-1 %*% xMat
  void getPCGofSigmaAndCovaMat(arma::mat t_xMat,              // matrix with dim of n(J-1) x p
                               arma::mat& t_iSigma_xMat,      // matrix with dim of n(J-1) x p
                               string t_excludechr);
  
  // use PCG to calculate xVec = Sigma^-1 %*% yVec
  void getPCGofSigmaAndVector(arma::vec t_y1Vec,    // vector with length of n(J-1)
                              arma::vec& t_xVec,    // vector with length of n(J-1)
                              string t_excludechr);
  
  // update parameters (except tau) until converge
  void updateParaConv(string t_excludechr);
  void updatePara(string t_excludechr);
  
  // 
  void updateEps();
  void updateEpsOneStep();
  
  // 
  void updateTau();
  
  // yMat = Sigma %*% xMat
  arma::mat getSigmaxMat(arma::mat t_xMat,   // matrix: n x (J-1) 
                         string t_excludechr);
  
  // used in getPCGofSigmaAndVector()
  arma::cube getInvBlockDiagSigma();
  
  arma::mat solverBlockDiagSigma(arma::cube& InvBlockDiagSigma,   // (J-1) x (J-1) x n
                                 arma::mat& xMat);                 // n x (J-1)
  
  
  // functions if m_LOCO = TRUE
  Rcpp::List getLOCO(string t_excludechr,
                     std::vector<int> t_indexSNPs);

  arma::rowvec getVarOneSNP(arma::vec GVec,
                            string excludechr,
                            Rcpp::List objP);
  
  double getVarP(arma::vec t_adjGVec, string t_excludechr);
  
  // functions if LOCO = FALSE
  arma::mat getVarRatio(arma::mat t_GMatRatio, string t_excludechr);
  
  void setPOLMMInner(arma::mat t_Cova,
                     arma::Col<int> t_yVec,     // should be from 1 to J
                     arma::vec t_beta,
                     arma::vec t_bVec,
                     arma::vec t_eps,           // 
                     double t_tau);
  
  void setControlList(Rcpp::List t_controlList);
  void setArray();
  
  arma::sp_mat setZMat_sp();
    
public:
  
  void setPOLMMObj(bool t_flagSparseGRM,       // if 1, then use SparseGRM, otherwise, use DenseGRM
                   bool t_flagGMatRatio,       // if 1, then use GMatRatio, otherwise, extract from Plink files
                   PlinkClass* t_ptrPlinkObj,
                   DenseGRMClass* t_ptrDenseGRMObj,
                   arma::mat t_Cova,
                   arma::Col<int> t_yVec,     // should be from 1 to J
                   arma::vec t_beta,
                   arma::vec t_bVec,
                   arma::vec t_eps,           // 
                   double t_tau,
                   arma::mat t_GMatRatio,
                   Rcpp::List t_SparseGRM,
                   Rcpp::List t_controlList);
  
  arma::vec getKinbVecPOLMM(arma::vec t_bVec, string t_excludeChr);
  
  void setNullModel(std::string,
                    arma::ivec,
                    arma::mat, arma::Col<int>, Rcpp::List, arma::vec, arma::vec, arma::vec, double, arma::mat, Rcpp::List);
  void fitPOLMM();
  
  Rcpp::List getPOLMM();
  
  void closeGenoObj();
  
};

}

#endif
