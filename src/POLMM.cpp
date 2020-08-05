
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <string>

#include "DenseGRM.hpp"
#include "Plink.hpp"
#include "SubFunc.hpp"
#include "POLMM.hpp"

namespace POLMM {

using namespace Rcpp;
using namespace std;
using namespace Plink;
using namespace DenseGRM;

void POLMMClass::setPOLMMObj(bool t_flagSparseGRM,       // if 1, then use SparseGRM, otherwise, use DenseGRM
                             bool t_flagGMatRatio,       // if 1, then use GMatRatio, otherwise, extract from Plink files
                             PlinkClass* t_ptrPlinkObj,
                             DenseGRMClass* t_ptrDenseGRMObj,
                             arma::mat t_Cova,
                             arma::Col<int> t_yVec,     // should be from 1 to J
                             arma::vec t_beta,
                             arma::vec t_bVec,
                             arma::vec t_eps,           // 
                             double t_tau,
                             arma::mat t_GMatRatio,     // only used if m_LOCO = FALSE
                             Rcpp::List t_SparseGRM,    // results of function getKinMatList()
                             Rcpp::List t_controlList)
{
  setControlList(t_controlList);
  setPOLMMInner(t_Cova, t_yVec, t_beta,  t_bVec,  t_eps,  t_tau);
  
  m_ptrPlinkObj = t_ptrPlinkObj;
  m_ptrDenseGRMObj = t_ptrDenseGRMObj;
  
  // if t_flagSparseGRM = 1, then use "SparseGRM" methods, otherwise, use "DenseGRM" methods
  m_flagSparseGRM = t_flagSparseGRM;
  if(m_flagSparseGRM){
    m_SparseGRM = t_SparseGRM;
    arma::sp_mat temp = m_SparseGRM["none"];
    m_SparseGRM_all = temp;
    m_ZMat_sp = setZMat_sp();
    m_M = 0;
  }else{
    m_M = t_ptrDenseGRMObj->getM();
  }
  
  // if t_flagGMatRatio = 1, then use "GMatRatio", otherwise, extract "GMatRatio" from Plink files
  m_flagGMatRatio = t_flagGMatRatio;
  if(m_flagGMatRatio){
    m_GMatRatio = t_GMatRatio;
  }
  
  getTraceRandMat();
}

void POLMMClass::fitPOLMM()
{
  // initial vector
  arma::vec t1  = getTime();
  updateMats();
  
  // start iteration
  cout << "Start iteration ....." << endl;
  
  for(m_iter = 0; m_iter < m_maxiter; m_iter ++){
    
    // update fixed effect coefficients
    updateParaConv("none");
    
    // update tau
    double tau0 = m_tau;
    updateTau();
    
    if(std::isnan(m_tau))
      stop("Parameter tau is NA.");
    
    cout << "iter: " << m_iter << endl;
    cout << "beta: " << endl << m_beta << endl;
    cout << "tau: " << m_tau << endl << endl;
    
    double diffTau = abs(m_tau - tau0) / (abs(m_tau) + abs(tau0) + m_tolTau);
    if(diffTau < m_tolTau)
      break;
  }
  
  if(m_LOCO){
    
    // turn on LOCO option
    Rcpp::StringVector chrVec = m_ptrPlinkObj->getChrVec();
    Rcpp::StringVector uniqchr = unique(chrVec);
    
    if(m_flagSparseGRM){
      Rcpp::StringVector uniqchr_sp = m_SparseGRM.names(); 
      uniqchr = intersect(uniqchr, uniqchr_sp);
    }
    
    cout << "uniqchr is " << uniqchr << endl;
    
    for(int i = 0; i < uniqchr.size(); i ++){
      
      
      string excludechr = string(uniqchr(i));
      cout << endl << "Leave One Chromosome Out: Chr " << excludechr << endl;
      
      updateParaConv(excludechr);

      m_GMatRatio = m_ptrPlinkObj->getGMat(100, excludechr, m_minMafVarRatio, m_maxMissingVarRatio);
      arma::mat VarRatioMat = getVarRatio(m_GMatRatio, excludechr);
      double VarRatio = arma::mean(VarRatioMat.col(4));
      
      Rcpp::List temp = List::create(Named("muMat") = m_muMat,
                                     Named("iRMat") = m_iRMat,
                                     Named("VarRatioMat") = VarRatioMat,
                                     Named("VarRatio") = VarRatio);
      
      m_LOCOList[excludechr] = temp;
    }
    
  }else{
    
    // turn off LOCO option
    if(!m_flagGMatRatio){
      m_GMatRatio = m_ptrPlinkObj->getGMat(100, "none", m_minMafVarRatio, m_maxMissingVarRatio);
    }
    arma::mat VarRatioMat = getVarRatio(m_GMatRatio, "none");
    double VarRatio = arma::mean(VarRatioMat.col(4));
    
    Rcpp::List temp = List::create(Named("muMat") = m_muMat,
                                   Named("iRMat") = m_iRMat,
                                   Named("VarRatioMat") = VarRatioMat,
                                   Named("VarRatio") = VarRatio);
    m_LOCOList["LOCO=F"] = temp;
  }
  
  // complete null POLMM fitting 
  arma::vec t2  = getTime();
  printTime(t1, t2, "fit the null POLMM.");
}

arma::mat POLMMClass::getVarRatio(arma::mat t_GMatRatio, string t_excludechr)
{
  cout << "Start estimating variance ratio...." << endl;
  Rcpp::List objP = getobjP(m_Cova, m_yMat, m_muMat, m_iRMat);
  
  arma::vec GVec(m_n);
  arma::rowvec VarOneSNP(5);
  
  arma::mat VarRatioMat(m_nSNPsVarRatio, 5);
  arma::mat newVarRatio(10, 5);
  
  int index = 0;
  int indexTot = 0;
  while(index < m_nSNPsVarRatio){
    GVec = t_GMatRatio.col(index);
    VarOneSNP = getVarOneSNP(GVec, t_excludechr, objP);
    VarRatioMat.row(index) = VarOneSNP;
    index++;
    indexTot++;
  }
  
  arma::vec VarRatio = VarRatioMat.col(4);
  double CV = calCV(VarRatio);
  cout << "nSNPs for CV: " << index << endl;
  cout << "CV: " << CV << endl;
  
  while(CV > m_CVcutoff && VarRatioMat.n_rows <= 100){
    int indexTemp = 0;
    while(indexTemp < 10){
      indexTot++;
      GVec = t_GMatRatio.col(indexTot);
      VarOneSNP = getVarOneSNP(GVec, t_excludechr, objP);
      newVarRatio.row(indexTemp) = VarOneSNP;
      index++;
      indexTemp++;
    }
    VarRatioMat.insert_rows(0, newVarRatio);
    arma::vec VarRatio = VarRatioMat.col(4);
    CV = calCV(VarRatio);
    cout << "nSNPs for CV: " << index << endl;
    cout << "CV: " << CV << endl;
  }
  return(VarRatioMat);
}

arma::rowvec POLMMClass::getVarOneSNP(arma::vec GVec,
                                      string excludechr,
                                      Rcpp::List objP)
{
  arma::rowvec VarOut(5);
  double AF = sum(GVec) / GVec.size() / 2;
  if(AF > 0.5)
    AF = 1 - AF;
  
  Rcpp::List adjGList = outputadjGFast(GVec, objP);
  arma::vec adjGVec = adjGList["adjGVec"];
  double Stat = adjGList["Stat"];
  double VarW = adjGList["VarW"];
  double VarP = getVarP(adjGVec, excludechr);
  
  VarOut(0) = AF;
  VarOut(1) = Stat;
  VarOut(2) = VarW;
  VarOut(3) = VarP;
  VarOut(4) = VarP/VarW;
  return(VarOut);
}

// update parameters (except tau) until converge
void POLMMClass::updateParaConv(string t_excludechr)
{
  for(int iter = 0; iter < m_maxiter; iter ++){
    
    arma::vec beta0 = m_beta;
    
    // update beta and bVec
    updatePara(t_excludechr);
    updateMats();
    
    // update eps (cutpoints)
    updateEps();
    updateMats();
    
    cout << "beta: " << endl << m_beta << endl;
    
    double diffBeta = max(abs(m_beta - beta0)/(abs(m_beta) + abs(beta0) + m_tolBeta));
    cout << "diffBeta:\t" << diffBeta << endl << endl;
    if(diffBeta < m_tolBeta)
      break;
  }
}

void POLMMClass::updatePara(string t_excludechr)
{
  getPCGofSigmaAndCovaMat(m_CovaMat, m_iSigma_CovaMat, t_excludechr);
  arma::vec YVec = convert1(m_YMat, m_n, m_J);
  getPCGofSigmaAndVector(YVec, m_iSigma_YVec, t_excludechr); 
  
  // update beta
  arma::mat XSigmaX = inv(m_CovaMat.t() * m_iSigma_CovaMat);
  arma::vec Cova_iSigma_YVec = m_CovaMat.t() * m_iSigma_YVec;
  m_beta = XSigmaX * Cova_iSigma_YVec;
  m_iSigmaX_XSigmaX = m_iSigma_CovaMat * XSigmaX;
  
  // update bVec
  arma::vec Z_iSigma_YVec = ZMat(m_iSigma_YVec);
  arma::vec Z_iSigma_Xbeta = ZMat(m_iSigma_CovaMat * m_beta);
  arma::vec tempVec = Z_iSigma_YVec - Z_iSigma_Xbeta;
  m_bVec = m_tau * getKinbVecPOLMM(tempVec, t_excludechr);
}

// use PCG to calculate iSigma_xMat = Sigma^-1 %*% xMat
void POLMMClass::getPCGofSigmaAndCovaMat(arma::mat t_xMat,              // matrix with dim of n(J-1) x p
                                         arma::mat& t_iSigma_xMat,      // matrix with dim of n(J-1) x p
                                         string t_excludechr)
{
  int p1 = t_xMat.n_cols;
  for(int i = 0; i < p1; i++){
    arma::vec y1Vec = t_xMat.col(i);
    arma::vec iSigma_y1Vec = t_iSigma_xMat.col(i);
    getPCGofSigmaAndVector(y1Vec, iSigma_y1Vec, t_excludechr);
    t_iSigma_xMat.col(i) = iSigma_y1Vec;
  }
}

void POLMMClass::updateTau()
{
  cout << "Start updating tau..." << endl;
  arma::vec YVec = convert1(m_YMat, m_n, m_J);
  getPCGofSigmaAndCovaMat(m_CovaMat, m_iSigma_CovaMat, "none");
  getPCGofSigmaAndVector(YVec, m_iSigma_YVec, "none"); 
  m_iSigmaX_XSigmaX = m_iSigma_CovaMat * inv(m_CovaMat.t() * m_iSigma_CovaMat);
  arma::vec PYVec = m_iSigma_YVec - m_iSigmaX_XSigmaX * (m_CovaMat.t() * m_iSigma_YVec);
  arma::vec ZPYVec = ZMat(PYVec);
  arma::vec VPYVec = tZMat(getKinbVecPOLMM(ZPYVec, "none"));
  
  getPCGofSigmaAndVector(VPYVec, m_iSigma_VPYVec, "none");
  arma::vec PVPYVec = m_iSigma_VPYVec - m_iSigmaX_XSigmaX * (m_CovaMat.t() * m_iSigma_VPYVec);
  double YPVPY = as_scalar(YVec.t() * PVPYVec);
  double YPVPVPY = as_scalar(VPYVec.t() * PVPYVec);
  // The below is to calculate trace
  getPCGofSigmaAndCovaMat(m_V_TRM, m_iSigma_V_TRM, "none");
  double tracePV = 0;
  int m = m_TraceRandMat.n_cols;
  for(int i = 0; i < m; i++){
    arma::vec iSigma_V_TRM_col = m_iSigma_V_TRM.col(i);
    arma::vec P_V_TRM_col = iSigma_V_TRM_col - m_iSigmaX_XSigmaX * (m_CovaMat.t() * iSigma_V_TRM_col);
    tracePV += as_scalar(m_TraceRandMat.col(i).t() * P_V_TRM_col);
  }
  tracePV /= m;
  // final step
  double deriv = 0.5 * YPVPY - 0.5 * tracePV;
  double AI = 0.5 * YPVPVPY;
  double dtau = deriv / AI;
  double tau0 = m_tau;
  m_tau = tau0 + dtau;
  while(m_tau < 0){
    dtau = dtau / 2;
    m_tau = tau0 + dtau;
  }
  if(m_tau < 1e-4){
    m_tau = 0;
  }
}

// use PCG to calculate xVec = Sigma^-1 %*% yVec
void POLMMClass::getPCGofSigmaAndVector(arma::vec t_y1Vec,    // vector with length of n(J-1)
                                        arma::vec& t_xVec,    // vector with length of n(J-1)
                                        string t_excludechr)
{
  // if(m_flagSparseGRM){
    // setSigmaMat_sp();
    // t_xVec = spsolve(SigmaMat_sp, t_y1Vec);
  // }else{
    arma::mat xMat = convert2(t_xVec, m_n, m_J);
    arma::mat y1Mat = convert2(t_y1Vec, m_n, m_J);
    // r2Vec and z2Vec are for the current step; r1Vec and z1Vec are for the previous step
    int iter = 0;
    arma::mat r2Mat = y1Mat - getSigmaxMat(xMat, t_excludechr);  // n x (J-1): r0 = y1Mat- Sigma %*% xMat
    double meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
    if(meanL2 <= m_tolPCG){
      // do nothing, xMat is already close to (Sigma)^-1 %*% y1Mat
    }else{
      iter++;
      
      arma::cube InvBlockDiagSigma = getInvBlockDiagSigma();
      arma::mat z2Mat = solverBlockDiagSigma(InvBlockDiagSigma, r2Mat);
      
      //
      arma::mat z1Mat, r1Mat;
      double beta1 = 0;
      arma::mat pMat = z2Mat;
      arma::mat ApMat = getSigmaxMat(pMat, t_excludechr);
      double alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
      xMat = xMat + alpha * pMat;
      r1Mat = r2Mat;
      z1Mat = z2Mat;
      r2Mat = r1Mat - alpha * ApMat;
      
      meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
      
      while (meanL2 > m_tolPCG && iter < m_maxiterPCG){
        iter++;
        
        //  z2Mat = minvMat % r2Mat;
        z2Mat = solverBlockDiagSigma(InvBlockDiagSigma, r2Mat);
        //
        beta1 = getInnerProd(z2Mat, r2Mat) / getInnerProd(z1Mat, r1Mat);
        pMat = z2Mat + beta1 * pMat;
        ApMat = getSigmaxMat(pMat, t_excludechr);
        alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
        xMat = xMat + alpha * pMat;
        r1Mat = r2Mat;
        z1Mat = z2Mat;
        r2Mat = r1Mat - alpha * ApMat;
        meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
      }
    }
    
    t_xVec = convert1(xMat, m_n, m_J);
    if (iter >= m_maxiterPCG){
      cout << "pcg did not converge. You may increase maxiter number." << endl;
    }
    if(m_showInfo)
      cout << "iter from getPCG1ofSigmaAndVector " << iter << endl; 
  // }
}

double POLMMClass::getVarP(arma::vec t_adjGVec,
                           string t_excludechr)
{
  arma::vec adjGVecLong = tZMat(t_adjGVec);
  arma::vec iSigmaGVec(m_n * (m_J-1), arma::fill::zeros);
  getPCGofSigmaAndVector(adjGVecLong, iSigmaGVec, t_excludechr);
  double VarP = as_scalar(adjGVecLong.t() * (iSigmaGVec - m_iSigmaX_XSigmaX * (m_CovaMat.t() * iSigmaGVec)));
  return(VarP);
}

arma::mat POLMMClass::solverBlockDiagSigma(arma::cube& InvBlockDiagSigma,   // (J-1) x (J-1) x n
                                           arma::mat& xMat)                 // n x (J-1)
{
  arma::mat outMat(m_n, m_J-1);
  for(int i = 0; i < m_n; i++){
    outMat.row(i) = xMat.row(i) * InvBlockDiagSigma.slice(i); // could invert matrix?? be careful!
  }
  return(outMat);
}

// yMat = Sigma %*% xMat
arma::mat POLMMClass::getSigmaxMat(arma::mat t_xMat,   // matrix: n x (J-1) 
                                   string t_excludechr)
{
  arma::mat iR_xMat = m_iRMat % t_xMat;
  arma::mat iPsi_iR_xMat = getiPsixMat(iR_xMat);
  arma::mat yMat = m_iRMat % iPsi_iR_xMat;
  if(m_tau == 0){}
  else{
    arma::vec tZ_xMat = getRowSums(t_xMat);  // rowSums(xMat): n x 1
    arma::vec V_tZ_xMat = getKinbVecPOLMM(tZ_xMat, t_excludechr);
    yMat.each_col() += m_tau * V_tZ_xMat;
  }
  return(yMat);
}

// outMat = iPsiMat %*% xMat, iPsiMat is determined by muMat
arma::mat POLMMClass::getiPsixMat(arma::mat t_xMat)   // matrix with dim of n x (J-1)
{
  arma::mat iPsi_xMat(m_n, m_J-1);
  for(int i = 0; i < m_n; i ++){   // loop for samples
    double sumx = sum(t_xMat.row(i));
    for(int j = 0; j < m_J-1; j ++){
      iPsi_xMat(i,j) = sumx / m_muMat(i, m_J-1) + t_xMat(i,j) / m_muMat(i,j);
    }
  }
  return(iPsi_xMat);
}

// used in getPCGofSigmaAndVector()
arma::cube POLMMClass::getInvBlockDiagSigma()
{
  // get diagonal elements of GRM
  arma::vec DiagGRM;
  if(m_flagSparseGRM){
    DiagGRM = m_tau * m_SparseGRM_all.diag();
  }else{
    arma::vec* pDiagStdGeno = m_ptrDenseGRMObj->getDiagStdGeno();
    DiagGRM = m_tau * (*pDiagStdGeno);   // n x 1
  }
  
  double temp;
  //
  arma::cube InvBlockDiagSigma(m_J-1, m_J-1, m_n, arma::fill::zeros);
  for(int i = 0; i < m_n; i++){
    for(int j2 = 0; j2 < m_J-1; j2++){
      for(int j1 = 0; j1 < m_J-1; j1++){
        temp = m_iRMat(i,j2) * (1 / m_muMat(i, m_J-1)) * m_iRMat(i,j1) + DiagGRM(i);
        if(j2 == j1){
          temp += m_iRMat(i,j2) * (1 / m_muMat(i,j2)) * m_iRMat(i,j1); 
        }
        InvBlockDiagSigma(j2, j1, i) = temp;
      }
    }
    InvBlockDiagSigma.slice(i) = inv(InvBlockDiagSigma.slice(i));
  }
  return(InvBlockDiagSigma);
}

// update (m_eta, m_WMat, m_muMat, m_mMat, m_nuMat, m_iRMat, m_YMat) based on (m_beta, m_bVec, m_eps)
void POLMMClass::updateMats()
{
  // update (m_eta)
  m_eta = m_Cova * m_beta + m_bVec;
  
  // update (m_WMat, m_muMat, m_mMat, m_nuMat)
  double tmpExp, tmpnu0, tmpnu1;
  for(int i = 0; i < m_n; i ++){  // loop for samples
    tmpnu0 = 0;  // eps_0 = -Inf
    for(int j = 0; j < m_J-1; j ++){  // loop from eps_1 to eps_{J-1}
      tmpExp = exp(m_eps(j) - m_eta(i));
      tmpnu1 = tmpExp / (1 + tmpExp);
      m_muMat(i,j) = tmpnu1 - tmpnu0;
      m_WMat(i,j) = tmpnu1 * (1 - tmpnu1);
      m_mMat(i,j) = tmpnu1 + tmpnu0 - 1;
      m_nuMat(i,j) = tmpnu1;
      tmpnu0 = tmpnu1;
    }
    int j = m_J-1;      // eps_J = Inf
    tmpnu1 = 1;
    m_muMat(i,j) = tmpnu1 - tmpnu0;
    m_WMat(i,j) = tmpnu1 * (1 - tmpnu1);
    m_mMat(i,j) = tmpnu1 + tmpnu0 - 1;
    m_nuMat(i,j) = tmpnu1;
  }
  
  // update (iRMat)
  for(int i = 0; i < m_n; i ++){
    for(int j = 0; j < m_J-1; j ++){
      m_iRMat(i,j) = 1 / (m_mMat(i,j) - m_mMat(i, m_J-1));
    }
  }
  
  // update (YMat)
  arma::mat xMat = m_yMat.cols(0, m_J-2) - m_muMat.cols(0, m_J-2);
  arma::mat iPsi_xMat = getiPsixMat(xMat);
  for(int i = 0; i < m_n; i++){  // loop for samples
    for(int j = 0; j < m_J-1; j++)
      m_YMat(i,j) = m_eta(i) + (m_iRMat(i,j) * iPsi_xMat(i,j));
  }
}

// set up m_TraceRandMat (TRM) and m_V_TRM, only used once at setPOLMMObj()
void POLMMClass::getTraceRandMat()
{
  arma::vec t1  = getTime();
  for(int itrace = 0; itrace < m_tracenrun; itrace++)
  {
    arma::vec uVec = nb(m_n * (m_J-1));
    uVec = uVec * 2 - 1;
    m_TraceRandMat.col(itrace) = uVec;
    arma::vec ZuVec = ZMat(uVec);
    // m_V_TRM.col(itrace) = tZMat(getKinbVecPOLMM(ZuVec, "none"));
    arma::vec tempVec = getKinbVecPOLMM(ZuVec, "none");
    m_V_TRM.col(itrace) = tZMat(tempVec);
  }

  arma::vec t2  = getTime();
  std::string info = "calculate " + std::to_string(m_tracenrun) + " genKinbVec()";
  printTime(t1, t2, info);
}

arma::vec POLMMClass::getKinbVecPOLMM(arma::vec t_bVec, string t_excludeChr)
{
  arma::vec KinbVec;
  if(m_flagSparseGRM){
    arma::sp_mat temp = m_SparseGRM[t_excludeChr];
    KinbVec = temp * t_bVec;
  }else{
    KinbVec = getKinbVec(t_bVec, m_ptrDenseGRMObj, t_excludeChr, m_grainSize);
  }
  Rcpp::checkUserInterrupt();
  return KinbVec;
}

void POLMMClass::setPOLMMInner(arma::mat t_Cova,
                               arma::Col<int> t_yVec,     // should be from 1 to J
                               arma::vec t_beta,
                               arma::vec t_bVec,
                               arma::vec t_eps,           // 
                               double t_tau)
{
  m_n = t_Cova.n_rows;
  m_p = t_Cova.n_cols;
  m_J = max(t_yVec);
  
  m_CovaMat = getCovaMat(t_Cova, m_n, m_J, m_p);
  m_yVec = t_yVec;
  m_yMat = getyMat();
  
  m_Cova = t_Cova;
  m_beta = t_beta;
  m_bVec = t_bVec;
  m_eps = t_eps;
  m_tau = t_tau;
  
  m_initParList = List::create(Named("beta") = m_beta,
                             Named("bVec") = m_bVec,
                             Named("eps") = m_eps,
                             Named("tau") = m_tau);
  
  setArray();
  set_seed(m_seed);
}

// yMat: matrix with dim of n x J
arma::mat POLMMClass::getyMat()
{
  arma::mat yMat(m_n, m_J, arma::fill::zeros);
  for(int i = 0; i < m_n; i++)
    yMat(i, m_yVec(i)-1) = 1;
  return(yMat);
}

void POLMMClass::setArray()
{
  m_WMat.zeros(m_n, m_J);    
  m_muMat.zeros(m_n, m_J);   
  m_mMat.zeros(m_n, m_J);
  m_nuMat.zeros(m_n, m_J);
  m_iRMat.zeros(m_n, m_J-1);
  m_YMat.zeros(m_n, m_J-1);
  m_iSigma_CovaMat.zeros(m_n * (m_J-1), m_p);
  //
  m_eta.zeros(m_n);
  m_iSigma_YVec.zeros(m_n * (m_J-1));
  m_iSigma_VPYVec.zeros(m_n * (m_J-1));
  //
  m_TraceRandMat.zeros(m_n * (m_J-1), m_tracenrun);
  m_V_TRM.zeros(m_n * (m_J-1), m_tracenrun);
  m_iSigma_V_TRM.zeros(m_n * (m_J-1), m_tracenrun);
}

void POLMMClass::setControlList(Rcpp::List t_controlList)
{
  m_controlList = t_controlList;
  m_memoryChunk = t_controlList["memoryChunk"];
  m_minMafGRM = t_controlList["minMafGRM"];
  m_maxMissingGRM = t_controlList["maxMissingGRM"];
  m_maxiter = t_controlList["maxiter"];
  m_maxiterPCG = t_controlList["maxiterPCG"]; 
  m_maxiterEps = t_controlList["maxiterEps"];
  m_tolBeta = t_controlList["tolBeta"]; 
  m_tolTau = t_controlList["tolTau"]; 
  m_tolPCG = t_controlList["tolPCG"]; 
  m_tolEps = t_controlList["tolEps"];
  m_tracenrun = t_controlList["tracenrun"]; 
  m_seed = t_controlList["seed"];
  m_minMafVarRatio = t_controlList["minMafVarRatio"];
  m_maxMissingVarRatio = t_controlList["maxMissingVarRatio"];
  m_nSNPsVarRatio = t_controlList["nSNPsVarRatio"];
  m_CVcutoff = t_controlList["CVcutoff"];
  m_LOCO = t_controlList["LOCO"];
  m_grainSize = t_controlList["grainSize"];
  m_showInfo = t_controlList["showInfo"];
}

arma::sp_mat POLMMClass::setZMat_sp()
{
  arma::umat locations(2, m_n * (m_J-1));
  arma::urowvec l1 = arma::linspace<arma::urowvec>(0, m_n * (m_J-1) - 1, m_n * (m_J-1));
  arma::urowvec l2 = l1 / (m_J-1);
  locations.row(0) = l1;
  locations.row(1) = l2;
  arma::vec values(m_n * (m_J-1));
  values.fill(1);
  arma::sp_mat ZMat_sp(locations, values);
  return ZMat_sp;
}

// sum up each (J-1) elements: n(J-1) x 1 -> n x 1
arma::vec POLMMClass::ZMat(arma::vec t_xVec)
{
  arma::vec y1Vec(m_n, arma::fill::zeros);
  int index = 0;
  for(int i = 0; i < m_n; i ++){
    for(int j = 0; j < m_J-1; j ++){
      y1Vec(i) += t_xVec(index);
      index++;
    }
  }
  return(y1Vec);
}

// duplicate each element for (J-1) times: n x 1 -> n(J-1) x 1 
arma::vec POLMMClass::tZMat(arma::vec t_xVec)
{
  arma::vec y1Vec(m_n * (m_J-1));
  int index = 0;
  for(int i = 0; i < m_n; i ++){
    for(int j = 0; j < m_J-1; j ++){
      y1Vec(index) = t_xVec(i);
      index++;
    }
  }
  return(y1Vec);
}

void POLMMClass::updateEps()
{
  for(int iter = 0; iter < m_maxiterEps; iter ++){
    
    arma::vec eps0 = m_eps;
    
    updateEpsOneStep();
    updateMats();
    
    double diffeps = max(abs(m_eps - eps0)/(abs(m_eps) + abs(eps0) + m_tolEps));
    
    if(diffeps < m_tolEps){
      
      cout << "UpdateEps iter: " << iter << endl;
      cout << "eps: " << endl << m_eps << endl;
      break;
    }
  }
}

// need update in case that eps(k+1) < eps(k)
void POLMMClass::updateEpsOneStep()
{
  // the first eps is fixed at 0
  arma::vec d1eps(m_J-2, arma::fill::zeros);
  arma::mat d2eps(m_J-2, m_J-2, arma::fill::zeros);
  double temp1, temp2, temp3;
  
  for(int k = 1; k < m_J-1; k++){
    for(int i = 0; i < m_n; i++){
      temp1 = m_yMat(i, k) / m_muMat(i, k) - m_yMat(i, k+1) / m_muMat(i, k+1);
      temp2 = - m_yMat(i, k) / m_muMat(i, k) / m_muMat(i, k) - m_yMat(i, k+1) / m_muMat(i, k+1) / m_muMat(i, k+1);
      d1eps(k - 1) += m_WMat(i, k) * temp1;
      d2eps(k - 1, k - 1) += m_WMat(i, k) * (1 - 2 * m_nuMat(i, k)) * temp1 + m_WMat(i, k) * m_WMat(i, k) * temp2;
      if(k < m_J-2){
        temp3 = m_WMat(i,k) * m_WMat(i, k+1) * m_yMat(i, k+1) / m_muMat(i, k+1) / m_muMat(i, k+1);
        d2eps(k-1, k) += temp3;
        d2eps(k, k-1) += temp3;
      }
    }
  }
  
  arma::vec deps = -1 * inv(d2eps) * d1eps;
  for(int k = 1; k < m_J-1; k ++){
    m_eps(k) += deps(k-1);
  }
}

Rcpp::List POLMMClass::getPOLMM()
{
  Rcpp::List outList = List::create(Named("N") = m_n,              // number of samples
                                    Named("M") = m_M,              // number of SNPs in Plink file
                                    Named("controlList") = m_controlList,
                                    Named("iter") = m_iter,
                                    Named("eta") = m_eta,          // X %*% beta + bVec
                                    Named("yVec") = m_yVec,        // matrix with dim of n x 1: observation
                                    Named("Cova") = m_Cova,        // matrix with dim of n(J-1) x p: covariates
                                    Named("muMat") = m_muMat,      // matrix with dim of n x J: probability
                                    Named("YMat") = m_YMat,        // matrix with dim of n x (J-1): working variables
                                    Named("beta") = m_beta,        // parameter for covariates
                                    Named("bVec") = m_bVec,        // terms of random effect 
                                    Named("tau") = m_tau,          // variance component
                                    Named("eps") = m_eps,          // cutpoints
                                    Named("LOCOList") = m_LOCOList);         
  
  return(outList);
}

void POLMMClass::closeGenoObj()
{
  if(!m_flagSparseGRM)
    m_ptrDenseGRMObj->closeDenseGRMObj();
}

}