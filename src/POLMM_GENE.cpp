
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <string>

#include "SubFunc.hpp"
#include "POLMM_GENE.hpp"

namespace POLMMGENE{

using namespace Rcpp;
using namespace std;

Rcpp::List POLMMGENEClass::getStatVarS(arma::mat t_GMat)
{
  int n = t_GMat.n_rows;  // number of samples
  int q = t_GMat.n_cols;  // number of markers in the region
  
  arma::vec StatVec(q);
  arma::mat VarSMat(q, q);
  arma::mat adjGMat(n, q);
  arma::mat ZPZ_adjGMat(n, q);
  
  for(int i = 0; i < q; i++){
    
    // get adjusted GVec
    arma::vec GVec = t_GMat.col(i);
    arma::vec adjGVec = getadjGFast(GVec, m_objP["XXR_Psi_RX_new"], m_objP["XR_Psi_R_new"], m_objP["n"], m_objP["J"], m_objP["p"]);
    adjGMat.col(i) = adjGVec;
    
    // get Stat for the marker
    double Stat = getStatFast(adjGVec, m_objP["RymuVec"], m_objP["n"]);
    StatVec(i) = Stat;
    
    // get t(Z) %*% P %*% Z %*% adjGVec
    arma::vec ZPZ_adjGVec = get_ZPZ_adjGVec(adjGVec, m_excludechr);
    ZPZ_adjGMat.col(i) = ZPZ_adjGVec;
  }
  
  for(int i = 0; i < q; i++){
    for(int j = 0; j < i; j++){
      double cov = as_scalar(adjGMat.col(i).t() * ZPZ_adjGMat.col(j));
      VarSMat(i, j) = VarSMat(j,i) = cov;
    }
    double var = as_scalar(adjGMat.col(i).t() * ZPZ_adjGMat.col(i));
    VarSMat(i, i) = var;
  }
  
  Rcpp::List OutList = List::create(Named("StatVec") = StatVec,             
                                    Named("VarSMat") = VarSMat);
  return OutList;
}

void POLMMGENEClass::setPOLMMGENEobj(int t_maxiterPCG,
                                     double t_tolPCG,
                                     arma::mat t_Cova,
                                     arma::Col<int> t_yVec,     // should be from 1 to J
                                     double t_tau,
                                     Rcpp::List t_SparseGRM,    // results of function getKinMatList()
                                     Rcpp::List t_LOCOList,
                                     arma::vec t_eta)
{
  m_flagSparseGRM = true;
  m_maxiterPCG = t_maxiterPCG; 
  m_tolPCG = t_tolPCG;
  m_tau = t_tau;
  
  m_Cova = t_Cova;
  m_yVec = t_yVec;
  
  m_n = m_Cova.n_rows;
  m_p = m_Cova.n_cols;
  m_J = max(m_yVec);
  
  m_CovaMat = getCovaMat(t_Cova, m_n, m_J, m_p);
  m_yMat = getyMat();
  
  m_SparseGRM = t_SparseGRM;
  m_LOCOList = t_LOCOList;
  m_eta = t_eta;
}

void POLMMGENEClass::setPOLMMGENEchr(Rcpp::List t_LOCOList, string t_excludechr)
{
  m_excludechr = t_excludechr;
  arma::sp_mat temp = m_SparseGRM[t_excludechr];
  m_SparseGRM_all = temp;
  
  Rcpp::List LOCO = t_LOCOList[t_excludechr];
  
  arma::mat temp1 = LOCO["muMat"];
  m_muMat = temp1;
  arma::mat temp2 = LOCO["iRMat"];
  m_iRMat = temp2;

  m_objP = getobjP(m_Cova, m_yMat, m_muMat, m_iRMat);
  
  // the below is to make m_iSigmaX_XSigmaX
  arma::mat iSigma_CovaMat(m_n * (m_J-1), m_p);
  getPCGofSigmaAndCovaMat(m_CovaMat, iSigma_CovaMat, t_excludechr);
  
  arma::mat xMat = m_yMat.cols(0, m_J-2) - m_muMat.cols(0, m_J-2);
  arma::mat iPsi_xMat = getiPsixMat(xMat);
  arma::mat YMat(m_n, m_J-1, arma::fill::zeros);
  for(int i = 0; i < m_n; i++){  // loop for samples
    for(int j = 0; j < m_J-1; j++)
      YMat(i,j) = m_eta(i) + (m_iRMat(i,j) * iPsi_xMat(i,j));
  }
  
  arma::vec YVec = convert1(YMat, m_n, m_J);
  arma::vec iSigma_YVec(m_n * (m_J-1));
  getPCGofSigmaAndVector(YVec, iSigma_YVec, t_excludechr); 
  
  arma::mat XSigmaX = inv(m_CovaMat.t() * iSigma_CovaMat);
  m_iSigmaX_XSigmaX = iSigma_CovaMat * XSigmaX;
  
}

void POLMMGENEClass::closePOLMMGENEobj()
{
  // add something later
}

// sum up each (J-1) elements: n(J-1) x 1 -> n x 1
arma::vec POLMMGENEClass::ZMat(arma::vec t_xVec)
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
arma::vec POLMMGENEClass::tZMat(arma::vec t_xVec)
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

arma::vec POLMMGENEClass::getKinbVecPOLMM(arma::vec t_bVec, string t_excludeChr)
{
  arma::vec KinbVec;
  if(m_flagSparseGRM){
    arma::sp_mat temp = m_SparseGRM[t_excludeChr];
    KinbVec = temp * t_bVec;
  }else{
    // KinbVec = getKinbVec(t_bVec, m_ptrDenseGRMObj, t_excludeChr, m_grainSize);
    // for POLMM_GENE, we set m_flatSparseGRM = T, this is not used
  }
  Rcpp::checkUserInterrupt();
  return KinbVec;
}

// outMat = iPsiMat %*% xMat, iPsiMat is determined by muMat
arma::mat POLMMGENEClass::getiPsixMat(arma::mat t_xMat)   // matrix with dim of n x (J-1)
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

// yMat = Sigma %*% xMat
arma::mat POLMMGENEClass::getSigmaxMat(arma::mat t_xMat,   // matrix: n x (J-1) 
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

// used in getPCGofSigmaAndVector()
arma::cube POLMMGENEClass::getInvBlockDiagSigma()
{
  // get diagonal elements of GRM
  arma::vec DiagGRM;
  
  if(m_flagSparseGRM){
    DiagGRM = m_tau * m_SparseGRM_all.diag();
  }else{
    // arma::vec* pDiagStdGeno = m_ptrDenseGRMObj->getDiagStdGeno();
    // DiagGRM = m_tau * (*pDiagStdGeno);   // n x 1
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

arma::mat POLMMGENEClass::solverBlockDiagSigma(arma::cube& InvBlockDiagSigma,   // (J-1) x (J-1) x n
                                               arma::mat& xMat)                 // n x (J-1)
{
  arma::mat outMat(m_n, m_J-1);
  for(int i = 0; i < m_n; i++){
    outMat.row(i) = xMat.row(i) * InvBlockDiagSigma.slice(i); // could invert matrix?? be careful!
  }
  return(outMat);
}

// use PCG to calculate iSigma_xMat = Sigma^-1 %*% xMat
void POLMMGENEClass::getPCGofSigmaAndCovaMat(arma::mat t_xMat,              // matrix with dim of n(J-1) x p
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

// use PCG to calculate xVec = Sigma^-1 %*% yVec
void POLMMGENEClass::getPCGofSigmaAndVector(arma::vec t_y1Vec,    // vector with length of n(J-1)
                                            arma::vec& t_xVec,    // vector with length of n(J-1)
                                            string t_excludechr)
{
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

arma::vec POLMMGENEClass::get_ZPZ_adjGVec(arma::vec t_adjGVec,
                                          string t_excludechr)
{
  arma::vec adjGVecLong = tZMat(t_adjGVec);  // that is Z %*% adjGVec
  arma::vec iSigmaGVec(m_n * (m_J-1), arma::fill::zeros);
  getPCGofSigmaAndVector(adjGVecLong, iSigmaGVec, t_excludechr);  // iSigmaGVec = Sigma^-1 %*% Z %*% adjGVec
  
  arma::vec PZ_adjGVec = iSigmaGVec - m_iSigmaX_XSigmaX * (m_CovaMat.t() * iSigmaGVec);
  arma::vec ZPZ_adjGVec = ZMat(PZ_adjGVec);
  
  return(ZPZ_adjGVec);
}

// yMat: matrix with dim of n x J
arma::mat POLMMGENEClass::getyMat()
{
  arma::mat yMat(m_n, m_J, arma::fill::zeros);
  for(int i = 0; i < m_n; i++)
    yMat(i, m_yVec(i)-1) = 1;
  return(yMat);
}

}

