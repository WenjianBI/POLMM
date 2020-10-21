
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <string>

#include "SubFunc.hpp"
#include "POLMM_GENE.hpp"
#include "ER_SPA.hpp"

namespace POLMMGENE{

using namespace Rcpp;
using namespace std;

Rcpp::List POLMMGENEClass::getStatVarS(arma::mat t_GMat, 
                                       double t_NonZero_cutoff,
                                       double t_StdStat_cutoff)
{
  int n = t_GMat.n_rows;       // number of samples
  int q = t_GMat.n_cols;       // number of markers in the region
  
  arma::vec StatVec(q);        // vector of statistics
  arma::mat VarSMat(q, q);     // variance matrix (after adjusting for relatedness)
  arma::mat adjGMat(n, q);     // adjusted genotype vector
  arma::mat ZPZ_adjGMat(n, q); // t(Z) %*% P %*% Z %*% adjGMat
  arma::vec NonZeroVec(q);     // vector of NonZero number [for method selection]
  arma::vec StdStatVec(q);     // vector of standardized statistics [for method selection]
  arma::vec VarWVec(q);        // vector of variance (without adjusting for relatedness) [for fastSPA]
  arma::vec Ratio0Vec(q);      // vector of Ratio0 [for fastSPA]
  arma::uvec idxERVec(q, arma::fill::zeros);       // 0 or 1: vector of indicator for markers to ER 
  arma::uvec idxSPAVec(q, arma::fill::zeros);      // 0 or 1: vector of indicator for markers to SPA [for fastSPA]
  
  // loop for all markers
  for(int i = 0; i < q; i++){
    
    // get adjusted GVec each marker
    arma::vec GVec = t_GMat.col(i);
    arma::vec adjGVec = getadjGFast(GVec, m_objP["XXR_Psi_RX_new"], m_objP["XR_Psi_R_new"], 
                                    m_objP["n"], m_objP["p"]);
    adjGMat.col(i) = adjGVec;
    
    // get Stat for each marker
    double Stat = getStatFast(adjGVec, m_objP["RymuVec"]);
    StatVec(i) = Stat;
    
    // get t(Z) %*% P %*% Z %*% adjGVec for each marker
    arma::vec ZPZ_adjGVec = get_ZPZ_adjGVec(adjGVec, m_excludechr);
    ZPZ_adjGMat.col(i) = ZPZ_adjGVec;
    double VarS = as_scalar(adjGVec.t() * ZPZ_adjGVec);
    VarSMat(i, i) = VarS;
    
    // use cutoff to determine normal approximation, SPA, and ER
    arma::uvec idxNonZero = arma::find(GVec > 0.1); // we have changed (GVec < 0.2) to 0
    int NonZero = idxNonZero.size(); 
    double StdStat = std::abs(Stat) / sqrt(VarS);
    NonZeroVec(i) = NonZero;
    StdStatVec(i) = StdStat;
    
    // Efficient Resampling (ER)
    // if(NonZero <= t_NonZero_cutoff){  // updated on 10-15-2020
    if(NonZero <= t_NonZero_cutoff && StdStat > t_StdStat_cutoff){
      idxERVec(i) = 1;
    }
    
    // Saddlepoint approximation (SPA)
    if(NonZero > t_NonZero_cutoff && StdStat > t_StdStat_cutoff){
      idxSPAVec(i) = 1;
      // check getVarWFast()
      double VarW = 0;
      double VarW1 = 0;
      for(int j = 0; j < n; j++){
        double temp = m_RPsiRVec(j) * adjGVec(j) * adjGVec(j);
        VarW += temp;
        if(GVec(j) > 0.1)
          VarW1 += temp;
      }
      double VarW0 = VarW - VarW1;
      double Ratio0 = VarW0 / VarW;
      VarWVec(i) = VarW;
      Ratio0Vec(i) = Ratio0;
    }
  }
  
  for(int i = 0; i < q; i++){
    for(int j = 0; j < i; j++){
      double cov = as_scalar(adjGMat.col(i).t() * ZPZ_adjGMat.col(j));
      VarSMat(i, j) = VarSMat(j,i) = cov;
    }
  }
  
  arma::uvec whichSPA = arma::find(idxSPAVec == 1);
  VarWVec = VarWVec.elem(whichSPA);
  Ratio0Vec = Ratio0Vec.elem(whichSPA);
  adjGMat = adjGMat.cols(whichSPA);
  
  Rcpp::List OutList = List::create(Named("StatVec") = StatVec,
                                    Named("VarSMat") = VarSMat,
                                    Named("NonZeroVec") = NonZeroVec,
                                    Named("StdStatVec") = StdStatVec,
                                    // the following 2 elements are to determine which markers should use fastSPA or ER
                                    Named("idxERVec") = idxERVec,
                                    Named("idxSPAVec") = idxSPAVec,
                                    // the following 5 elements are only for fastSPA, determined by idxSPAVec
                                    Named("muMat") = m_muMat,
                                    Named("iRMat") = m_iRMat,
                                    Named("VarWVec") = VarWVec,           
                                    Named("Ratio0Vec") = Ratio0Vec,
                                    Named("adjGMat") = adjGMat);
  return OutList;
}


POLMMGENEClass::POLMMGENEClass(int t_maxiterPCG,
                               double t_tolPCG,
                               arma::mat t_Cova,
                               arma::uvec t_yVec,         // should be from 1 to J
                               double t_tau,
                               Rcpp::List t_SparseGRM,    // results of function getKinMatList()
                               Rcpp::List t_LOCOList,
                               arma::vec t_eta,
                               int t_nMaxNonZero)
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
  
  m_SeqMat = makeSeqMat(t_nMaxNonZero, m_J);
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
  arma::vec temp3 = m_objP["RPsiR"];
  m_RPsiRVec = temp3;
  
  // the below is to make m_iSigmaX_XSigmaX
  arma::mat iSigma_CovaMat(m_n * (m_J-1), m_p);
  getPCGofSigmaAndCovaMat(m_CovaMat, iSigma_CovaMat, t_excludechr);
  
  // arma::mat xMat = m_yMat.cols(0, m_J-2) - m_muMat.cols(0, m_J-2);
  // arma::mat iPsi_xMat = getiPsixMat(xMat);
  // arma::mat YMat(m_n, m_J-1, arma::fill::zeros);
  // for(int i = 0; i < m_n; i++){  // loop for samples
  //   for(int j = 0; j < m_J-1; j++)
  //     YMat(i,j) = m_eta(i) + (m_iRMat(i,j) * iPsi_xMat(i,j));
  // }
  // 
  // arma::vec YVec = convert1(YMat, m_n, m_J);
  // arma::vec iSigma_YVec(m_n * (m_J-1));
  // getPCGofSigmaAndVector(YVec, iSigma_YVec, t_excludechr); 
  
  arma::mat XSigmaX = inv(m_CovaMat.t() * iSigma_CovaMat);
  m_iSigmaX_XSigmaX = iSigma_CovaMat * XSigmaX;
  
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

void POLMMGENEClass::check_ZPZ_adjGVec(arma::vec t_adjGVec)
{
  arma::vec adjGVecLong = tZMat(t_adjGVec);  // that is Z %*% adjGVec
  arma::vec iSigmaGVec(m_n * (m_J-1), arma::fill::zeros);
  getPCGofSigmaAndVector(adjGVecLong, iSigmaGVec, m_excludechr);  // iSigmaGVec = Sigma^-1 %*% Z %*% adjGVec
  
  arma::vec PZ_adjGVec = iSigmaGVec - m_iSigmaX_XSigmaX * (m_CovaMat.t() * iSigmaGVec);
  arma::vec ZPZ_adjGVec = ZMat(PZ_adjGVec);
  
  std::cout << "m_excludechr\t" << m_excludechr << std::endl;
  std::cout << "adjGVecLong.head(10)\t" << adjGVecLong.head(10) << std::endl;
  std::cout << "iSigmaGVec.head(10)\t" << iSigmaGVec.head(10) << std::endl;
  std::cout << "m_CovaMat.head_rows(10)\t" << m_CovaMat.head_rows(10) << std::endl;
  std::cout << "m_iSigmaX_XSigmaX.head_rows(10)\t" << m_iSigmaX_XSigmaX.head_rows(10) << std::endl;
  std::cout << "PZ_adjGVec.head(10)\t" << PZ_adjGVec.head(10) << std::endl;
  std::cout << "ZPZ_adjGVec.head(10)\t" << ZPZ_adjGVec.head(10) << std::endl;
}

double POLMMGENEClass::checkError(){
  double errorCHR = m_iSigmaX_XSigmaX(0,0);
  return errorCHR;
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


double POLMMGENEClass::getPvalERinClass(arma::vec t_GVec)
{
  // here, we use 0.1 to be consistent with the imputation step, imputation and allele flip should be finished in previous steps
  arma::uvec idxVec = arma::find(t_GVec > 0.1);  
  int n1 = idxVec.size();
  arma::umat SeqMat = updateSeqMat(m_SeqMat, n1, m_J);
  
  // only use subjects in idxVec
  arma::vec GVec = t_GVec.elem(idxVec);
  arma::mat muMat = m_muMat.rows(idxVec);
  arma::uvec yVec = m_yVec.elem(idxVec); // from 1 to J
  arma::mat iRMat = m_iRMat.rows(idxVec);
  
  double pvalER = getPvalER(yVec, GVec, muMat, iRMat);
  return pvalER;
}


}

