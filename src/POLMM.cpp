
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "POLMM.hpp"
#include "util.hpp"


namespace POLMM {

POLMMClass::POLMMClass(arma::mat t_muMat,
                       arma::mat t_iRMat,
                       arma::mat t_Cova,
                       arma::vec t_yVec,
                       arma::sp_mat t_SparseGRM,
                       double t_tau,
                       bool t_printPCGInfo,
                       double t_tolPCG,
                       int t_maxiterPCG)
{
  m_muMat = t_muMat;
  m_iRMat = t_iRMat;
  m_Cova = t_Cova;
  
  m_n = m_muMat.n_rows;
  m_J = m_muMat.n_cols;
  m_p = m_Cova.n_cols;
  
  m_CovaMat = getCovaMat(m_Cova, m_n, m_J, m_p);       // n(J-1) x p
  
  m_SparseGRM = t_SparseGRM;
  m_tau = t_tau;
  m_printPCGInfo = t_printPCGInfo;
  m_tolPCG = t_tolPCG;
  m_maxiterPCG = t_maxiterPCG;
  
  m_InvBlockDiagSigma = getInvBlockDiagSigma();
  
  // output for Step 2
  arma::mat XR_Psi_R(m_p, m_n * (m_J-1));                // p x n(J-1)
  for(int k = 0; k < m_p; k++){
    arma::mat xMat = Vec2Mat(m_CovaMat.col(k), m_n, m_J);
    arma::vec temp = Mat2Vec(getPsixMat(xMat / m_iRMat) / m_iRMat, m_n, m_J);
    XR_Psi_R.row(k) = temp.t();
  }
  
  m_XXR_Psi_RX = m_Cova * inv(XR_Psi_R * m_CovaMat);               // (n x p) * (p x p) = n x p
  
  // sum each (J-1) rows to 1 row: p x n(J-1) -> p x n
  m_XR_Psi_R = sumCols(XR_Psi_R, m_J);      // p x n
  
  arma::mat yMat(m_n, m_J, arma::fill::zeros);
  for(int i = 0; i < m_n; i++)
    yMat(i, t_yVec(i)-1) = 1;
  
  arma::mat ymuMat = yMat - m_muMat;                    // n x J
  arma::mat RymuMat = ymuMat.cols(0, m_J-2) / t_iRMat;    // n x (J-1): R %*% (y - mu)
  m_RymuVec = sumCols(RymuMat, m_J);            // n x 1
  
  arma::mat iSigma_CovaMat(m_n * (m_J-1), m_p);
  getPCGofSigmaAndCovaMat(m_CovaMat, iSigma_CovaMat);
  
  arma::mat XSigmaX = inv(m_CovaMat.t() * iSigma_CovaMat);
  m_iSigmaX_XSigmaX = iSigma_CovaMat * XSigmaX;
}

arma::vec POLMMClass::getadjGFast(arma::vec t_GVec)
{
  // To increase computational efficiency when lots of GVec elements are 0
  arma::vec XR_Psi_RG(m_p, arma::fill::zeros);
  for(int i = 0; i < m_n; i++){
    if(t_GVec(i) != 0){
      XR_Psi_RG += m_XR_Psi_R.col(i) * t_GVec(i);
    }
  }
  arma::vec adjGVec = t_GVec - m_XXR_Psi_RX * XR_Psi_RG;
  return adjGVec;
}

double POLMMClass::getStatFast(arma::vec t_adjGVec)         // n x 1
{
  double Stat = 0;
  for(int i = 0; i < m_n; i++){
    if(t_adjGVec(i) != 0){
      Stat += t_adjGVec(i) * m_RymuVec(i);
    }
  }
  return Stat;
}

arma::vec POLMMClass::get_ZPZ_adjGVec(arma::vec t_adjGVec)
{
  arma::vec adjGVecLong = Vec2LongVec(t_adjGVec, m_n, m_J);  // that is Z %*% adjGVec
  
  arma::vec iSigmaGVec(m_n * (m_J-1), arma::fill::zeros);
  
  getPCGofSigmaAndVector(adjGVecLong, iSigmaGVec);  // iSigmaGVec = Sigma^-1 %*% Z %*% adjGVec
  
  arma::vec PZ_adjGVec = iSigmaGVec - m_iSigmaX_XSigmaX * (m_CovaMat.t() * iSigmaGVec);
  arma::vec ZPZ_adjGVec = LongVec2Vec(PZ_adjGVec, m_n, m_J);
  
  return ZPZ_adjGVec;
}

// use PCG to calculate iSigma_xMat = Sigma^-1 %*% xMat
void POLMMClass::getPCGofSigmaAndCovaMat(arma::mat t_xMat,              // matrix with dim of n(J-1) x p
                                         arma::mat& t_iSigma_xMat)      // matrix with dim of n(J-1) x p
{
  int p1 = t_xMat.n_cols;
  for(int i = 0; i < p1; i++){
    arma::vec y1Vec = t_xMat.col(i);
    arma::vec iSigma_y1Vec = t_iSigma_xMat.col(i);
    getPCGofSigmaAndVector(y1Vec, iSigma_y1Vec);
    t_iSigma_xMat.col(i) = iSigma_y1Vec;
  }
}

// use PCG to calculate xVec = Sigma^-1 %*% yVec
void POLMMClass::getPCGofSigmaAndVector(arma::vec t_y1Vec,    // vector with length of n(J-1)
                                        arma::vec& t_xVec)    // vector with length of n(J-1)
{
  arma::mat xMat = Vec2Mat(t_xVec, m_n, m_J);
  arma::mat y1Mat = Vec2Mat(t_y1Vec, m_n, m_J);
  
  // r2Vec and z2Vec are for the current step; r1Vec and z1Vec are for the previous step
  int iter = 0;
  arma::mat r2Mat = y1Mat - getSigmaxMat(xMat);  // n x (J-1): r0 = y1Mat- Sigma %*% xMat
  double meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
  if(meanL2 <= m_tolPCG){
    // do nothing, xMat is already close to (Sigma)^-1 %*% y1Mat
  }else{
    
    iter++;
    arma::mat z2Mat = solverBlockDiagSigma(r2Mat);
    //
    arma::mat z1Mat, r1Mat;
    double beta1 = 0;
    arma::mat pMat = z2Mat;
    arma::mat ApMat = getSigmaxMat(pMat);
    double alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
    xMat = xMat + alpha * pMat;
    r1Mat = r2Mat;
    z1Mat = z2Mat;
    r2Mat = r1Mat - alpha * ApMat;
    
    meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
    
    while (meanL2 > m_tolPCG && iter < m_maxiterPCG){
      
      iter++;
      //  z2Mat = minvMat % r2Mat;
      z2Mat = solverBlockDiagSigma(r2Mat);
      //
      beta1 = getInnerProd(z2Mat, r2Mat) / getInnerProd(z1Mat, r1Mat);
      pMat = z2Mat + beta1 * pMat;
      ApMat = getSigmaxMat(pMat);
      alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
      xMat = xMat + alpha * pMat;
      r1Mat = r2Mat;
      z1Mat = z2Mat;
      r2Mat = r1Mat - alpha * ApMat;
      meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
    }
  }
  
  t_xVec = Mat2Vec(xMat, m_n, m_J);
  if (iter >= m_maxiterPCG){
    std::cout << "pcg did not converge. You may increase maxiter number." << std::endl;
  }
  if(m_printPCGInfo)
    std::cout << "iter from getPCG1ofSigmaAndVector " << iter << std::endl; 
  // }
}

// yMat = Sigma %*% xMat
arma::mat POLMMClass::getSigmaxMat(arma::mat& t_xMat)   // matrix: n x (J-1) 
{
  arma::mat iR_xMat = m_iRMat % t_xMat;
  arma::mat iPsi_iR_xMat = getiPsixMat(iR_xMat);
  arma::mat yMat = m_iRMat % iPsi_iR_xMat;
  if(m_tau == 0){}
  else{
    arma::vec tZ_xMat = arma::sum(t_xMat, 1);  // rowSums(xMat): n x 1
    arma::vec V_tZ_xMat = m_SparseGRM * tZ_xMat;
    yMat.each_col() += m_tau * V_tZ_xMat;
  }
  return yMat;
}

// outMat = iPsiMat %*% xMat, iPsiMat is determined by muMat
arma::mat POLMMClass::getiPsixMat(arma::mat t_xMat)   // matrix with dim of n x (J-1)
{
  arma::mat iPsi_xMat(m_n, m_J-1);
  for(int i = 0; i < m_n; i ++){   // loop for samples
    double sumx = arma::sum(t_xMat.row(i));
    for(int j = 0; j < m_J-1; j++){
      iPsi_xMat(i,j) = sumx / m_muMat(i, m_J-1) + t_xMat(i,j) / m_muMat(i,j);
    }
  }
  return iPsi_xMat;
}

// outMat = PsiMat %*% xMat, PsiMat is determined by muMat
arma::mat POLMMClass::getPsixMat(arma::mat t_xMat)   // matrix: n x (J-1)
{
  arma::mat Psi_xMat(m_n, m_J-1);
  // loop for samples
  for(int i = 0; i < m_n; i++){
    arma::rowvec muVec(m_J-1);
    for(int j = 0; j < m_J-1; j++){
      Psi_xMat(i,j) = m_muMat(i,j) * t_xMat(i,j);
      muVec(j) = m_muMat(i,j);
    }
    double sum_mu_x = sum(Psi_xMat.row(i));
    Psi_xMat.row(i) -= muVec * sum_mu_x; 
  }
  return Psi_xMat;
}

// used in getPCGofSigmaAndVector()
arma::cube POLMMClass::getInvBlockDiagSigma()
{
  // get diagonal elements of GRM
  arma::vec DiagGRM;
  DiagGRM = m_tau * m_SparseGRM.diag();
  
  arma::cube InvBlockDiagSigma(m_J-1, m_J-1, m_n, arma::fill::zeros);
  for(int i = 0; i < m_n; i++){
    for(int j2 = 0; j2 < m_J-1; j2++){
      for(int j1 = 0; j1 < m_J-1; j1++){
        double temp = m_iRMat(i,j2) * (1 / m_muMat(i, m_J-1)) * m_iRMat(i,j1) + DiagGRM(i);
        if(j2 == j1){
          temp += m_iRMat(i,j2) * (1 / m_muMat(i,j2)) * m_iRMat(i,j1); 
        }
        InvBlockDiagSigma(j2, j1, i) = temp;
      }
    }
    InvBlockDiagSigma.slice(i) = inv(InvBlockDiagSigma.slice(i));
  }
  return InvBlockDiagSigma;
}

arma::mat POLMMClass::solverBlockDiagSigma(arma::mat& t_xMat)     // n x (J-1)
{
  arma::mat outMat(m_n, m_J-1);
  for(int i = 0; i < m_n; i++){
    outMat.row(i) = t_xMat.row(i) * m_InvBlockDiagSigma.slice(i); // could invert matrix?? be careful!
  }
  return outMat;
}

}

// make a global variable for future usage
// static POLMM::POLMMClass* ptr_gPOLMMobj = NULL;

// ptr_gPOLMMobj = new POLMM::POLMMClass(t_muMat,
//                                       t_iRMat,
//                                       t_Cova,
//                                       t_yVec,
//                                       t_SparseGRM,
//                                       t_tau,
//                                       t_printPCGInfo,
//                                       t_tolPCG,
//                                       t_maxiterPCG);

// // [[Rcpp::export]]
// void setPOLMMobjInR(arma::mat t_muMat,
//                     arma::mat t_iRMat,
//                     arma::mat t_Cova,
//                     arma::vec t_yVec,
//                     Rcpp::List t_SPmatR,    // output of makeSPmatR()
//                     double t_tau,
//                     bool t_printPCGInfo,
//                     double t_tolPCG,
//                     int t_maxiterPCG)
// {
//   arma::umat locations = t_SPmatR["locations"];
//   arma::vec values = t_SPmatR["values"];
//   arma::sp_mat SparseGRM = arma::sp_mat(locations, values);
//   ptr_gPOLMMobj = new POLMM::POLMMClass(t_muMat,
//                                         t_iRMat,
//                                         t_Cova,
//                                         t_yVec,
//                                         SparseGRM,
//                                         t_tau,
//                                         t_printPCGInfo,
//                                         t_tolPCG,
//                                         t_maxiterPCG);
// }


