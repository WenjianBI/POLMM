
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
  setRPsiR();
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

Rcpp::List POLMMClass::MAIN_SPA(double t_Stat,
                                arma::vec t_adjGVec,
                                arma::vec t_K1roots,
                                double t_VarP,
                                double t_VarW,
                                double t_Ratio0,
                                arma::uvec t_posG1)
{
  Rcpp::List resSPA = fastSaddle_Prob(t_Stat, t_VarP, t_VarW, t_Ratio0, t_K1roots,
                                      t_adjGVec.elem(t_posG1), m_muMat.rows(t_posG1), m_iRMat.rows(t_posG1));
  return resSPA;
}

void POLMMClass::setRPsiR()   
{
  // arma::cube RPsiR(J-1, J-1, n, arma::fill::zeros);
  arma::vec RPsiRVec(m_n, arma::fill::zeros);
  arma::mat muRMat = m_muMat.cols(0, m_J-2) / m_iRMat;
  for(int i = 0; i < m_n; i++){
    for(int j1 = 0; j1 < m_J-1; j1++){
      // RPsiR(j1,j1,i) += muRMat(i,j1) / iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      RPsiRVec(i) += muRMat(i,j1) / m_iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      for(int j2 = j1+1; j2 < m_J-1; j2++){
        // RPsiR(j1,j2,i) -= muRMat(i,j1) * muRMat(i,j2);
        RPsiRVec(i) -= 2* muRMat(i,j1) * muRMat(i,j2);
      }
    }
  }
  m_RPsiR = RPsiRVec;
}

arma::vec POLMMClass::getVarWVec(arma::vec adjGVec)
{
  arma::vec VarWVec = m_RPsiR % pow(adjGVec, 2);
  return VarWVec;
}

double K0(double t_x,
          arma::mat t_muMat,     // N x (J-1)
          arma::mat t_cMat,      // N x (J-1)
          double t_m1)           // sum(muMat * cMat)
{
  arma::mat temp1Mat = - t_muMat + t_muMat % exp(t_cMat * t_x);
  arma::vec temp1Vec = log(1 + arma::sum(temp1Mat, 1));   // arma::sum(Mat, 1) is rowSums()
  double y = sum(temp1Vec) - t_m1 * t_x;
  return y;
}

arma::vec K12(double t_x,
              arma::mat t_muMat,
              arma::mat t_cMat,
              double t_m1)
{
  arma::mat temp0Mat = t_muMat % exp(t_cMat * t_x);
  arma::mat temp1Mat = - t_muMat + temp0Mat;
  arma::mat temp2Mat = temp0Mat % t_cMat;
  arma::mat temp3Mat = temp2Mat % t_cMat;
  
  arma::vec temp1Vec = 1 + arma::sum(temp1Mat, 1);
  arma::vec temp2Vec = arma::sum(temp2Mat, 1);
  arma::vec temp3Vec = arma::sum(temp3Mat, 1);
  
  arma::vec yVec(2);
  yVec(0) = sum(temp2Vec / temp1Vec) - t_m1;
  // yMat[i,2] = sum((temp3Vec*temp1Vec-temp2Vec^2)/temp1Vec^2, na.rm=TRUE);
  yVec(1) = sum((temp3Vec % temp1Vec - pow(temp2Vec, 2)) / pow(temp1Vec, 2));
  
  return yVec;
}

Rcpp::List fastgetroot_K1(double t_Stat,
                          double t_initX,
                          double t_Ratio0,
                          arma::mat t_muMat,
                          arma::mat t_cMat,
                          double t_m1)
{
  double x = t_initX;
  double K1 = 0;
  double K2 = 0;
  double diffX = arma::datum::inf;
  bool converge = true;
  double tol = 0.001;
  int maxiter = 100;
  int iter = 0;
  
  for(iter = 0; iter < maxiter; iter ++){
    double oldX = x;
    double oldDiffX = diffX;
    double oldK1 = K1;
    
    arma::vec K12Vec = K12(x, t_muMat, t_cMat, t_m1);
    
    K1 = K12Vec(0) - t_Stat + t_Ratio0 * x;
    K2 = K12Vec(1) + t_Ratio0;
    
    diffX = -1 * K1 / K2;
    
    std::cout << "K1:\t" << K1 << std::endl;
    std::cout << "K2:\t" << K2 << std::endl;
    std::cout << "diffX:\t" << diffX << std::endl;
    
    if(!std::isfinite(K1)){
      // checked it on 07/05:
      // if the solution 'x' tends to infinity, 'K2' tends to 0, and 'K1' tends to 0 very slowly.
      // then we can set the one sided p value as 0 (instead of setting converge = F)
      x = arma::sign(t_Stat) * arma::datum::inf;
      K2 = 0;
      break;
    }
    
    if(arma::sign(K1) != arma::sign(oldK1)){
      while(std::abs(diffX) > std::abs(oldDiffX) - tol){
        diffX = diffX / 2;
      }
    }
    if(std::abs(diffX) < tol) break;
    
    x = oldX + diffX;
    
  }
  
  if(iter == maxiter - 1) 
    converge = false;
  
  Rcpp::List yList = Rcpp::List::create(Rcpp::Named("root") = x,
                                        Rcpp::Named("iter") = iter,
                                        Rcpp::Named("converge") = converge,
                                        Rcpp::Named("K2") = K2);
  return yList;
}

double fastGet_Saddle_Prob(double t_Stat,
                           double t_zeta,
                           double t_K2,
                           double t_Ratio0,
                           arma::mat t_muMat,
                           arma::mat t_cMat,
                           double t_m1,          // sum(muMat * cMat)
                           bool t_lowerTail)
{
  double k1 = K0(t_zeta, t_muMat, t_cMat, t_m1) + 1/2 * pow(t_zeta, 2) * t_Ratio0;
  double k2 = t_K2;
  double pval = 0;
  if(std::isfinite(k1) && std::isfinite(k2))
  {
    double w = arma::sign(t_zeta) * sqrt(2 * (t_zeta * t_Stat - k1));
    double v = t_zeta * sqrt(t_K2);
    
    double Z = w + 1/w * log(v/w);
    pval = arma::normcdf(arma::sign(t_lowerTail-0.5) * Z);
  }
  
  return pval;
}

// add partial normal approximation to speed up the SPA
Rcpp::List fastSaddle_Prob(double t_Stat,
                           double t_VarP,
                           double t_VarW,
                           double t_Ratio0,      // Ratio of variance (G==0)
                           arma::vec t_K1roots,  // 2 x 1
                           arma::vec t_adjGVec1, // N1 x 1, where N1 is length(G!=0)
                           arma::mat t_muMat1,   // N1 x J
                           arma::mat t_iRMat1)   // N1 x (J-1)
{
  int J = t_muMat1.n_cols;
  t_muMat1 = t_muMat1.cols(0, J-2);
  double adjStat = t_Stat / sqrt(t_VarP);
  int N1 = t_muMat1.n_rows;
  
  double sqrtVarW = sqrt(t_VarW);
  
  arma::mat cMat(N1, J-1);
  for(int i = 0; i < N1; i ++){
    for(int j = 0; j < J-1; j ++){
      cMat(i,j) = t_adjGVec1(i) * t_iRMat1(i,j) / sqrtVarW;
    }
  }
  
  double m1 = arma::accu(t_muMat1 % cMat);
  
  Rcpp::List outUni1 = fastgetroot_K1(std::abs(adjStat), std::min(t_K1roots(0), 5.0), 
                                      t_Ratio0, t_muMat1, cMat, m1);
  Rcpp::List outUni2 = fastgetroot_K1(-1 * std::abs(adjStat), std::max(t_K1roots(1), -5.0), 
                                      t_Ratio0, t_muMat1, cMat, m1);
  
  bool converge = false;
  double pval = 0; 
  arma::vec K1roots;
  
  bool outUnit1Converge = outUni1["converge"];
  bool outUnit2Converge = outUni2["converge"];
  
  if(outUnit1Converge == true && outUnit2Converge == true){
    
    double p1 = fastGet_Saddle_Prob(std::abs(adjStat), outUni1["root"], 
                                    outUni1["K2"], t_Ratio0, t_muMat1, cMat, m1, false);
    
    double root = outUni1["root"];
    double K2 = outUni1["K2"];
    std::cout << "outUni1:\t" << root << "\t" << K2 << std::endl;
    std::cout << "p1:\t" << p1 << std::endl;
    
    double p2 = fastGet_Saddle_Prob(-1 * std::abs(adjStat), outUni2["root"], 
                                    outUni2["K2"], t_Ratio0, t_muMat1, cMat, m1, false);
    
    root = outUni2["root"];
    K2 = outUni2["K2"];
    std::cout << "outUni2:\t" << root << "\t" << K2 << std::endl;
    std::cout << "p2:\t" << p2 << std::endl;
    
    pval = p1 + p2;
    
    converge = true;
    K1roots = {outUni1["root"], outUni2["root"]};
  }else{
    std::cout << "SPA does not converge, use normal approximation p value." << std::endl;
    pval = 2 * arma::normcdf(-1 * std::abs(adjStat));
    K1roots = t_K1roots;
  }
  
  if(! std::isfinite(pval)){
    std::cout << "SPA does not give a valid p value, use normal approximation p value." << std::endl;
    pval = 2 * arma::normcdf(-1 * std::abs(adjStat));
    converge = false;
  }
  
  Rcpp::List yList = Rcpp::List::create(Rcpp::Named("pval") = pval,
                                        Rcpp::Named("converge") = converge,
                                        Rcpp::Named("K1roots") = K1roots);
  return yList;
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

// get RPsiP: (J-1) x (J-1) x n 
// Only used in getVarWFast(): 

