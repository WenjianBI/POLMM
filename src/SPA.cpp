// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include "SPA.hpp"
#include <math.h>

// arma::vec MAIN_TEST(arma::vec t_GVec,
//                     arma::mat t_XR_Psi_R_new,
//                     arma::mat t_XXR_Psi_RX_new,
//                     arma::mat t_RymuVec,
//                     arma::vec t_RPsiR,
//                     double t_r,
//                     arma::vec t_K1roots,
//                     double t_SPAcutoff,
//                     arma::mat t_iRMat,
//                     arma::mat t_muMat)
// {
//   arma::uvec posG1 = arma::find(t_GVec != 0);
//   arma::mat XR_Psi_RG1 = t_XR_Psi_R_new.cols(posG1) * t_GVec.elem(posG1);
//   arma::vec adjGVec = t_GVec - t_XXR_Psi_RX_new * XR_Psi_RG1;
//   double Stat = sum(adjGVec % t_RymuVec);
//   arma::vec VarWVec = t_RPsiR % pow(adjGVec, 2);
//   double VarW = sum(VarWVec);
//   double VarP = VarW * t_r;
//   double z = std::abs(Stat) / sqrt(VarP);
//   double beta = Stat / VarP;
//   double pvalNorm = 2 * arma::normcdf(-1 * z);
//   double pvalSPA = pvalNorm;
//   arma::vec K1roots = t_K1roots;
//   
//   if(z > t_SPAcutoff){
//     double VarW1 = sum(VarWVec(posG1));
//     double VarW0 = VarW - VarW1;
//     double Ratio0 = VarW0 / VarW;
//     
//     Rcpp::List resSPA = fastSaddle_Prob(Stat, VarP, VarW, Ratio0, K1roots,
//                                         adjGVec.elem(posG1), t_muMat.rows(posG1), t_iRMat.rows(posG1));
//     pvalSPA = resSPA["pval"];
//     arma::vec K1rootsTemp = resSPA["K1roots"];
//     K1roots = K1rootsTemp;
//   }
//   
//   arma::vec resTest = {pvalNorm, pvalSPA, beta};;
//   return(resTest);
// }


double K0(double t_x,
          arma::mat t_muMat,     // N x (J-1)
          arma::mat t_cMat,      // N x (J-1)
          double t_m1)           // sum(muMat * cMat)
{
  arma::mat temp1Mat = - t_muMat + t_muMat % exp(t_cMat * t_x);
  arma::vec temp1Vec = log(1 + arma::sum(temp1Mat, 1));   // arma::sum(Mat, 1) is rowSums()
  double y = sum(temp1Vec) - t_m1 * t_x;
  return(y);
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
  yVec(1) = sum(temp2Vec / temp1Vec) - t_m1;
  // yMat[i,2] = sum((temp3Vec*temp1Vec-temp2Vec^2)/temp1Vec^2, na.rm=TRUE);
  yVec(2) = sum((temp3Vec * temp1Vec - pow(temp2Vec,2)) / pow(temp1Vec,2));
  
  return(yVec);
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
    K1 = K12Vec(1) - t_Stat + t_Ratio0 * x;
    K2 = K12Vec(2) + t_Ratio0;
    
    diffX = -1 * K1 / K2;
    if(std::isfinite(K1)){
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
  return(yList);
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
  
  return(pval);
}

// add partial normal approximation to speed up the SPA
Rcpp::List fastSaddle_Prob(double t_Stat,
                           double t_VarP,
                           double t_VarW,
                           double t_Ratio0,      // Ratio of variance (G==0)
                           arma::vec t_K1roots,
                           arma::vec t_adjGVec1, // N1 x 1, where N1 is length(G!=0)
                           arma::mat t_muMat1,   // N1 x (J-1)
                           arma::mat t_iRMat1)   // N1 x (J-1)
{
  double adjStat = t_Stat / sqrt(t_VarP);
  int N1 = t_muMat1.n_rows;
  int J = t_muMat1.n_cols + 1;
  
  double sqrtVarW = sqrt(t_VarW);
  
  arma::mat cMat(N1, J-1);
  for(int i = 0; i < N1; i ++){
    for(int j = 0; j < J-1; j ++){
      cMat(i,j) = t_adjGVec1(i) * t_iRMat1(i,j) / sqrtVarW;
    }
  }
  
  double m1 = arma::accu(t_muMat1 * cMat);
  
  Rcpp::List outUni1 = fastgetroot_K1(std::abs(adjStat), std::min(t_K1roots[1], 5.0), 
                                      t_Ratio0, t_muMat1, cMat, m1);
  Rcpp::List outUni2 = fastgetroot_K1(-1 * std::abs(adjStat), std::max(t_K1roots[2], -5.0), 
                                      t_Ratio0, t_muMat1, cMat, m1);
  
  bool converge = false;
  double pval = 0; 
  arma::vec K1roots;
  
  bool outUnit1Converge = outUni1["converge"];
  bool outUnit2Converge = outUni2["converge"];
  
  if(outUnit1Converge == true & outUnit2Converge == true){
    
    double p1 = fastGet_Saddle_Prob(std::abs(adjStat), outUni1["root"], 
                                    outUni1["K2"], t_Ratio0, t_muMat1, cMat, m1, false);
    
    double p2 = fastGet_Saddle_Prob(-1 * std::abs(adjStat), outUni2["root"], 
                                    outUni2["K2"], t_Ratio0, t_muMat1, cMat, m1, false);
    
    pval = p1 + p2;
    
    converge = true;
    K1roots = {outUni1["root"], outUni2["root"]};
  }else{
    Rprintf("SPA does not converge, use normal approximation p value.");
    pval = 2 * arma::normcdf(-1 * std::abs(adjStat));
    K1roots = t_K1roots;
  }
  
  if(! std::isfinite(pval)){
    Rprintf("SPA does not give a valid p value, use normal approximation p value.");
    pval = 2 * arma::normcdf(-1 * std::abs(adjStat));
    converge = false;
  }
  
  Rcpp::List yList = Rcpp::List::create(Rcpp::Named("pval") = pval,
                                        Rcpp::Named("converge") = converge,
                                        Rcpp::Named("K1roots") = K1roots);
  return(yList);
}


