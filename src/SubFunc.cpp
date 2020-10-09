
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <sys/time.h>

#include "SubFunc.hpp"

using namespace Rcpp;
using namespace std;

int countLine(string t_file)
{
  int count = 0;
  string junk;
  string errInfo = "ERROR: Unable to open file: " + t_file; 
  
  ifstream ifile(t_file);
  if(!ifile.is_open())
    stop(errInfo);
  
  while(getline(ifile, junk))
    count++;
  
  return count;
}

arma::vec getTime(){
  arma::vec Time(2, arma::fill::zeros);
  struct timeval time;
  Time(0) = 0;
  if(!gettimeofday(&time,NULL))
    Time(0) = (double)time.tv_sec + (double)time.tv_usec * .000001;
  Time(1) = (double)clock() / CLOCKS_PER_SEC;
  return Time; 
}

void printTime(arma::vec t1, arma::vec t2, std::string message){
  double wallTime = t2(0) - t1(0);
  double cpuTime = t2(1) - t1(1);
  if(wallTime < 60){
    Rprintf ("It took %f seconds (%f CPU seconds) to %s.\n", 
             wallTime, cpuTime, message.c_str());
  }else if(wallTime < 3600){
    Rprintf ("It took %f minutes (%f CPU minutes) to %s.\n", 
             wallTime/60, cpuTime/60, message.c_str());
  }else{
    Rprintf ("It took %f hours (%f CPU hours) to %s.\n", 
             wallTime/3600, cpuTime/3600, message.c_str());
  }
}

double getinvStd(double t_freq)
{
  double Std = sqrt(2 * t_freq * (1-t_freq));
  if(Std == 0)
    return 0;
  else
    return 1/Std;
}

//http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

double calCV(arma::vec t_xVec){
  int n = t_xVec.size();
  double Mean = arma::mean(t_xVec);
  double Sd = arma::stddev(t_xVec);
  double CV = (Sd/Mean)/n;
  return CV;
}

arma::vec nb(int n){
  return(Rcpp::rbinom(n,1,0.5));
}

// C++ version of which(). Note: start from 0, not 1 
std::vector<int> whichCPP(Rcpp::StringVector strVec, 
                          std::string strValue)
{
  std::vector<int> indexVec;
  for(int i = 0; i < strVec.size(); i++){
    if(std::string(strVec(i))==strValue)
      indexVec.push_back(i);
  }
  return(indexVec);
}

// duplicate each row for (J-1) times: n x p -> n(J-1) x p
arma::mat getCovaMat(arma::mat Cova,   // matrix: n x p
                     int n, int J, int p)      
{
  arma::mat CovaMat(n*(J-1), p);
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      CovaMat.row(index) = Cova.row(i);
      index++;
    }
  }
  return(CovaMat);
}

// convert: n(J-1) x 1 -> n x (J-1)
arma::mat convert2(arma::vec xVec, // n(J-1) x 1 
                   int n, int J)
{
  arma::mat xMat(n,(J-1));
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      xMat(i,j) = xVec(index);
      index++;
    }
  }
  return(xMat);
}

// convert: n x (J-1) -> n(J-1) x 1
arma::vec convert1(arma::mat xMat, // matrix: n x (J-1)
                   int n, int J) 
{
  arma::vec xVec(n*(J-1));
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      xVec(index) = xMat(i,j);
      index++;
    }
  }
  return(xVec);
}

// sum up each row: n1 x n2 matrix -> n1 x 1 vector
arma::vec getRowSums(arma::mat t_xMat)
{
  int n1 = t_xMat.n_rows;
  int n2 = t_xMat.n_cols;
  arma::vec y1Vec(n1, arma::fill::zeros);
  for(int i = 0; i < n1; i++){
    for(int j = 0; j < n2; j++){
      y1Vec(i) += t_xMat(i,j);
    }
  }
  return(y1Vec);
}

// get a list for p value calculation in step 2
// [[Rcpp::export]]
Rcpp::List getobjP(arma::mat t_Cova,     // matrix: n x p
                   arma::mat t_yMat,
                   arma::mat t_muMat,    // matrix: n x J
                   arma::mat t_iRMat)    // matrix: n x (J-1)
{
  int n = t_muMat.n_rows;
  int J = t_muMat.n_cols;
  int p = t_Cova.n_cols;
  
  // output for Step 2
  arma::mat XR_Psi_R(p, n*(J-1));                // p x n(J-1)
  arma::mat CovaMat = getCovaMat(t_Cova, n, J, p); // n(J-1) x p
  for(int k = 0; k < p; k++){
    arma::mat xMat = convert2(CovaMat.col(k), n, J);
    arma::vec temp = convert1(getPsixMat(xMat / t_iRMat, t_muMat) / t_iRMat, n, J);
    XR_Psi_R.row(k) = temp.t();
  }
  // arma::mat XXR_Psi_RX = CovaMat * inv(XR_Psi_R * CovaMat);             // (n(J-1) x p) * (p x p) = n(J-1) x p
  arma::mat XXR_Psi_RX_new = t_Cova * inv(XR_Psi_R * CovaMat);               // (n x p) * (p x p) = n x p
  
  // sum each (J-1) rows to 1 row: p x n(J-1) -> p x n
  arma::mat XR_Psi_R_new = sumCols(XR_Psi_R, J);      // p x n
  arma::mat ymuMat = t_yMat - t_muMat;                    // n x J
  arma::mat RymuMat = ymuMat.cols(0, J-2) / t_iRMat;    // n x (J-1): R %*% (y - mu)
  arma::mat RymuVec = sumCols(RymuMat, J);            // n x 1
  // arma::cube RPsiR = getRPsiR(muMat, iRMat, n, J, p); // (J-1) x (J-1) x n 
  arma::vec RPsiR = getRPsiR(t_muMat, t_iRMat, n, J, p); // (J-1) x (J-1) x n 
  
  Rcpp::List objP = List::create(Named("n")=n,
                                 Named("J")=J,
                                 Named("p")=p,
                                 Named("XXR_Psi_RX_new") = XXR_Psi_RX_new,
                                 Named("XR_Psi_R_new") = XR_Psi_R_new,           
                                 Named("RymuVec") = RymuVec,
                                 Named("RPsiR") = RPsiR,
                                 Named("muMat") = t_muMat,
                                 Named("iRMat") = t_iRMat);
  return(objP);
}

// outMat = PsiMat %*% xMat, PsiMat is determined by muMat
arma::mat getPsixMat(arma::mat t_xMat,    // matrix: n x (J-1)
                     arma::mat t_muMat)   // matrix: n x J
{
  int n = t_muMat.n_rows;
  int J = t_muMat.n_cols;
  
  arma::mat Psi_xMat(n, J-1);
  // loop for samples
  for(int i = 0; i < n; i++){
    arma::rowvec muVec(J-1);
    for(int j = 0; j < J-1; j++){
      Psi_xMat(i,j) = t_muMat(i,j) * t_xMat(i,j);
      muVec(j) = t_muMat(i,j);
    }
    double sum_mu_x = sum(Psi_xMat.row(i));
    Psi_xMat.row(i) -= muVec * sum_mu_x; 
  }
  return(Psi_xMat);
}

// sum each (J-1) cols to 1 col: p x n(J-1) -> p x n (OR) p x (J-1) -> p x 1
arma::mat sumCols(arma::mat t_xMat,
                  int J)
{
  int n = t_xMat.n_cols / (J-1);
  int p = t_xMat.n_rows;
  arma::mat outMat(p, n, arma::fill::zeros);
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      outMat.col(i) += t_xMat.col(index);
      index++;
    }
  }
  return(outMat);
}

// get RPsiP: (J-1) x (J-1) x n 
// Only used in getVarWFast(): 
arma::vec getRPsiR(arma::mat t_muMat,
                   arma::mat t_iRMat,
                   int t_n, int t_J, int t_p)   
{
  // arma::cube RPsiR(J-1, J-1, n, arma::fill::zeros);
  arma::vec RPsiRVec(t_n, arma::fill::zeros);
  arma::mat muRMat = t_muMat.cols(0, t_J-2) / t_iRMat;
  for(int i = 0; i < t_n; i++){
    for(int j1 = 0; j1 < t_J-1; j1++){
      // RPsiR(j1,j1,i) += muRMat(i,j1) / iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      RPsiRVec(i) += muRMat(i,j1) / t_iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      for(int j2 = j1+1; j2 < t_J-1; j2++){
        // RPsiR(j1,j2,i) -= muRMat(i,j1) * muRMat(i,j2);
        RPsiRVec(i) -= 2* muRMat(i,j1) * muRMat(i,j2);
      }
    }
  }
  // return(RPsiR);
  return(RPsiRVec);
}

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat)
{
  double innerProd = arma::accu(x1Mat % x2Mat);
  return(innerProd);
}

// [[Rcpp::export]]
Rcpp::List outputadjGFast(arma::vec GVec,
                          Rcpp::List objP)
{
  arma::vec adjGVec = getadjGFast(GVec, objP["XXR_Psi_RX_new"], objP["XR_Psi_R_new"], objP["n"], objP["p"]);
  double Stat = getStatFast(adjGVec, objP["RymuVec"]);
  double VarW = getVarWFast(adjGVec, objP["RPsiR"]);
  Rcpp::List outList = List::create(Named("adjGVec")=adjGVec,
                                    Named("Stat")=Stat,           
                                    Named("VarW")=VarW);
  
  return(outList);
}

arma::vec getadjGFast(arma::vec GVec,
                      arma::mat XXR_Psi_RX_new,   // XXR_Psi_RX_new ( n x p )
                      arma::mat XR_Psi_R_new,     // XR_Psi_R_new ( p x n ), sum up XR_Psi_R ( p x n(J-1) ) for each subject 
                      int n, int p)
{
  // To increase computational efficiency when lots of GVec elements are 0
  arma::vec XR_Psi_RG1(p, arma::fill::zeros);
  for(int i = 0; i < n; i++){
    if(GVec(i) != 0){
      XR_Psi_RG1 += XR_Psi_R_new.col(i) * GVec(i);
    }
  }
  
  arma::vec adjGVec = GVec - XXR_Psi_RX_new * XR_Psi_RG1;
  return(adjGVec);
}

double getStatFast(arma::vec GVec,         // n x 1
                   arma::vec RymuVec)      // n x 1: row sum of the n x (J-1) matrix R %*% (yMat - muMat)
{
  int n = GVec.size();
  double Stat = 0;
  for(int i = 0; i < n; i++){
    if(GVec(i) != 0){
      Stat += GVec(i) * RymuVec(i);
    }
  }
  return(Stat);
}

double getVarWFast(arma::vec adjGVec,  // n x 1
                   arma::vec RPsiRVec) // n x 1
{
  int n = adjGVec.size();
  double VarW = 0;
  for(int i = 0; i < n; i++){
    VarW += RPsiRVec(i) * adjGVec(i) * adjGVec(i);
  }
  return(VarW);
}

// yMat: matrix with dim of n x J
// [[Rcpp::export]]
arma::mat getyMatR(arma::mat yVec, int n, int J)
{
  arma::mat yMat(n, J, arma::fill::zeros);
  for(int i = 0; i < n; i++)
    yMat(i, yVec(i)-1) = 1;
  return(yMat);
}
