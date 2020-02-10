#include "genoClassExt.hpp"

#include <RcppArmadillo.h>

#include <vector>
#include <cmath>
#include <RcppArmadilloExtensions/sample.h> // sample

arma::vec nb(int n){
  return(Rcpp::rbinom(n,1,0.5));
}

double get_cpu_time(){
  return (double)clock() / CLOCKS_PER_SEC;
}

double calCV(arma::vec xVec){
  int n = xVec.size();
  double Mean = arma::mean(xVec);
  double Sd = arma::stddev(xVec);
  double CV = (Sd/Mean)/n;
  return(CV);
}

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat)
{
  double innerProd = arma::accu(x1Mat % x2Mat);
  return(innerProd);
}

std::vector<int> whichCPP(Rcpp::StringVector chrVec, std::string excludechr)
{
  std::vector<int> indexSNPs;
  for(int i = 0; i < chrVec.size(); i++){
    if(std::string(chrVec(i))==excludechr){
      indexSNPs.push_back(i);
    }
  }
  return(indexSNPs);
}

arma::mat getadjGFast(arma::vec GVec,
                      arma::mat XXR_Psi_RX,     
                      arma::mat XR_Psi_R_new, // XR_Psi_R_new ( p x n ), sum up XR_Psi_R ( p x n(J-1) ) for each subject 
                      int n, int J, int p)
{
  arma::mat adjGMat(n, J-1);
  // To increase computational efficiency when lots of GVec elements are 0
  arma::vec XR_Psi_RG1(p, arma::fill::zeros);
  for(int i = 0; i < n; i++){
    if(GVec(i)!=0){
      XR_Psi_RG1 += XR_Psi_R_new.col(i) * GVec(i);
    }
  }
  
  arma::vec TempVec = XXR_Psi_RX * XR_Psi_RG1;   // n(J-1) x 1
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      adjGMat(i,j) = GVec(i) - TempVec(index);
      index++;
    }
  }
  return(adjGMat);
}

double getStatFast(arma::vec GVec,         // n x 1
                   arma::mat RymuMat,      // n x (J-1): R %*% (y - mu)
                   int n, int J)
{
  double Stat = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      Stat += GVec(i) * RymuMat(i,j);
    }
  }
  return(Stat);
}

double getVarWFast(arma::mat adjGMat,  // n x (J-1)
                   arma::cube RPsiR,   // (J-1) x (J-1) x n
                   int n, int J)
{
  double VarW = 0;
  for(int i = 0; i < n; i++){
    for(int j1 = 0; j1 < J-1; j1++){
      VarW += RPsiR(j1,j1,i) * adjGMat(i,j1) * adjGMat(i,j1);
      for(int j2 = j1+1; j2 < J-1; j2++){
        VarW += 2 * RPsiR(j1,j2,i) * adjGMat(i,j1) * adjGMat(i,j2);
      }
    }
  }
  return(VarW);
}

// [[Rcpp::export]]
Rcpp::List outputadjGFast(arma::vec GVec,
                          Rcpp::List objP)
{
  arma::mat adjGMat = getadjGFast(GVec, objP["XXR_Psi_RX"], objP["XR_Psi_R_new"], objP["n"], objP["J"], objP["p"]);
  double Stat = getStatFast(GVec, objP["RymuMat"], objP["n"], objP["J"]);
  double VarW = getVarWFast(adjGMat, objP["RPsiR"], objP["n"], objP["J"]);
  
  Rcpp::List outList = List::create(Named("adjGMat")=adjGMat,
                                    Named("Stat")=Stat,           
                                    Named("VarW")=VarW);
  
  return(outList);
}


class nullModelClass
{
private:
  
  // used as dimension
  int n, J, p, M;
  genoClass* ptrGeno;
  Rcpp::List initParList, controlList;
  
  // constant vector and matrix over the process
  arma::Col<int> yVec; // should start from 1, not 0
  arma::mat yMat, Cova, CovaMat, GMatRatio;
  
  // parameters
  arma::vec beta, bVec, eps;
  arma::vec beta0, eps0;
  double tau, tau0;
  double diffBeta, difftau, diffeps;
  
  // arguments
  int maxiter, maxiterPCG, maxiterEps, tracenrun, seed, nSNPsVarRatio;
  double tolBeta, tolTau, tolPCG, tolEps, minMafVarRatio, CVcutoff;
  bool LOCO;
  Rcpp::List LOCOList;
  
  // working vectors/matrix
  arma::mat WMat, muMat, mMat, nuMat, iRMat, YMat, iSigma_CovaMat, iSigmaX_XSigmaX;
  arma::vec eta, iSigma_YVec, iSigma_VPYVec;
  
  arma::mat TraceRandMat, V_TRM, iSigma_V_TRM;
  
  // functions in setNullModel()
  void setParList(Rcpp::List);    // set up parameter list
  void setcontrolList(Rcpp::List);
  void setArray();               // set up dimension 
  arma::mat getCovaMat();
  arma::mat getyMat();
  void getTraceRandMat();
  arma::vec ZMat(arma::vec), tZMat(arma::vec);
  
  // functions in fitNullModel()
  void updateMats();
  arma::mat getiPsixMat(arma::mat);
  arma::vec getRowSums(arma::mat);
  arma::vec convert1(arma::mat xMat);
  arma::mat convert2(arma::vec xMat);
  void getPCGofSigmaAndCovaMat(arma::mat, arma::mat&, string);
  void getPCGofSigmaAndVector(arma::vec, arma::vec&, string);
  void updateParaConv(string);
  void updatePara(string);
  void updateEpsOneStep();
  void updateEps();
  arma::mat getSigmaxMat(arma::mat, string);
  arma::cube getInvBlockDiagSigma();
  arma::mat solverBlockDiagSigma(arma::cube&, arma::mat& xMat);
  void updateTau();
  
  // functions if LOCO = TRUE
  arma::mat getPsixMat(arma::mat);
  Rcpp::List getLOCO(string, std::vector<int>);
  arma::mat getVarRatio(std::vector<int>, string, Rcpp::List);
  arma::rowvec getVarOneSNP(arma::vec, string, Rcpp::List);
  double getVarP(arma::mat, string);
  Rcpp::List getobjP();
  void getRPsiR(arma::cube&);
  
  // functions if LOCO = FALSE
  arma::mat getVarRatio(arma::mat, Rcpp::List);
  
public:
  
  void setNullModel(arma::mat, arma::Col<int>, genoClass*, arma::vec, arma::vec, arma::vec, double, arma::mat, Rcpp::List);
  void fitNullModel();
  Rcpp::List getNullModel();
  
};

void nullModelClass::fitNullModel()
{
  // initial vector
  updateMats();
  
  // start iteration
  cout << "Start iteration ....." << endl;
  int iter = 0;
  while(iter < maxiter)
  {
    cout << "iter " << iter << endl;
    
    // update other parameters until converge
    updateParaConv("n");
    
    // update tau
    tau0 = tau;
    updateTau();
    
    if(isnan(tau))
      stop("Parameter tau is NA.");
    
    cout << "beta: " << endl << beta << endl;
    cout << "tau: " << tau << endl << endl;
    iter++;
    
    double diffTau = abs(tau - tau0)/(abs(tau) + abs(tau0) + tolTau);
    if(diffTau < tolTau)
      break;
  }
  
  if(LOCO){
    Rcpp::StringVector chrVec = ptrGeno->getchrVec();
    Rcpp::StringVector uniqchr = unique(chrVec);
    cout << uniqchr << endl;
    for(int i = 0; i < uniqchr.size(); i++){
      string excludechr = string(uniqchr(i));
      std::vector<int> indexSNPs = whichCPP(chrVec, excludechr);
      indexSNPs = Rcpp::RcppArmadillo::sample(indexSNPs, indexSNPs.size(), FALSE);
      cout << endl << "Leave One Chromosome Out: Chr " << excludechr << endl;
      cout << "Number of SNPs in this chromosome:\t" << indexSNPs.size() << endl;
      Rcpp::List tempLOCO = getLOCO(excludechr, indexSNPs);
      LOCOList[excludechr] = tempLOCO;
    }
  }else{
    Rcpp::List objP = getobjP();
    // output variance matrix
    arma::mat VarRatioMat = getVarRatio(GMatRatio, objP);
    double VarRatio = arma::mean(VarRatioMat.col(4));
    
    Rcpp::List temp = List::create(Named("objP")=objP,
                                   Named("VarRatioMat")=VarRatioMat,
                                   Named("VarRatio")=VarRatio);
    LOCOList["LOCO=F"] = temp;
  }
}

// outMat = PsiMat %*% xMat, PsiMat is determined by muMat
arma::mat nullModelClass::getPsixMat(arma::mat xMat)    // matrix with dim of n x (J-1)
{
  arma::mat Psi_xMat(n, J-1);
  // loop for samples
  for(int i = 0; i < n; i++){
    arma::rowvec muVec(J-1);
    for(int j = 0; j < J-1; j++){
      Psi_xMat(i,j) = muMat(i,j) * xMat(i,j);
      muVec(j) = muMat(i,j);
    }
    double sum_mu_x = sum(Psi_xMat.row(i));
    Psi_xMat.row(i) -= muVec * sum_mu_x; 
  }
  return(Psi_xMat);
}

void nullModelClass::getRPsiR(arma::cube& RPsiR)   // (J-1) x (J-1) x n
{
  arma::mat muRMat = muMat.cols(0, J-2) / iRMat;
  for(int i = 0; i < n; i++){
    for(int j1 = 0; j1 < J-1; j1++){
      RPsiR(j1,j1,i) += muRMat(i,j1) / iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      for(int j2 = j1+1; j2 < J-1; j2++){
        RPsiR(j1,j2,i) -= muRMat(i,j1) * muRMat(i,j2);
        // RPsiR(j2,j1,i) -= muRMat(i,j1) * muRMat(i,j2); // v4.1: RPsiR is only used in function getVarWFast, 
      }
    }
  }
}

// get a list for p value calculation in step 2
Rcpp::List nullModelClass::getobjP()
{
  // output for Step 2
  arma::mat XR_Psi_R(p, n*(J-1));
  arma::mat xMat(n, J-1);
  arma::vec temp(n*(J-1));
  for(int k = 0; k < p; k++){
    xMat = convert2(CovaMat.col(k));
    temp = convert1(getPsixMat(xMat / iRMat) / iRMat);
    XR_Psi_R.row(k) = temp.t();
  }
  arma::mat XXR_Psi_RX = CovaMat * inv(XR_Psi_R * CovaMat);
  iSigmaX_XSigmaX = iSigma_CovaMat * inv(CovaMat.t() * iSigma_CovaMat);
  
  arma::cube RPsiR(J-1, J-1, n, arma::fill::zeros);
  getRPsiR(RPsiR);
  
  arma::mat XR_Psi_R_new(p, n, arma::fill::zeros);  // p x n
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      XR_Psi_R_new.col(i) += XR_Psi_R.col(index);
      index++;
    }
  }
  
  arma::mat ymuMat = yMat - muMat;
  arma::mat RymuMat = ymuMat.cols(0, J-2) / iRMat; // n x (J-1): R %*% (y - mu)
  
  Rcpp::List objP = List::create(Named("n")=n,
                                 Named("J")=J,
                                 Named("p")=p,
                                 Named("XXR_Psi_RX")=XXR_Psi_RX,
                                 Named("RPsiR")=RPsiR,
                                 Named("XR_Psi_R_new")=XR_Psi_R_new,           
                                 Named("RymuMat")=RymuMat,
                                 Named("muMat")=muMat,
                                 Named("iRMat")=iRMat);
  
  return(objP);
}

double nullModelClass::getVarP(arma::mat adjGMat,
                               string excludechr)
{
  arma::vec adjGVec = convert1(adjGMat);
  arma::vec iSigmaGVec(n*(J-1), arma::fill::zeros);
  getPCGofSigmaAndVector(adjGVec, iSigmaGVec, excludechr);
  double VarP = as_scalar(adjGVec.t() * (iSigmaGVec - iSigmaX_XSigmaX * (CovaMat.t() * iSigmaGVec)));
  return(VarP);
}

arma::rowvec nullModelClass::getVarOneSNP(arma::vec GVec,
                                          string excludechr,
                                          Rcpp::List objP)
{
  arma::rowvec VarOut(5);
  double AF = sum(GVec) / GVec.size() / 2;
  if(AF > 0.5)
    AF = 1 - AF;
  
  Rcpp::List adjGList = outputadjGFast(GVec, objP);
  arma::mat adjGMat = adjGList["adjGMat"];
  double Stat = adjGList["Stat"];
  double VarW = adjGList["VarW"];
  double VarP = getVarP(adjGMat, excludechr);
  
  VarOut(0) = AF;
  VarOut(1) = Stat;
  VarOut(2) = VarW;
  VarOut(3) = VarP;
  VarOut(4) = VarP/VarW;
  return(VarOut);
}

arma::mat nullModelClass::getVarRatio(arma::mat GMatRatio,
                                      Rcpp::List objP)
{
  arma::vec GVec(n);
  arma::rowvec VarOneSNP(5);
  
  arma::mat VarRatioMat(nSNPsVarRatio, 5);
  arma::mat newVarRatio(10, 5);
  
  // Rcpp::List adjGList = outputadjGFast(GVec, objP);
  
  int index = 0;
  int indexTot = 0;
  while(index < nSNPsVarRatio){
    indexTot++;
    GVec = GMatRatio.col(index);
    VarOneSNP = getVarOneSNP(GVec, "n", objP);
    VarRatioMat.row(index) = VarOneSNP;
    index++;
  }
  
  arma::vec VarRatio = VarRatioMat.col(4);
  double CV = calCV(VarRatio);
  cout << "nSNPs for CV: " << index << endl;
  cout << "CV: " << CV << endl;
  
  while(CV > CVcutoff && VarRatioMat.n_rows < 100){
    int indexTemp = 0;
    while(indexTemp < 10){
      indexTot++;
      GVec = GMatRatio.col(indexTot);
      VarOneSNP = getVarOneSNP(GVec, "n", objP);
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

arma::mat nullModelClass::getVarRatio(std::vector<int> indexSNPs,
                                      string excludechr,
                                      Rcpp::List objP)
{
  arma::vec alleleFreqVec = ptrGeno->getalleleFreqVec();
  
  arma::vec GVec(n);
  arma::rowvec VarOneSNP(5);
  
  arma::mat VarRatioMat(nSNPsVarRatio, 5);
  arma::mat newVarRatio(10, 5);
  
  // Rcpp::List adjGList = outputadjGFast(GVec, objP);
  
  int index = 0;
  int indexTot = 0;
  while(index < nSNPsVarRatio){
    int m = indexSNPs[indexTot];
    indexTot++;
    if(alleleFreqVec(m) > minMafVarRatio && alleleFreqVec(m) < 1-minMafVarRatio){
      ptrGeno->getOneMarker(m, 0, &GVec);
      VarOneSNP = getVarOneSNP(GVec, excludechr, objP);
      VarRatioMat.row(index) = VarOneSNP;
      index++;
    }
  }
  
  arma::vec VarRatio = VarRatioMat.col(4);
  double CV = calCV(VarRatio);
  cout << "nSNPs for CV: " << index << endl;
  cout << "CV: " << CV << endl;
  
  while(CV > CVcutoff && VarRatioMat.n_rows < 100){
    int indexTemp = 0;
    while(indexTemp < 10){
      int m = indexSNPs[indexTot]; // round double to int
      indexTot++;
      if(alleleFreqVec(m) > minMafVarRatio && alleleFreqVec(m) < 1-minMafVarRatio){
        ptrGeno->getOneMarker(m, 0, &GVec);
        VarOneSNP = getVarOneSNP(GVec, excludechr, objP);
        newVarRatio.row(indexTemp) = VarOneSNP;
        index++;
        indexTemp++;
      }
    }
    VarRatioMat.insert_rows(0, newVarRatio);
    arma::vec VarRatio = VarRatioMat.col(4);
    CV = calCV(VarRatio);
    cout << "nSNPs for CV: " << index << endl;
    cout << "CV: " << CV << endl;
  }
  return(VarRatioMat);
}

Rcpp::List nullModelClass::getLOCO(string excludechr,
                                   std::vector<int> indexSNPs)
{
  // estimate LOCO parameter
  updateParaConv(excludechr);
  
  // get a list for p value calculation in step 2
  Rcpp::List objP = getobjP();
  
  // output variance matrix
  arma::mat VarRatioMat = getVarRatio(indexSNPs, excludechr, objP);
  double VarRatio = arma::mean(VarRatioMat.col(4));
  
  Rcpp::List outLOCO = List::create(Named("objP")=objP,
                                    Named("VarRatioMat")=VarRatioMat,
                                    Named("VarRatio")=VarRatio);
  return(outLOCO);
}

void nullModelClass::setNullModel(arma::mat CovaR,
                                  arma::Col<int> yVecR,     // should be from 1 to J
                                  genoClass* ptrGenoInput,
                                  arma::vec betaR,
                                  arma::vec bVecR,
                                  arma::vec epsR,
                                  double tauR,
                                  arma::mat GMatRatioR,
                                  Rcpp::List controlListR)
{
  n = CovaR.n_rows;
  p = CovaR.n_cols;
  J = max(yVecR);
  M = ptrGenoInput->getM();
  
  Cova = CovaR;
  CovaMat = getCovaMat();
  yVec = yVecR;
  yMat = getyMat();
  ptrGeno = ptrGenoInput;
  
  beta = betaR;
  bVec = bVecR;
  eps = epsR;
  tau = tauR;
  GMatRatio = GMatRatioR;
  
  initParList = List::create(Named("betaR")=betaR,
                             Named("bVecR")=bVecR,
                             Named("epsR")=epsR,
                             Named("tauR")=tauR);
  
  setcontrolList(controlListR);
  setArray();
  set_seed(seed);
  getTraceRandMat();
}

arma::mat nullModelClass::getCovaMat()
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

// yMat: matrix with dim of n x J
arma::mat nullModelClass::getyMat()
{
  arma::mat yMat(n, J, arma::fill::zeros);
  for(int i = 0; i < n; i++){
    yMat(i, yVec(i)-1) = 1;
  }
  return(yMat);
}

void nullModelClass::setcontrolList(Rcpp::List controlListR)
{
  controlList = controlListR;
  maxiter = controlList["maxiter"];
  maxiterPCG = controlList["maxiterPCG"]; 
  maxiterEps = controlList["maxiterEps"];
  tolBeta = controlList["tolBeta"]; 
  tolTau = controlList["tolTau"]; 
  tolPCG = controlList["tolPCG"]; 
  tolEps = controlList["tolEps"];
  tracenrun = controlList["tracenrun"]; 
  seed = controlList["seed"];
  minMafVarRatio = controlList["minMafVarRatio"];
  nSNPsVarRatio = controlList["nSNPsVarRatio"];
  CVcutoff = controlList["CVcutoff"];
  LOCO = controlList["LOCO"];
}

void nullModelClass::setArray()
{
  WMat.zeros(n,J);    
  muMat.zeros(n,J);   
  mMat.zeros(n,J);
  nuMat.zeros(n,J);
  iRMat.zeros(n,J-1);
  YMat.zeros(n,J-1);
  iSigma_CovaMat.zeros(n*(J-1), p);
  //
  eta.zeros(n);
  iSigma_YVec.zeros(n*(J-1));
  iSigma_VPYVec.zeros(n*(J-1));
  //
  TraceRandMat.zeros(n*(J-1), tracenrun);
  V_TRM.zeros(n*(J-1), tracenrun);
  iSigma_V_TRM.zeros(n*(J-1), tracenrun);
}

void nullModelClass::getTraceRandMat()
{
  for(int itrace = 0; itrace < tracenrun; itrace++){
    cout << "itrace:\t" << itrace << endl;
    arma::vec uVec = nb(n*(J-1));
    uVec = uVec*2 - 1;
    TraceRandMat.col(itrace) = uVec;
    
    double cpuin  = get_cpu_time();
    V_TRM.col(itrace) = tZMat(getKinbVec(ZMat(uVec), ptrGeno, "n"));
    double cpuout  = get_cpu_time();
    cout << "CPU Time  in getKinbVec = " << cpuout - cpuin  << endl;
  }
}

// convert from n(J-1) to n by summing up each (J-1) elements
arma::vec nullModelClass::ZMat(arma::vec xVec)
{
  arma::vec y1Vec(n, arma::fill::zeros);
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      y1Vec(i) += xVec(index);
      index++;
    }
  }
  return(y1Vec);
}

// convert from n to n(J-1) by duplicating each element for (J-1) times
arma::vec nullModelClass::tZMat(arma::vec xVec)
{
  arma::vec y1Vec(n*(J-1));
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      y1Vec(index) = xVec(i);
      index++;
    }
  }
  return(y1Vec);
}


void nullModelClass::updateTau()
{
  arma::vec YVec = convert1(YMat);
  getPCGofSigmaAndCovaMat(CovaMat, iSigma_CovaMat, "n");
  getPCGofSigmaAndVector(YVec, iSigma_YVec, "n"); 
  arma::mat Cova_iSigma_CovaMat = CovaMat.t() * iSigma_CovaMat;
  arma::vec Cova_iSigma_YVec = CovaMat.t() * iSigma_YVec;
  arma::mat iSigmaXXiSigmaX = iSigma_CovaMat * inv(Cova_iSigma_CovaMat);
  arma::vec PYVec = iSigma_YVec - iSigmaXXiSigmaX * Cova_iSigma_YVec;
  arma::vec VPYVec = tZMat(getKinbVec(ZMat(PYVec), ptrGeno, "n"));
  getPCGofSigmaAndVector(VPYVec, iSigma_VPYVec, "n"); 
  arma::vec Cova_iSigma_VPYVec = CovaMat.t() * iSigma_VPYVec;
  arma::vec PVPYVec = iSigma_VPYVec - iSigmaXXiSigmaX * Cova_iSigma_VPYVec;
  double YPVPY = as_scalar(YVec.t() * PVPYVec);
  double YPVPVPY = as_scalar(VPYVec.t() * PVPYVec);
  // The below is to calculate trace
  getPCGofSigmaAndCovaMat(V_TRM, iSigma_V_TRM, "n");
  double tracePV = 0;
  int m = TraceRandMat.n_cols;
  for(int i = 0; i < m; i++){
    arma::vec iSigma_V_TRM_col = iSigma_V_TRM.col(i);
    arma::vec P_V_TRM_col = iSigma_V_TRM_col - iSigmaXXiSigmaX * (CovaMat.t() * iSigma_V_TRM_col);
    tracePV += as_scalar(TraceRandMat.col(i).t() * P_V_TRM_col);
  }
  tracePV /= m;
  // final step
  double deriv = 0.5 * YPVPY - 0.5 * tracePV;
  double AI = 0.5 * YPVPVPY;
  double dtau = deriv / AI;
  tau0 = tau;
  tau = tau0 + dtau;
  while(tau < 0){
    dtau = dtau / 2;
    tau = tau0 + dtau;
  }
  if(tau < 1e-4){
    tau = 0;
  }
}

// update (eta, WMat, muMat, mMat, nuMat, iRMat, YMat) based on (beta, bVec, eps)
void nullModelClass::updateMats()
{
  // update (eta)
  eta = Cova * beta + bVec;
  // update (WMat, muMat, mMat, nuMat)
  double tmpExp, tmpnu0, tmpnu1;
  for(int i = 0; i < n; i++){  // loop for samples
    tmpnu0 = 0;  // eps_0 = -Inf
    for(int j = 0; j < J-1; j++){ // loop from eps_1 to eps_{J-1}
      tmpExp = exp(eps(j) - eta(i));
      tmpnu1 = tmpExp / (1+tmpExp);
      muMat(i,j) = tmpnu1 - tmpnu0;
      WMat(i,j) = tmpnu1 * (1 - tmpnu1);
      mMat(i,j) = tmpnu1 + tmpnu0 - 1;
      nuMat(i,j) = tmpnu1;
      tmpnu0 = tmpnu1;
    }
    int j = J-1;      // eps_J = Inf
    tmpnu1 = 1;
    muMat(i,j) = tmpnu1 - tmpnu0;
    WMat(i,j) = tmpnu1 * (1 - tmpnu1);
    mMat(i,j) = tmpnu1 + tmpnu0 - 1;
    nuMat(i,j) = tmpnu1;
  }
  // update (iRMat)
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      iRMat(i,j) = 1 / (mMat(i,j) - mMat(i,J-1));
    }
  }
  // update (YMat)
  arma::mat xMat = yMat.cols(0,J-2) - muMat.cols(0,J-2);
  arma::mat iPsi_xMat = getiPsixMat(xMat);
  for(int i = 0; i < n; i++){  // loop for samples
    for(int j = 0; j < J-1; j++)
      YMat(i,j) = eta(i) + (iRMat(i,j) * iPsi_xMat(i,j));
  }
}

// outMat = iPsiMat %*% xMat, iPsiMat is determined by muMat
arma::mat nullModelClass::getiPsixMat(arma::mat xMat)   // matrix with dim of n x (J-1)
{
  arma::mat iPsi_xMat(n, J-1);
  for(int i = 0; i < n; i++){   // loop for samples
    double sumx = sum(xMat.row(i));
    for(int j = 0; j < J-1; j++){
      iPsi_xMat(i,j) = sumx / muMat(i,J-1) + xMat(i,j) / muMat(i,j);
    }
  }
  return(iPsi_xMat);
}

// update parameters (except tau) until converge
void nullModelClass::updateParaConv(string excludechr)
{
  int iter = 0;
  while(iter < maxiter)
  {
    beta0 = beta;
    // update beta and bVec
    updatePara(excludechr);
    updateMats();
    
    // update eps (cutpoints)
    updateEps();
    updateMats();
    
    cout << "beta: " << endl << beta << endl;
    iter++;
    
    diffBeta = max(abs(beta - beta0)/(abs(beta) + abs(beta0) + tolBeta));
    if(diffBeta < tolBeta)
      break;
  }
}

void nullModelClass::updateEps()
{
  int iter = 0;
  while(iter < maxiterEps)
  {
    arma::vec eps0 = eps;
    updateEpsOneStep();
    updateMats();
    
    diffeps = max(abs(eps - eps0));
    if(diffeps < tolEps) break;
    iter++;
  }
  cout << "UpdateEps iter: " << iter << endl;
  cout << "eps: " << endl << eps << endl;
}

// need update in case that eps(k+1) < eps(k)
void nullModelClass::updateEpsOneStep()
{
  // the first eps is fixed at 0
  arma::vec d1eps(J-2, arma::fill::zeros);
  arma::mat d2eps(J-2, J-2, arma::fill::zeros);
  double temp1, temp2, temp3;
  
  for(int k = 1; k < J-1; k++){
    for(int i = 0; i < n; i++){
      temp1 = yMat(i,k)/muMat(i,k) - yMat(i,k+1)/muMat(i,k+1);
      temp2 = - yMat(i,k)/muMat(i,k)/muMat(i,k) - yMat(i,k+1)/muMat(i,k+1)/muMat(i,k+1);
      d1eps(k-1) += WMat(i,k) * temp1;
      d2eps(k-1, k-1) += WMat(i,k) * (1-2*nuMat(i,k)) * temp1 + WMat(i,k)*WMat(i,k)*temp2;
      if(k < J-2){
        temp3 = WMat(i,k) * WMat(i,k+1) * yMat(i,k+1)/muMat(i,k+1)/muMat(i,k+1);
        d2eps(k-1, k) += temp3;
        d2eps(k, k-1) += temp3;
      }
    }
  }
  
  arma::vec deps = -1 * inv(d2eps) * d1eps;
  for(int k = 1; k < J-1; k++){
    eps(k) += deps(k-1);
  }
}

void nullModelClass::updatePara(string excludechr)
{
  getPCGofSigmaAndCovaMat(CovaMat, iSigma_CovaMat, excludechr);
  arma::vec YVec = convert1(YMat);
  getPCGofSigmaAndVector(YVec, iSigma_YVec, excludechr); 
  
  // update beta
  arma::mat Cova_iSigma_CovaMat = CovaMat.t() * iSigma_CovaMat;
  arma::vec Cova_iSigma_YVec = CovaMat.t() * iSigma_YVec;
  beta = inv(Cova_iSigma_CovaMat) * Cova_iSigma_YVec;
  
  // update bVec
  arma::vec Z_iSigma_YVec = ZMat(iSigma_YVec);
  arma::vec Z_iSigma_Xbeta = ZMat(iSigma_CovaMat * beta);
  bVec = tau * getKinbVec(Z_iSigma_YVec - Z_iSigma_Xbeta, ptrGeno, excludechr);
}

// use PCG to calculate iSigma_xMat = Sigma^-1 %*% xMat
void nullModelClass::getPCGofSigmaAndCovaMat(arma::mat xMat,             // matrix with dim of n(J-1) x p
                                             arma::mat& iSigma_xMat,      // matrix with dim of n(J-1) x p
                                             string excludechr)
{
  int p1 = xMat.n_cols;
  for(int i = 0; i < p1; i++){
    arma::vec y1Vec = xMat.col(i);
    arma::vec iSigma_y1Vec = iSigma_xMat.col(i);
    getPCGofSigmaAndVector(y1Vec, iSigma_y1Vec, excludechr);
    iSigma_xMat.col(i) = iSigma_y1Vec;
  }
}

// use PCG to calculate xVec = Sigma^-1 %*% yVec
void nullModelClass::getPCGofSigmaAndVector(arma::vec y1Vec,    // vector with length of n(J-1)
                                            arma::vec& xVec,   // vector with length of n(J-1)
                                            string excludechr)
{
  arma::mat xMat = convert2(xVec);
  arma::mat y1Mat = convert2(y1Vec);
  // r2Vec and z2Vec are for the current step; r1Vec and z1Vec are for the previous step
  int iter = 0;
  arma::mat r2Mat = y1Mat - getSigmaxMat(xMat, excludechr);  // n x (J-1): r0 = y1Mat- Sigma %*% xMat
  double meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(n*(J-1));
  if(meanL2 <= tolPCG){
    // do nothing, xMat is already close to (Sigma)^-1 %*% y1Mat
  }else{
    iter++;
    
    arma::cube InvBlockDiagSigma = getInvBlockDiagSigma();
    arma::mat z2Mat = solverBlockDiagSigma(InvBlockDiagSigma, r2Mat);
    
    //
    arma::mat z1Mat, r1Mat;
    double beta1 = 0;
    arma::mat pMat = z2Mat;
    arma::mat ApMat = getSigmaxMat(pMat, excludechr);
    double alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
    xMat = xMat + alpha * pMat;
    r1Mat = r2Mat;
    z1Mat = z2Mat;
    r2Mat = r1Mat - alpha * ApMat;
    
    meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(n*(J-1));
    
    while (meanL2 > tolPCG && iter < maxiterPCG){
      iter++;
      
      //  z2Mat = minvMat % r2Mat;
      z2Mat = solverBlockDiagSigma(InvBlockDiagSigma, r2Mat);
      //
      beta1 = getInnerProd(z2Mat, r2Mat) / getInnerProd(z1Mat, r1Mat);
      pMat = z2Mat + beta1 * pMat;
      ApMat = getSigmaxMat(pMat, excludechr);
      alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
      xMat = xMat + alpha * pMat;
      r1Mat = r2Mat;
      z1Mat = z2Mat;
      r2Mat = r1Mat - alpha * ApMat;
      meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(n*(J-1));
    }
  }
  
  xVec = convert1(xMat);
  if (iter >= maxiterPCG){
    cout << "pcg did not converge. You may increase maxiter number." << endl;
  }
  cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
}

arma::mat nullModelClass::solverBlockDiagSigma(arma::cube& InvBlockDiagSigma,      // (J-1) x (J-1) x n
                                               arma::mat& xMat)                 // n x (J-1)
{
  arma::mat outMat(n, J-1);
  for(int i = 0; i < n; i++){
    outMat.row(i) = xMat.row(i) * InvBlockDiagSigma.slice(i); // could invert matrix?? be careful!
  }
  return(outMat);
}


// used in getPCGofSigmaAndVector()
arma::cube nullModelClass::getInvBlockDiagSigma()
{
  // get diagonal elements of GRM
  arma::vec* pDiagStdGeno = ptrGeno->getDiagStdGeno();
  arma::vec DiagGRM = tau * (*pDiagStdGeno) / M;   // n x 1
  
  double temp;
  //
  arma::cube InvBlockDiagSigma(J-1, J-1, n, arma::fill::zeros);
  for(int i = 0; i < n; i++){
    for(int j2 = 0; j2 < J-1; j2++){
      for(int j1 = 0; j1 < J-1; j1++){
        temp = iRMat(i,j2) * (1 / muMat(i,J-1)) * iRMat(i,j1) + DiagGRM(i);
        if(j2 == j1){
          temp += iRMat(i,j2) * (1 / muMat(i,j2)) * iRMat(i,j1); 
        }
        InvBlockDiagSigma(j2, j1, i) = temp;
      }
    }
    InvBlockDiagSigma.slice(i) = inv(InvBlockDiagSigma.slice(i));
  }
  return(InvBlockDiagSigma);
}

// convert from n x (J-1) matrix to n(J-1) vector
arma::vec nullModelClass::convert1(arma::mat xMat) // matrix with dim of n x (J-1)
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

// convert n(J-1) vector to n x (J-1) matrix
arma::mat nullModelClass::convert2(arma::vec xVec) 
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

// yMat1 = Sigma %*% xMat
arma::mat nullModelClass::getSigmaxMat(arma::mat xMat,   // matrix with dim of n x (J-1)
                                       string excludechr)
{
  arma::mat iR_xMat = iRMat % xMat;
  arma::mat iPsi_iR_xMat = getiPsixMat(iR_xMat);
  arma::mat yMat1 = iRMat % iPsi_iR_xMat;
  if(tau == 0){}
  else{
    arma::vec tZ_xMat = getRowSums(xMat);  // rowSums(xMat): n x 1
    arma::vec V_tZ_xMat = getKinbVec(tZ_xMat, ptrGeno, excludechr);
    yMat1.each_col() += tau * V_tZ_xMat;
  }
  return(yMat1);
}

arma::vec nullModelClass::getRowSums(arma::mat xMat)
{
  int n1 = xMat.n_rows;
  int n2 = xMat.n_cols;
  arma::vec y1Vec(n1, arma::fill::zeros);
  for(int i = 0; i < n1; i++){
    for(int j = 0; j < n2; j++){
      y1Vec(i) += xMat(i,j);
    }
  }
  return(y1Vec);
}


Rcpp::List nullModelClass::getNullModel()
{
  Rcpp::List outList = List::create(Named("N")=n,              // number of samples
                                    Named("M")=M,              // number of SNPs in Plink file
                                    Named("controlList")=controlList,
                                    Named("eta")=eta,          // X %*% beta + bVec
                                    Named("yVec")=yVec,        // matrix with dim of n x J: observation
                                    Named("Cova")=Cova,        // matrix with dim of n(J-1) x p: covariates
                                    Named("muMat")=muMat,      // matrix with dim of n x J: probability
                                    Named("YMat")=YMat,        // matrix with dim of n x (J-1): working variables
                                    Named("beta")=beta,        // parameter for covariates
                                    Named("bVec")=bVec,        // terms of random effect 
                                    Named("tau")=tau,          // variance component
                                    Named("eps")=eps,          // cutpoints
                                    Named("LOCOList")=LOCOList);         
  
  return(outList);
}


////////////////// main function //////////////

// [[Rcpp::export]]
Rcpp::List fitNullcpp(std::string Plink,
                      arma::vec posSampleInPlink,
                      arma::mat CovaR,
                      arma::Col<int> yVecR, // should be from 1 to J
                      arma::vec betaR,
                      arma::vec bVecR,
                      arma::vec epsR,
                      double tauR,
                      arma::mat GMatRatioR,
                      Rcpp::List controlListR)
{
  // read in genotype from Plink file for GRM
  genoClass obj;
  float memoryChunk = controlListR["memoryChunk"];
  obj.setGenoObj(Plink, posSampleInPlink, memoryChunk);
  genoClass* ptrGenoInput = &obj;
  
  // Null model fitting
  nullModelClass objm;
  objm.setNullModel(CovaR, yVecR, ptrGenoInput, betaR,  bVecR,  epsR,  tauR, GMatRatioR, controlListR);
  objm.fitNullModel();
  
  Rcpp::List outList = objm.getNullModel();
  return(outList);
}


