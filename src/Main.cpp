
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <string>

#include "DenseGRM.hpp"
#include "Plink.hpp"
#include "SubFunc.hpp"
#include "POLMM.hpp"
#include "POLMM_GENE.hpp"

using namespace Rcpp;
using namespace std;
using namespace Plink;
using namespace DenseGRM;
using namespace POLMM;
using namespace POLMMGENE;

Rcpp::List getKinMatList(Rcpp::List KinMatListR)
{
  int nKin = KinMatListR.size();
  Rcpp::CharacterVector NameKin = KinMatListR.names();
  Rcpp::List KinMatList_sp;
  for(int i = 0; i < nKin; i ++){
    string excludeChr = string(NameKin[i]);
    Rcpp::List KinMatTemp = KinMatListR[excludeChr];
    arma::umat locations = KinMatTemp["locations"];
    arma::vec values = KinMatTemp["values"];
    int n = KinMatTemp["nSubj"];
    // make a sparse matrix
    arma::sp_mat KinMat(locations, values, n, n);
    KinMatList_sp[excludeChr] = KinMat;
  }
  return KinMatList_sp;
}


// [[Rcpp::export]]
Rcpp::List fitPOLMMcpp(bool t_flagSparseGRM,       // if 1, then use SparseGRM, otherwise, use DenseGRM
                       bool t_flagGMatRatio,       // if 1, then use GMatRatio, otherwise, extract from Plink files
                       std::string t_bimfile,
                       std::string t_famfile,
                       std::string t_bedfile, 
                       arma::ivec t_posSampleInPlink,
                       arma::mat t_Cova,
                       arma::Col<int> t_yVec,     // should be from 1 to J
                       arma::vec t_beta,
                       arma::vec t_bVec,
                       arma::vec t_eps,           // 
                       double t_tau,
                       arma::mat t_GMatRatio,     // only used if m_LOCO = FALSE
                       Rcpp::List t_SparseGRM,
                       Rcpp::List t_controlList)
{
  // Plink and DenseGRM class object
  PlinkClass PlinkObj;
  DenseGRMClass DenseGRMObj;
  
  Rcpp::List KinMatList;
  if(!t_flagGMatRatio || !t_flagSparseGRM){
    PlinkObj.setPlinkObj(t_bimfile, t_famfile, t_bedfile, t_posSampleInPlink);
  }
  
  if(!t_flagSparseGRM){
    double memoryChunk = t_controlList["memoryChunk"];
    double minMafGRM = t_controlList["minMafGRM"];
    double maxMissingGRM = t_controlList["maxMissingGRM"];
    DenseGRMObj.setDenseGRMObj(&PlinkObj, memoryChunk, minMafGRM, maxMissingGRM);
  }else{
    KinMatList = getKinMatList(t_SparseGRM);
  }
  
  // POLMM class object
  POLMMClass POLMMObj;
  POLMMObj.setPOLMMObj(t_flagSparseGRM,       // if 1, then use SparseGRM, otherwise, use DenseGRM
                       t_flagGMatRatio,       // if 1, then use GMatRatio, otherwise, extract from Plink files
                       &PlinkObj,
                       &DenseGRMObj,
                       t_Cova,
                       t_yVec,     // should be from 1 to J
                       t_beta,
                       t_bVec,
                       t_eps,           // 
                       t_tau,
                       t_GMatRatio,     // only used if m_LOCO = FALSE
                       KinMatList,
                       t_controlList);
  
  // Null model fitting
  if(!t_controlList["onlyCheckTime"])
    POLMMObj.fitPOLMM();
  
  Rcpp::List outList = POLMMObj.getPOLMM();
  
  POLMMObj.closeGenoObj();
  
  return(outList);
}

// make a global variable for region- and Gene-based association test
static POLMMGENEClass* ptr_POLMMGENEobj = NULL;

// [[Rcpp::export]]
void setPOLMMGENEobj(int t_maxiterPCG,
                     double t_tolPCG,
                     arma::mat t_Cova,
                     arma::uvec t_yVec,     // should be from 1 to J
                     double t_tau,
                     Rcpp::List t_SparseGRM,    // results of function getKinMatList()
                     Rcpp::List t_LOCOList,
                     arma::vec t_eta,
                     int t_nMaxNonZero)
{
  Rcpp::List KinMatList = getKinMatList(t_SparseGRM);
  
  ptr_POLMMGENEobj = new POLMMGENEClass(t_maxiterPCG, 
                                        t_tolPCG,
                                        t_Cova,
                                        t_yVec,
                                        t_tau,
                                        KinMatList,
                                        t_LOCOList,
                                        t_eta,
                                        t_nMaxNonZero);
}

// [[Rcpp::export]]
void closePOLMMGENEobj()
{
  delete ptr_POLMMGENEobj;
}

// [[Rcpp::export]]
void setPOLMMGENEchr(Rcpp::List t_LOCOList, std::string t_excludechr)
{
  ptr_POLMMGENEobj->setPOLMMGENEchr(t_LOCOList, t_excludechr);
}

// [[Rcpp::export]]
Rcpp::List getStatVarS(arma::mat t_GMat,
                       double t_NonZero_cutoff,
                       double t_StdStat_cutoff)
{
  Rcpp::List OutList = ptr_POLMMGENEobj->getStatVarS(t_GMat,
                                                     t_NonZero_cutoff,
                                                     t_StdStat_cutoff);
  return(OutList);
}

// [[Rcpp::export]]
double getPvalERtoR(arma::vec t_GVec)
{
  double PvalER = ptr_POLMMGENEobj->getPvalERinClass(t_GVec);
  return PvalER;
}



