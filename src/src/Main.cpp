
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "PLINK.hpp"
#include "POLMM.hpp"
#include "util.hpp"
#include "Main.hpp"

// need to pre-define "ptr_gPLINKobj" and "ptr_gPOLMMobj"

static PLINK::PlinkClass* ptr_gPLINKobj = NULL;
static POLMM::POLMMClass* ptr_gPOLMMobj = NULL;

// [[Rcpp::export]]
void setPLINKobjInR(std::string t_bimFile,
                    std::string t_famFile,
                    std::string t_bedFile,
                    std::vector<std::string> t_SampleInModel)
{
  ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
                                        t_famFile,
                                        t_bedFile,
                                        t_SampleInModel);
  
  int n = ptr_gPLINKobj->getN();
  std::cout << "n:\t" << n << std::endl;
}

// [[Rcpp::export]]
void setPOLMMobjInR(arma::mat t_muMat,
                    arma::mat t_iRMat,
                    arma::mat t_Cova,
                    arma::vec t_yVec,
                    Rcpp::List t_SPmatR,    // output of makeSPmatR()
                    double t_tau,
                    bool t_printPCGInfo,
                    double t_tolPCG,
                    int t_maxiterPCG)
{
  arma::umat locations = t_SPmatR["locations"];
  arma::vec values = t_SPmatR["values"];
  arma::sp_mat SparseGRM = arma::sp_mat(locations, values);
  ptr_gPOLMMobj = new POLMM::POLMMClass(t_muMat,
                                        t_iRMat,
                                        t_Cova,
                                        t_yVec,
                                        SparseGRM,
                                        t_tau,
                                        t_printPCGInfo,
                                        t_tolPCG,
                                        t_maxiterPCG);
}

// [[Rcpp::export]]
Rcpp::List MAIN_REGION(std::vector<std::string> t_MarkerReqstd,
                       double t_NonZero_cutoff,
                       double t_StdStat_cutoff,
                       int t_maxMarkers,
                       std::string t_outputFile,
                       double t_missingRate_cutoff,
                       double t_maxMAF_cutoff)
{
  // extract information from global variable ptr_gPLINKobj
  int n = ptr_gPLINKobj->getN();
  std::vector<uint32_t> posMarkerInPlink = ptr_gPLINKobj->getPosMarkerInPlink(t_MarkerReqstd);
  int q = posMarkerInPlink.size();         // number of markers in the region
  
  // set up output
  std::vector<std::string> a1Vec; 
  std::vector<std::string> a2Vec; 
  std::vector<std::string> markerVec;
  std::vector<uint32_t> pdVec;
  std::vector<double> freqVec; 
  std::vector<double> StatVec;
  std::vector<bool> flipVec;
  arma::fmat VarSMat;
  
  //
  arma::fmat adjGMat(n, t_maxMarkers);       // adjusted genotype vector 
  arma::fmat ZPZ_adjGMat(n, t_maxMarkers);   // t(Z) %*% P %*% Z %*% adjGMat
  
  int indexPassingQC = 0;
  int indexChunkSave = 0;
  
  // loop for all markers
  for(int i = 0; i < q; i++){

    uint32_t posMarker = posMarkerInPlink.at(i);
    double freq, missingRate;
    std::vector<uint32_t> posMissingGeno;
    std::string a1, a2, marker;
    uint32_t pd;
    uint8_t chr;
    bool flip = false;
    bool flagTrueGeno = true;

    arma::vec GVec = ptr_gPLINKobj->getOneMarker(posMarker, freq, missingRate, posMissingGeno,
                                                 a1, a2, marker, pd, chr, flagTrueGeno);
    double MAF = std::min(freq, 1 - freq);

    // Quality Control (QC) based on missing rate and allele frequency
    if(missingRate > t_missingRate_cutoff || MAF > t_maxMAF_cutoff)
      continue;

    // push back to output
    a1Vec.push_back(a1);
    a2Vec.push_back(a2);
    markerVec.push_back(marker);
    pdVec.push_back(pd);

    if(missingRate != 0)
      imputeGeno(GVec, freq, posMissingGeno);

    if(freq > 0.5){
      GVec = 2 - GVec;
      flip = true;
    }

    freqVec.push_back(MAF);
    flipVec.push_back(flip);

    arma::vec adjGVec = ptr_gPOLMMobj->getadjGFast(GVec);
    
    double Stat = ptr_gPOLMMobj->getStatFast(adjGVec);
    
    StatVec.push_back(Stat);

    // get t(Z) %*% P %*% Z %*% adjGVec for each marker
    arma::vec ZPZ_adjGVec = ptr_gPOLMMobj->get_ZPZ_adjGVec(adjGVec);
    
    double VarS = as_scalar(adjGVec.t() * ZPZ_adjGVec);

    double StdStat = std::abs(Stat) / sqrt(VarS);

    if(StdStat > t_StdStat_cutoff){ // then use SPA/ER to correct p value
      // functions of SPA or ER
    }

    // insert adjGVec and ZPZ_adjGVec into a pre-defined matrix
    adjGMat.col(indexPassingQC) = arma::conv_to<arma::fvec>::from(adjGVec);
    ZPZ_adjGMat.col(indexPassingQC) = arma::conv_to<arma::fvec>::from(ZPZ_adjGVec);

    indexPassingQC++;
    
    if(indexPassingQC % t_maxMarkers == 0){
      adjGMat.save(t_outputFile + "_adjGMat" + std::to_string(indexChunkSave) + ".bin");
      ZPZ_adjGMat.save(t_outputFile + "_ZPZ_adjGMat" + std::to_string(indexChunkSave) + ".bin");
      indexPassingQC = 0;
      indexChunkSave++;
    }
  }

  // the region includes more markers than memory requested
  int LastChunkSave = 0;
  if((indexChunkSave > 0) & (indexPassingQC != 0)){ 
    adjGMat = adjGMat.cols(0, indexPassingQC - 1);
    ZPZ_adjGMat = ZPZ_adjGMat.cols(0, indexPassingQC - 1);
    adjGMat.save(t_outputFile + "_adjGMat" + std::to_string(indexChunkSave) + ".bin");
    ZPZ_adjGMat.save(t_outputFile + "_ZPZ_adjGMat" + std::to_string(indexChunkSave) + ".bin");
    indexChunkSave++;
    LastChunkSave = 1;
  }
  
  indexPassingQC = indexPassingQC + (indexChunkSave - LastChunkSave) * t_maxMarkers;
  
  // calculate variance-covariance matrix
  VarSMat.resize(indexPassingQC, indexPassingQC);    // variance matrix (after adjusting for relatedness)
  
  // not so many markers in the region, so everything is in memory
  if(indexChunkSave == 0)
  {
    VarSMat = getSymmMat(adjGMat, ZPZ_adjGMat, indexPassingQC);
  }
  
  // the region includes more markers than limitation, so everything is in hard-drive storage
  if(indexChunkSave > 0)
  {
    int first_row = 0;
    int first_col = 0;
    int last_row = t_maxMarkers - 1;
    int last_col = t_maxMarkers - 1;
    
    for(int index1 = 0; index1 < indexChunkSave; index1++)
    {
      adjGMat.load(t_outputFile + "_adjGMat" + std::to_string(index1) + ".bin");
      
      // off-diagonal sub-matrix
      for(int index2 = 0; index2 < index1; index2++)
      {
        ZPZ_adjGMat.load(t_outputFile + "_ZPZ_adjGMat" + std::to_string(index2) + ".bin");

        arma::fmat offVarSMat = getfmatMulti(adjGMat, ZPZ_adjGMat);
        VarSMat.submat(first_row, first_col, std::min(last_row, indexPassingQC - 1), std::min(last_col, indexPassingQC - 1)) = offVarSMat;
        VarSMat.submat(first_col, first_row, std::min(last_col, indexPassingQC - 1), std::min(last_row, indexPassingQC - 1)) = offVarSMat.t();
        first_col += t_maxMarkers;
        last_col += t_maxMarkers;
      }

      // // diagonal sub-matrix
      ZPZ_adjGMat.load(t_outputFile + "_ZPZ_adjGMat" + std::to_string(index1) + ".bin");
      arma::fmat diagVarSMat = getSymmMat(adjGMat, ZPZ_adjGMat, adjGMat.n_cols);
      VarSMat.submat(first_row, first_col, std::min(last_row, indexPassingQC - 1), std::min(last_col, indexPassingQC - 1)) = diagVarSMat;
      first_row += t_maxMarkers;
      last_row += t_maxMarkers;
      first_col = 0;
      last_col = t_maxMarkers - 1;
    }
  }
  
  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("StatVec") = StatVec,
                                          Rcpp::Named("VarSMat") = VarSMat,
                                          Rcpp::Named("markerVec") = markerVec,
                                          Rcpp::Named("freqVec") = freqVec,
                                          Rcpp::Named("a1Vec") = a1Vec,
                                          Rcpp::Named("a2Vec") = a2Vec,
                                          Rcpp::Named("pdVec") = pdVec,
                                          Rcpp::Named("flipVec") = flipVec);
  return OutList;
}

// yMat = t(xMat1) * xMat2; 
// if yMat is a symmetric matrix, we only need to calculate its lower triangular part
arma::fmat getSymmMat(arma::fmat& xMat1,  // n x p
                      arma::fmat& xMat2,  // n x p
                      int p)
{
  int p1 = xMat1.n_cols;
  if(p > p1){
    Rcpp::stop("check function getSymmMat() in Main.cpp.");
  }
  arma::fmat yMat(p, p);
  for(int i = 0; i < p; i++){
    arma::fvec xVec1 = xMat1.col(i);
    for(int j = 0; j < i; j++){
      float cov = as_scalar(xVec1.t() * xMat2.col(j));
      yMat(i, j) = yMat(j, i) = cov;
    }
    float var = as_scalar(xMat1.col(i).t() * xMat2.col(i));
    yMat(i, i) = var;
  }
  return yMat;
}

// yMat = t(xMat1) * xMat2; 
// even if yMat is a nonsymmetric matrix, if xMat1 and xMat2 is arma::fmat, sometimes the direct matrix multiplication cannot be used
arma::fmat getfmatMulti(arma::fmat& xMat1,  // n x p1
                        arma::fmat& xMat2)  // n x p2
{
  int p1 = xMat1.n_cols;
  int p2 = xMat2.n_cols;
  arma::fmat yMat(p1, p2);
  for(int i = 0; i < p1; i++){
    arma::fvec xVec1 = xMat1.col(i);
    for(int j = 0; j < p2; j++){
      float value = as_scalar(xVec1.t() * xMat2.col(j));
      yMat(i, j) = value;
    }
  }
  return yMat;
}