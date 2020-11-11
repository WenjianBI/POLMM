
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
  std::cout << "Setting Sparse GRM...." << std::endl;
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
  arma::mat VarSMat;
  
  //
  arma::mat adjGMat(t_maxMarkers, n);       // adjusted genotype vector 
  arma::mat ZPZ_adjGMat(n, t_maxMarkers);   // t(Z) %*% P %*% Z %*% adjGMat
  
  int indexPassingQC = 0;
  int indexChunkSave = 0;
  
  // loop for all markers
  for(int i = 0; i < q; i++){
    
    if(i % 1000 == 0){
      std::cout << "Completed " << i << "/" << q << " markers in the region." << std::endl;
      std::cout << "indexPassingQC:\t" << indexPassingQC << std::endl;
    }

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
    if((missingRate > t_missingRate_cutoff) || (MAF > t_maxMAF_cutoff) || (MAF == 0))
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
    adjGMat.row(indexPassingQC) = adjGVec.t();
    ZPZ_adjGMat.col(indexPassingQC) = ZPZ_adjGVec;

    indexPassingQC++;
    
    if(indexPassingQC % t_maxMarkers == 0){
      adjGMat.save(t_outputFile + "_adjGMat" + std::to_string(indexChunkSave) + ".bin");
      ZPZ_adjGMat.save(t_outputFile + "_ZPZ_adjGMat" + std::to_string(indexChunkSave) + ".bin");
      indexPassingQC = 0;
      std::cout << "Completed chunk "<< indexChunkSave << "!" << std::endl;
      indexChunkSave++;
    }
    Rcpp::checkUserInterrupt();
  }

  // total number of markers that pass QC from MAF and missing rate
  int nPassingQC = indexPassingQC + indexChunkSave * t_maxMarkers;
  
  if(indexPassingQC != 0){
    adjGMat = adjGMat.rows(0, indexPassingQC - 1);
    ZPZ_adjGMat = ZPZ_adjGMat.cols(0, indexPassingQC - 1);
  }
  
  // the region includes more markers than memory requested
  if((indexChunkSave > 0) & (indexPassingQC != 0)){ 
    adjGMat.save(t_outputFile + "_adjGMat" + std::to_string(indexChunkSave) + ".bin");
    ZPZ_adjGMat.save(t_outputFile + "_ZPZ_adjGMat" + std::to_string(indexChunkSave) + ".bin");
    std::cout << "Completed chunk "<< indexChunkSave << "!" << std::endl;
    indexChunkSave++;
  }
  
  // calculate variance-covariance matrix
  VarSMat.resize(nPassingQC, nPassingQC);    // variance matrix (after adjusting for relatedness)
  
  // not so many markers in the region, so everything is in memory
  if(indexChunkSave == 0)
  {
    VarSMat = adjGMat * ZPZ_adjGMat;
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
        std::cout << "Analyzing chunks (" << index1 << "/" << indexChunkSave - 1 << ", " << index2 << "/" << indexChunkSave - 1 << ")........" << std::endl;
        ZPZ_adjGMat.load(t_outputFile + "_ZPZ_adjGMat" + std::to_string(index2) + ".bin");

        arma::mat offVarSMat = adjGMat * ZPZ_adjGMat;
        VarSMat.submat(first_row, first_col, 
                       std::min(last_row, nPassingQC - 1), 
                       std::min(last_col, nPassingQC - 1)) = offVarSMat;
        VarSMat.submat(first_col, first_row, 
                       std::min(last_col, nPassingQC - 1), 
                       std::min(last_row, nPassingQC - 1)) = offVarSMat.t();
        first_col += t_maxMarkers;
        last_col += t_maxMarkers;
      }

      // // diagonal sub-matrix
      std::cout << "Analyzing chunks (" << index1 << "/" << indexChunkSave - 1 << ", " << index1 << "/" << indexChunkSave - 1 << ")........" << std::endl;
      ZPZ_adjGMat.load(t_outputFile + "_ZPZ_adjGMat" + std::to_string(index1) + ".bin");
      // arma::fmat diagVarSMat = getSymmMat(adjGMat, ZPZ_adjGMat, adjGMat.n_cols);
      arma::mat diagVarSMat = adjGMat * ZPZ_adjGMat;
      VarSMat.submat(first_row, first_col, 
                     std::min(last_row, nPassingQC - 1), 
                     std::min(last_col, nPassingQC - 1)) = diagVarSMat;
      first_row += t_maxMarkers;
      last_row += t_maxMarkers;
      first_col = 0;
      last_col = t_maxMarkers - 1;
      Rcpp::checkUserInterrupt();
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
