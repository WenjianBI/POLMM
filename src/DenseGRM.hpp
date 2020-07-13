
#ifndef DENSEGRM_H
#define DENSEGRM_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "Plink.hpp"

namespace DenseGRM {

using namespace Rcpp;
using namespace std;
using namespace Plink;

class DenseGRMClass{
private:
  
  // Basic information
  int m_N0, m_N;                  // number of samples, before and after filter
  int m_M0, m_M;                  // number of markers, before and after filter
  
  long long int m_numBytesofEachMarker0;
  long long int m_numBytesofEachMarker;
  
  // arrays to store genotype data to avoid large continuous memory usage
  std::vector< std::vector<unsigned char>* > m_genoVecofPointers;
  int m_numArrays;                  
  int m_numMarkersofEachArray;
  int m_numMarkersofLastArray;
  double m_numBytesReserve;         // unit: Gb
  
  // OneMarkerG4 --> bufferG4 --> bufferG1 --> OneMarkerG1
  arma::vec m_OneMarkerG1;
  
  // summary information for markers
  arma::vec m_DiagStdGeno;
  Rcpp::NumericVector m_freqVec; 
  Rcpp::NumericVector m_invStdVec;
  Rcpp::StringVector m_chrVec0;
  Rcpp::StringVector m_chrVec;
  Rcpp::List m_chrIndexLOCO;
  
  // functions used in DenseGRMClass
  void setArrays(PlinkClass* t_ptrPlinkObj, double t_memoryChunk);
  void setOneMarkerArray(int t_indexMarker);
  void setDiagStdGeno();
  void setchrIndexLOCO(Rcpp::StringVector t_chrVecNames);

  // insert geno (0,1,2,3) to specific pos (0,1,2,3) of address c (1 byte)
  void setGenotype(unsigned char* c, const int pos, const int geno) {
    (*c) |= (geno << (pos << 1));
  }
  
  void setStdGenoLookUpArr(double maf, double invsd, arma::vec& stdGenoLookUpArr){
    double maf2 = 2 * maf;
    stdGenoLookUpArr(0) = (2-maf2)*invsd; // HOM_ALT = 0x0;
    stdGenoLookUpArr(1) = 0;              // MISSING = 0x1;
    stdGenoLookUpArr(2) = (1-maf2)*invsd; // HET = 0x2;
    stdGenoLookUpArr(3) = (0-maf2)*invsd; // HOM_REF = 0x3;
  }
  
public:
  
  // setup DenseGRMClass
  void setDenseGRMObj(PlinkClass* t_ptrPlinkObj, 
                      double t_memoryChunk,     // unit is Gb
                      double t_minMafGRM, 
                      double t_maxMissingGRM);
  
  void closeDenseGRMObj();
  void getOneMarkerStd(size_t t_indexMarker, arma::vec* t_oneMarkerStd);
  
  int getN(){return m_N;}
  int getM(){return m_M;}
  arma::vec getFreqVec(){return m_freqVec;}
  arma::vec* getDiagStdGeno(){return &m_DiagStdGeno;}
  Rcpp::StringVector getChrVec(){return m_chrVec;}
  Rcpp::List getChrIndexLOCO(){return m_chrIndexLOCO;}
};

arma::vec getKinbVec(arma::vec t_bVec, DenseGRMClass* t_ptrDenseGRM, string t_excludeChr, int t_grainSize);

}

#endif
