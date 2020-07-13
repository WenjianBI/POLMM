
#ifndef PLINK_H
#define PLINK_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace Plink {

using namespace Rcpp;
using namespace std;

class PlinkClass{
private:
  
  // PLINK format
  const static unsigned char HOM_REF = 0x3;  // 0b11 ;
  const static unsigned char HET = 0x2;      // 0b10 ;
  const static unsigned char HOM_ALT = 0x0;  // 0b00 ;
  const static unsigned char MISSING = 0x1;  // 0b01 ;
  
  // PLINK files
  string m_bimfile, m_famfile, m_bedfile;
  arma::ivec m_posSampleInPlink;
  
  // key parameters of sample size, number of markers, et al.
  int m_N0, m_N, m_M0;
  long long int m_numBytesofEachMarker0;
  long long int m_numBytesofEachMarker;
  Rcpp::StringVector m_chrVec;
  
  // input file stream of .bed file
  ifstream m_ibedfile;
  
  // pipeline: OneMarkerG4 --> bufferG4 --> bufferG1 --> OneMarkerG1
  std::vector<unsigned char> m_OneMarkerG4;
  
  // setup m_chrVec from m_bimfile
  Rcpp::StringVector setChrVec();
  
  // extract geno (0,1,2,3) at specific pos (0,1,2,3) of address c (1 byte)  
  void getGenotype(unsigned char* c, const int pos, int& geno) {
    geno = ((*c) >> (pos << 1)) & 0b11; 
  }
  
public:
  
  // setup PlinkClass
  void setPlinkObj(string t_bimfile,
                   string t_famfile,
                   string t_bedfile,
                   arma::ivec t_posSampleInPlink);
  
  // get genotype of one marker (and frequency and missing rate)
  arma::vec getOneMarker(long long int t_posMarker, 
                         double& t_freq, 
                         double& t_missingRate);
  
  // impute genotype of one marker based on its frequency
  void imputeOneMarker(arma::vec& t_OneMarkerG1, 
                       double t_freq);
  
  // get GMat
  arma::mat getGMat(int t_nMarker, 
                    string t_chrName, 
                    double t_minMafVarRatio, 
                    double t_maxMissingVarRatio);
  
  // get chrVec
  Rcpp::StringVector getChrVec(){return m_chrVec;}
  int getN0(){return m_N0;}
  int getN(){return m_N;}
  int getM0(){return m_M0;}
  long long int getnumBytesofEachMarker0(){return m_numBytesofEachMarker0;}
  long long int getnumBytesofEachMarker(){return m_numBytesofEachMarker;}
  
};

}

#endif
