
#ifndef PLINK_HPP
#define PLINK_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace PLINK {

class PlinkClass{
private:
  
  // information from bim file
  uint32_t m_M0, m_M;
  std::map<std::string, uint8_t> m_chrMaps;
  std::vector<uint8_t> m_chr;                   // Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
  std::vector<std::string> m_MarkerInPlink;     // Variant identifier
  std::vector<float> m_gd;                      // Position in morgans or centimorgans (safe to use dummy value of '0')
  std::vector<uint32_t> m_pd;                   // Base-pair coordinate (1-based; limited to 2^31-2)
  std::vector<std::string> m_a1;                // Allele 1 (corresponding to clear bits in .bed; usually minor)
  std::vector<std::string> m_a2;                // Allele 2 (corresponding to set bits in .bed; usually major)
  
  // information from fam file
  std::vector<std::string> m_SampleInPlink;
  uint32_t m_N0, m_N;
  unsigned long long int m_numBytesofEachMarker0, m_numBytesofEachMarker;
  
  // input file stream of .bed file
  std::ifstream m_ibedFile;
  
  // PLINK files
  std::string m_bimFile, m_famFile, m_bedFile;
  std::vector<uint32_t> m_posSampleInPlink;
  
  // PLINK format
  const static unsigned char HOM_REF = 0x3;  // 0b11 ;
  const static unsigned char HET = 0x2;      // 0b10 ;
  const static unsigned char HOM_ALT = 0x0;  // 0b00 ;
  const static unsigned char MISSING = 0x1;  // 0b01 ;
  
  // or use "arma::datum::nan"
  std::map<int8_t, int8_t> m_genoMaps = {{3, 0},{2, 1},{0, 2},{1, -1}};
  
  
  // pipeline: OneMarkerG4 --> bufferG4 --> bufferG1 --> OneMarkerG1
  std::vector<unsigned char> m_OneMarkerG4;
  
  void setChrMaps();
  void readBimFile();
  void readFamFile();
  
  // extract geno (0,1,2,3) at specific pos (0,1,2,3) of address c (1 byte)  
  void getGenotype(unsigned char* c, const int pos, int& geno) {
    geno = ((*c) >> (pos << 1)) & 0b11; 
  }
  
public:
  
  PlinkClass(std::string t_bimFile,
             std::string t_famFile,
             std::string t_bedFile,
             std::vector<std::string> t_SampleInModel);
  
  // setup PlinkClass
  void setPlinkobj(std::string t_bimFile,
                   std::string t_famFile,
                   std::string t_bedFile);
  
  void setPosSampleInPlink(std::vector<std::string> t_SampleInModel);
  std::vector<uint32_t> getPosMarkerInPlink(std::vector<std::string> t_MarkerReqstd);
  
  // get genotype of one marker (and frequency and missing rate)
  arma::vec getOneMarker(unsigned long long int t_posMarker, 
                         double& t_freq, 
                         double& t_missingRate,
                         std::vector<uint32_t>& t_posMissingGeno,
                         std::string& t_a1,
                         std::string& t_a2,
                         std::string& t_marker,
                         uint32_t& t_pd,
                         uint8_t& t_chr,
                         bool t_flagTrueGeno);
  

  uint32_t getN0(){return m_N0;}
  uint32_t getN(){return m_N;}
  uint32_t getM0(){return m_M0;}
  uint32_t getM(){return m_M;}
  uint32_t getnumBytesofEachMarker0(){return m_numBytesofEachMarker0;}
  uint32_t getnumBytesofEachMarker(){return m_numBytesofEachMarker;}
  
};

}

#endif
