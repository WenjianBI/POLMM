// [[Rcpp::depends(BH)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <boost/algorithm/string.hpp>
#include "PLINK.hpp"

// make a global variable for future usage
// static PLINK::PlinkClass* ptr_gPLINKobj = NULL;

// ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
//                                       t_famFile,
//                                       t_bedFile,
//                                       t_SampleInModel);

namespace PLINK {

PlinkClass::PlinkClass(std::string t_bimFile,
                       std::string t_famFile,
                       std::string t_bedFile,
                       std::vector<std::string> t_SampleInModel)
{
  setPlinkobj(t_bimFile, t_famFile, t_bedFile);
  setPosSampleInPlink(t_SampleInModel);
}

void PlinkClass::setChrMaps()
{
  uint8_t last_chr_autosome = 22;
  // Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
  for(uint8_t index = 0; index <= 26; index++){
    m_chrMaps[std::to_string(index)] = index;
  }
  m_chrMaps["X"] = last_chr_autosome + 1;
  m_chrMaps["x"] = last_chr_autosome + 1;
  m_chrMaps["Y"] = last_chr_autosome + 2;
  m_chrMaps["y"] = last_chr_autosome + 2;
  m_chrMaps["XY"] = last_chr_autosome + 3;
  m_chrMaps["xy"] = last_chr_autosome + 3;
  m_chrMaps["MT"] = last_chr_autosome + 4;
  m_chrMaps["mt"] = last_chr_autosome + 4;
}

void PlinkClass::readBimFile()
{
  std::ifstream bim(m_bimFile);
  m_M0 = 0;
  std::string line;
  uint8_t chr_item;
  while(getline(bim, line)){
    m_M0++;
    std::vector<std::string> line_elements;
    boost::split(line_elements, line, boost::is_any_of("\t "));
    boost::replace_all(line_elements[line_elements.size() - 1], "\r", "");
    chr_item = m_chrMaps.at(line_elements[0]);
    m_chr.push_back(chr_item);
    m_MarkerInPlink.push_back(line_elements[1]);
    m_gd.push_back(std::stof(line_elements[2]));
    m_pd.push_back(std::stoi(line_elements[3]));
    std::transform(line_elements[4].begin(), line_elements[4].end(), line_elements[4].begin(), toupper);
    std::transform(line_elements[5].begin(), line_elements[5].end(), line_elements[5].begin(), toupper);
    m_a1.push_back(line_elements[4]);
    m_a2.push_back(line_elements[5]);

  }
  m_M = m_M0;
}

void PlinkClass::readFamFile()
{
  std::ifstream fam(m_famFile);
  m_N0 = 0;
  std::string line;
  while(getline(fam, line)){
    m_N0 ++;
    std::vector<std::string> line_elements;
    boost::split(line_elements, line, boost::is_any_of("\t "));
    // boost::replace_all(line_elements[line_elements.size() - 1], "\r", "");
    m_SampleInPlink.push_back(line_elements[1]);
  }
  m_N = m_N0;
  m_numBytesofEachMarker0 = (m_N0 + 3) / 4;
  m_OneMarkerG4.reserve(m_numBytesofEachMarker0);  
  m_OneMarkerG4.resize(m_numBytesofEachMarker0);
} 
    
// set PlinkClass by reading plink files
void PlinkClass::setPlinkobj(std::string t_bimFile,
                             std::string t_famFile,
                             std::string t_bedFile)
{
  m_bimFile = t_bimFile;
  m_famFile = t_famFile;
  m_bedFile = t_bedFile;
  
  setChrMaps();
  readBimFile();
  readFamFile();
  m_ibedFile.open(m_bedFile.c_str(), std::ios::binary);
  
  m_ibedFile.seekg(2);
  char magicNumber3;
  m_ibedFile.read(&magicNumber3, 1);
  
  if(magicNumber3 != 1)
    Rcpp::stop("The third magic number of the plink bed file is not 00000001. Please use SNP-major plink (plink version >= 1.9) files.");
}

void PlinkClass::setPosSampleInPlink(std::vector<std::string> t_SampleInModel)
{
  m_N = t_SampleInModel.size();
  m_numBytesofEachMarker = (m_N + 3) / 4;
  m_posSampleInPlink.clear();
  for(uint32_t i = 0; i < m_N; i++){
    std::string sample = t_SampleInModel.at(i);
    auto pos = std::find(m_SampleInPlink.begin(), m_SampleInPlink.end(), sample);
    if(pos != m_SampleInPlink.end()){
      m_posSampleInPlink.push_back(pos - m_SampleInPlink.begin());
    }else{
      Rcpp::stop("At least one subject requested is not in Plink file.");
    }
  }
}

std::vector<uint32_t> PlinkClass::getPosMarkerInPlink(std::vector<std::string> t_MarkerReqstd)
{
  int M = t_MarkerReqstd.size();
  std::vector<uint32_t> posMarkerInPlink;
  for(int i = 0; i < M; i++){
    std::string marker = t_MarkerReqstd.at(i);
    auto pos = std::find(m_MarkerInPlink.begin(), m_MarkerInPlink.end(), marker);
    if(pos != m_MarkerInPlink.end()){
      posMarkerInPlink.push_back(pos - m_MarkerInPlink.begin());
    }else{
      Rcpp::warning("Marker %s is not found in plink file.", marker);
    }
  }
  return posMarkerInPlink;
}

arma::vec PlinkClass::getOneMarker(uint32_t t_posMarker, 
                                   double& t_freq, 
                                   double& t_missingRate,
                                   std::vector<uint32_t>& t_posMissingGeno,
                                   std::string& t_a1,
                                   std::string& t_a2,
                                   std::string& t_marker,
                                   uint32_t& t_pd,
                                   uint8_t& t_chr,
                                   bool t_flagTrueGeno)
{
  int sum = 0;
  int numMissing = 0;
  arma::vec OneMarkerG1(m_N);
  
  m_ibedFile.seekg(3 + m_numBytesofEachMarker0 * t_posMarker);
  m_ibedFile.read((char*)(&m_OneMarkerG4[0]), m_numBytesofEachMarker0);
  
  t_posMissingGeno.clear();
  t_a1 = m_a1[t_posMarker];
  t_a2 = m_a2[t_posMarker];
  t_marker = m_MarkerInPlink[t_posMarker];
  t_pd = m_pd[t_posMarker];
  t_chr = m_chr[t_posMarker];
  
  for(uint32_t i = 0; i < m_N; i++){
    uint32_t ind = m_posSampleInPlink[i];             // C++ start from 0
    unsigned char bufferG4 = m_OneMarkerG4[ind/4];    // unsigned char: 1 byte for 4 genotypes (4 samples)
    int bufferG1;                                     // int: 1 genotype (1 sample)
    getGenotype(&bufferG4, ind%4, bufferG1);          // bufferG4 -> bufferG1
    switch(bufferG1){
    case HOM_REF: break;
    case HET: sum+=1; break;
    case HOM_ALT: sum+=2; break;
    case MISSING: numMissing++; break;  
    }
    // OneMarkerG1[i] = bufferG1;
    if(t_flagTrueGeno){
      OneMarkerG1[i] = m_genoMaps[bufferG1];
    }else{
      OneMarkerG1[i] = bufferG1;
    }
  }
  
  int count = m_N - numMissing;
  t_missingRate = (double)numMissing / (double)m_N;
  t_freq = (double)sum / (double)count;
  t_freq = t_freq / 2;
  return OneMarkerG1;
}

}

// make a global variable for future usage
// static PLINK::PlinkClass* ptr_gPLINKobj = NULL;

// ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
//                                       t_famFile,
//                                       t_bedFile,
//                                       t_SampleInModel);


// // [[Rcpp::export]]
// void setPLINKobjInR(std::string t_bimFile,
//                     std::string t_famFile,
//                     std::string t_bedFile,
//                     std::vector<std::string> t_SampleInModel)
// {
//   ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
//                                         t_famFile,
//                                         t_bedFile,
//                                         t_SampleInModel);
//   
//   int n = ptr_gPLINKobj->getN();
//   std::cout << "n:\t" << n << std::endl;
// }



