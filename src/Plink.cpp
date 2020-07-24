
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <RcppArmadilloExtensions/sample.h> // sample
#include <string>

#include "Plink.hpp"
#include "SubFunc.hpp"

namespace Plink {

using namespace Rcpp;
using namespace std;

// set PlinkClass by reading plink files
void PlinkClass::setPlinkObj(string t_bimfile,
                             string t_famfile,
                             string t_bedfile,
                             arma::ivec t_posSampleInPlink)
{
  m_bimfile = t_bimfile;
  m_famfile = t_famfile;
  m_bedfile = t_bedfile;
  
  m_N0 = countLine(m_famfile);
  m_M0 = countLine(m_bimfile);
  m_N = t_posSampleInPlink.size();
  
  m_numBytesofEachMarker0 = (m_N0 + 3) / 4;
  m_numBytesofEachMarker = (m_N + 3) / 4;
  
  m_posSampleInPlink = t_posSampleInPlink;
  m_chrVec = setChrVec();
  
  if((max(m_posSampleInPlink) > m_N0) || (min(m_posSampleInPlink) < 1)) 
    stop("posSampleInPlink should be between 1 and N");
  
  m_ibedfile.open(m_bedfile.c_str(), ios::binary);
  
  // get binray data from plink files
  m_OneMarkerG4.reserve(m_numBytesofEachMarker0);  
  m_OneMarkerG4.resize(m_numBytesofEachMarker0);
}

Rcpp::StringVector PlinkClass::setChrVec()
{
  ifstream ifile(m_bimfile);
  
  Rcpp::StringVector chrVec(m_M0);
  string temp;
  
  int counts = 0;
  while(getline(ifile, temp)){
    string chr = temp.substr(0, temp.find('\t'));
    chrVec[counts] = Rcpp::String(chr);
    counts ++;
  }
  
  return chrVec;
}

arma::mat PlinkClass::getGMat(int t_nMarker, 
                              string t_chrName, 
                              double t_minMafVarRatio, 
                              double t_maxMissingVarRatio)
{
  arma::mat GMat(m_N, t_nMarker);
  
  std::vector<int> indexSNPs(m_M0);
  std::iota (std::begin(indexSNPs), std::end(indexSNPs), 0);
  
  if(t_chrName != "none"){
    indexSNPs = whichCPP(m_chrVec, t_chrName);
  }
  
  indexSNPs = Rcpp::RcppArmadillo::sample(indexSNPs, indexSNPs.size(), FALSE);
  
  cout << "There are " << indexSNPs.size() << " markers in Plink files." << endl;
  
  double freq, missingRate;
  int posGMat = 0;
  unsigned int i;
  for(i = 0; i < indexSNPs.size(); i ++){
    arma::vec oneMarker = getOneMarker(indexSNPs[i], freq, missingRate);
    
    cout << "freq is " << freq << " and missingRate is "<< missingRate << "." << endl << endl;
    cout << "t_minMafVarRatio is " << t_minMafVarRatio << " and t_maxMissingVarRatio is "<< t_maxMissingVarRatio << "." << endl << endl;
    
    if(freq >= t_minMafVarRatio && freq <= 1 - t_minMafVarRatio && missingRate <= t_maxMissingVarRatio){
      
      // long long int posMarker = indexSNPs[i];
      // cout << "extract " << posMarker << "-th marker at chr " << m_chrVec[posMarker] << "." << endl;
      // cout << "freq is " << freq << " and missingRate is "<< missingRate << "." << endl << endl;
      
      imputeOneMarker(oneMarker, freq);
      GMat.col(posGMat) = oneMarker;
      posGMat ++;
    }
    if(posGMat >= t_nMarker){
      break;
    }
  }
  
  if(i == indexSNPs.size()){
    cout << "Only extract " << posGMat << " markers!" << endl;
    stop("Probably give an incorrect chromosome name!");
  }
  
  return GMat;
}

arma::vec PlinkClass::getOneMarker(long long int t_posMarker, 
                                   double& t_freq, 
                                   double& t_missingRate)
{
  int sum = 0;
  int numMissing = 0;
  arma::vec OneMarkerG1(m_N);
  
  m_ibedfile.seekg(3 + m_numBytesofEachMarker0 * t_posMarker);
  m_ibedfile.read((char*)(&m_OneMarkerG4[0]), m_numBytesofEachMarker0);
  
  for(int index = 0; index < m_N; index ++){
    int ind = m_posSampleInPlink[index] - 1;          // C++ start from 0
    unsigned char bufferG4 = m_OneMarkerG4[ind/4];    // unsigned char: 1 byte for 4 genotypes (4 samples)
    int bufferG1;                                     // int: 1 genotype (1 sample)
    getGenotype(&bufferG4, ind%4, bufferG1);          // bufferG4 -> bufferG1
    switch(bufferG1){
    case HOM_REF: break;
    case HET: sum+=1; break;
    case HOM_ALT: sum+=2; break;
    case MISSING: numMissing++; break;
    }
    OneMarkerG1[index] = bufferG1;
  }
  
  int count = m_N - numMissing;
  t_missingRate = (double)numMissing / (double)m_N;
  t_freq = (double)sum / (double)count;
  t_freq = t_freq / 2;
  return OneMarkerG1;
}

void PlinkClass::imputeOneMarker(arma::vec& t_OneMarkerG1,
                                 double t_freq)
{
  for(int index = 0; index < m_N; index ++){
    int bufferG1 = t_OneMarkerG1[index];
    double newG1 = 0;
    switch(bufferG1){
    case HOM_REF: break;
    case HET: newG1 = 1; break;
    case HOM_ALT: newG1 = 2; break;
    case MISSING: newG1 = 2 * t_freq; break;
    }
    t_OneMarkerG1[index] = newG1;
  }
}

}
