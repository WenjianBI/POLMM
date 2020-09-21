
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include <RcppArmadillo.h>

#include "DenseGRM.hpp"
#include "Plink.hpp"
#include "SubFunc.hpp"

namespace DenseGRM {

using namespace Rcpp;
using namespace RcppParallel;
using namespace std;
using namespace Plink;

void DenseGRMClass::setDenseGRMObj(PlinkClass* t_ptrPlinkObj, 
                                   double t_memoryChunk,     // unit is Gb
                                   double t_minMafGRM, 
                                   double t_maxMissingGRM)    
{
  setArrays(t_ptrPlinkObj, t_memoryChunk);
  
  cout << "Number of samples in Plink file:\t" << m_N0 << endl;
  cout << "Number of markers in Plink file:\t" << m_M0 << endl;
  cout << "Number of samples in model fitting:\t" << m_N << endl;
  
  cout << "numBytesofEachMarker in Plink file:\t" << m_numBytesofEachMarker0 << endl;
  cout << "numBytesofEachMarker in model fitting:\t" << m_numBytesofEachMarker << endl;
  cout << "numBytesReserve:\t" << m_numBytesReserve << " Gb" << endl << endl;
  
  cout << "numArrays:\t" << m_numArrays << endl;
  cout << "numMarkersofEachArray:\t" << m_numMarkersofEachArray << endl;
  cout << "numMarkersofLastArray:\t" << m_numMarkersofLastArray << endl << endl;
  
  ///////////////////////////////////////// MAIN PART //////////////////////////////
  
  double freq, invStd, missingRate;
  arma::vec t1  = getTime();
  m_M = 0;
  string bimTemp;
  for(long long int i = 0; i < m_M0; i ++)
  {
    // OneMarkerG4 --> OneMarkerG1
    m_OneMarkerG1 = t_ptrPlinkObj->getOneMarker(i, freq, missingRate);
    Rcpp::String chr = m_chrVec0[i];
    
    if(freq >= t_minMafGRM && freq <= 1-t_minMafGRM && missingRate <= t_maxMissingGRM){
      setOneMarkerArray(m_M);
      invStd = getinvStd(freq);
      m_freqVec.push_back(freq);
      m_invStdVec.push_back(invStd);
      m_chrVec.push_back(chr);
      m_M ++;
    }
    if((i+1) % 10000 == 0){
      cout << "Complete\t" << i+1 <<"\tSNPs!!!!" << endl;
      cout << "Allele Frequency:\t" << freq << endl;
    }
  }
  
  arma::vec t2  = getTime();
  printTime(t1, t2, "read Plink files");
  
  cout << endl << "Remove markers with MAF < " << t_minMafGRM << " and missing rate > " << t_maxMissingGRM << endl;
  cout << "Number of markers:\t" << m_M << endl << endl;
  
  Rcpp::IntegerVector chrVecCounts = table(m_chrVec);
  Rcpp::StringVector chrVecNames = chrVecCounts.names();
  setchrIndexLOCO(chrVecNames);
  
  for(int i = 0; i < chrVecCounts.size(); i++){
    cout << "Number of markers in chr " << chrVecNames(i) << ":\t" << chrVecCounts(i) << endl;
  }
  cout << endl;
  
  setDiagStdGeno();
  
}

void DenseGRMClass::setArrays(PlinkClass* t_ptrPlinkObj, double t_memoryChunk)
{
  m_N = t_ptrPlinkObj->getN();
  m_N0 = t_ptrPlinkObj->getN0();
  m_M0 = t_ptrPlinkObj->getM0();
  m_numBytesofEachMarker0 = t_ptrPlinkObj->getnumBytesofEachMarker0();
  m_numBytesofEachMarker = t_ptrPlinkObj->getnumBytesofEachMarker();
  m_chrVec0 = t_ptrPlinkObj->getChrVec();
  
  m_numBytesReserve = (double)(m_numBytesofEachMarker+2) * m_M0 / pow(10.0, 9.0);
  m_numMarkersofEachArray = floor(t_memoryChunk * pow(10.0, 9.0) / m_numBytesofEachMarker);
  if(m_numMarkersofEachArray > m_M0)
    m_numMarkersofEachArray = m_M0;
  m_numArrays = (m_M0 - 1) / m_numMarkersofEachArray + 1;
  m_numMarkersofLastArray = m_M0 - (m_numArrays - 1) * m_numMarkersofEachArray;
  
  m_OneMarkerG1.zeros(m_N);
  
  // set m_genoVecofPointers
  m_genoVecofPointers.resize(m_numArrays);
  for (int i = 0; i < m_numArrays-1 ; i++){
    m_genoVecofPointers[i] = new vector<unsigned char>;
    m_genoVecofPointers[i]->reserve(m_numMarkersofEachArray * m_numBytesofEachMarker);
  }
  m_genoVecofPointers[m_numArrays-1] = new vector<unsigned char>;
  m_genoVecofPointers[m_numArrays-1]->reserve(m_numMarkersofLastArray * m_numBytesofEachMarker);
}


// m_oneMarkerG1 --> m_genoVecofPointers
void DenseGRMClass::setOneMarkerArray(int t_indexMarker)
{
  int whichArray = t_indexMarker / m_numMarkersofEachArray;
  
  unsigned char bufferG4 = 0;
  for(int ind = 0; ind < m_N; ind++){
    int posInByte = ind % 4;
    int bufferG1 = m_OneMarkerG1[ind];
    setGenotype(&bufferG4, posInByte, bufferG1);
    if((posInByte == 3) || (ind == (m_N-1))){
      m_genoVecofPointers[whichArray]->push_back(bufferG4); // avoid large continuous memory usage
      bufferG4 = 0;
    }
  }
}

void DenseGRMClass::getOneMarkerStd(size_t t_indexMarker, arma::vec* t_oneMarkerStd)
{
  // avoid large continuous memory usage
  int whichArray = t_indexMarker / m_numMarkersofEachArray;
  int posMarker = t_indexMarker % m_numMarkersofEachArray;
  
  // set up start byte index and end byte index
  int startBtIdx = m_numBytesofEachMarker * posMarker;
  int endBtIdx = startBtIdx + m_numBytesofEachMarker;
  
  arma::vec stdGenoLookUpArr(4);
  setStdGenoLookUpArr(m_freqVec[t_indexMarker], m_invStdVec[t_indexMarker], stdGenoLookUpArr);
  std::vector<unsigned char>* genoPtr = m_genoVecofPointers[whichArray];
  
  int ind = 0;
  for(int BtIdx = startBtIdx; BtIdx < endBtIdx - 1; BtIdx ++){
    unsigned char bufferG4 = genoPtr->at(BtIdx); // unsigned char: 4 markers
    // for(posInByte = 0; (posInByte < 4) & (ind < N); posInByte++, ind++){ 
    // added on 2020/04/04 to avoid checking (ind < N) on each iteration
    for(unsigned char posInByte = 0; posInByte < 4; posInByte ++, ind ++){
      // getGenotype(unsigned char* c, const int pos, int& geno)
      int bufferG1 = (bufferG4 >> (posInByte << 1)) & 0b11;
      t_oneMarkerStd->at(ind) = stdGenoLookUpArr(bufferG1);
    }
  }
  
  // added on 2020/04/04 to avoid checking (ind < N) on each iteration
  unsigned char bufferG4 = genoPtr->at(endBtIdx - 1); // unsigned char: 4 markers
  for(unsigned char posInByte = 0; (posInByte < 4) & (ind < m_N); posInByte++, ind++){
    int bufferG1 = (bufferG4 >> (posInByte << 1)) & 0b11;
    t_oneMarkerStd->at(ind) = stdGenoLookUpArr(bufferG1);
  }
  
}

void DenseGRMClass::setDiagStdGeno()
{
  m_DiagStdGeno.zeros(m_N);
  arma::vec oneMarkerStd(m_N);
  
  // cycle for m_M markers
  for(int i = 0; i < m_M; i ++){
    getOneMarkerStd(i, &oneMarkerStd);  // write Standard Genotype to OneMarkerG1
    m_DiagStdGeno = m_DiagStdGeno + (oneMarkerStd) % (oneMarkerStd);
  }
  m_DiagStdGeno /= m_M;
}

void DenseGRMClass::closeDenseGRMObj()
{
  for (int i = 0; i < m_numArrays; i++){
    m_genoVecofPointers[i]->clear();	
    delete m_genoVecofPointers[i];
  }
  m_genoVecofPointers.clear();
  cout << "Close the genotype object!\n";
}


arma::Mat<int> makeChrIndex(Rcpp::String t_excludeChr, 
                            Rcpp::StringVector t_chrVec)
{
  int m = t_chrVec.length();
  arma::Mat<int> chrIndex;
  arma::Row<int> newChrIndex(2);
  int indexStart = 0;
  for(int i = 0; i < m; i ++){
    if(t_chrVec[i] == t_excludeChr){
      if(indexStart != i){
        newChrIndex(0) = indexStart;
        newChrIndex(1) = i;
        chrIndex.insert_rows(0, newChrIndex);
      }
      indexStart = i + 1;
    }
  }

  if(indexStart != m){
    newChrIndex(0) = indexStart;
    newChrIndex(1) = m;
    chrIndex.insert_rows(0, newChrIndex);
  }

  return(chrIndex);
}

void DenseGRMClass::setchrIndexLOCO(Rcpp::StringVector t_chrVecNames)
{
  arma::Mat<int> chrIndex = {0, m_M};
  Rcpp::List chrIndexLOCO = List::create(Named("none") = chrIndex);
  
  int uniqChrNum = t_chrVecNames.length();

  for(int i = 0; i < uniqChrNum; i ++){
    Rcpp::String excludeChr = t_chrVecNames[i];
    arma::Mat<int> chrIndex = makeChrIndex(excludeChr, m_chrVec);
    chrIndexLOCO.push_back(chrIndex, excludeChr);
  }
  
  m_chrIndexLOCO = chrIndexLOCO;
}

//http://gallery.rcpp.org/articles/parallel-inner-product/
struct getKinbVecParallel : public Worker
{
  // source vectors
  arma::vec bVec;
  unsigned int N;
  unsigned int M;
  DenseGRMClass* ptrDenseGRM;

  // product that I have accumulated
  arma::vec KinbVec;
  int counts;

  // constructors
  getKinbVecParallel(arma::vec bVec, DenseGRMClass* ptrDenseGRM)
    : bVec(bVec), ptrDenseGRM(ptrDenseGRM), counts(0)
  {
    M = ptrDenseGRM->getM();
    N = ptrDenseGRM->getN();
    KinbVec.zeros(N);
  }
  getKinbVecParallel(const getKinbVecParallel& getKinbVecParallel, Split)
    : bVec(getKinbVecParallel.bVec), ptrDenseGRM(getKinbVecParallel.ptrDenseGRM), counts(0)
  {
    N = getKinbVecParallel.N;
    M = getKinbVecParallel.M;
    KinbVec.zeros(getKinbVecParallel.N);
  }
  // process just the elements of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    arma::vec oneMarker(N);
    for(unsigned int i = begin; i < end; i ++){
      ptrDenseGRM->getOneMarkerStd(i, &oneMarker);
      // if we use float instead of double, the following step can be faster (x2),
      // but the package is not valid on Windows any more.
      KinbVec += oneMarker * arma::dot(oneMarker, bVec);
      counts ++;
    }
  }

  // join my value with that of another InnerProduct
  void join(const getKinbVecParallel & rhs) {
    KinbVec += rhs.KinbVec;
    counts += rhs.counts;
  }
};

arma::vec getKinbVec(arma::vec t_bVec, DenseGRMClass* t_ptrDenseGRM, string t_excludeChr, int t_grainSize)
{
  // int M = t_ptrDenseGRM->getM();
  Rcpp::List chrIndexLOCO = t_ptrDenseGRM->getChrIndexLOCO();
  
  arma::Mat<int> ChrIdx = chrIndexLOCO[t_excludeChr];
  
  // declare the InnerProduct instance that takes a pointer to the vector data
  getKinbVecParallel getKinbVecParallel(t_bVec, t_ptrDenseGRM);

  int nIdx = ChrIdx.n_rows;
  for(int i = 0; i < nIdx; i++){
    int idxStart = ChrIdx(i,0);
    int idxEnd = ChrIdx(i,1);
    // cout << idxStart << "\t" << idxEnd << endl;
    // call paralleReduce to start the work
    parallelReduce(idxStart, idxEnd, getKinbVecParallel, t_grainSize);
  }

  Rcpp::checkUserInterrupt();
  return getKinbVecParallel.KinbVec / getKinbVecParallel.counts;
}

}
