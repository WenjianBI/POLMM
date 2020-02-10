// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

// class to store genotype from plink file
class genoClass{
private:
  
  // PLINK format
  const static unsigned char HOM_REF = 0x3;  // 0b11 ;
  const static unsigned char HET = 0x2;      // 0b10 ;
  const static unsigned char HOM_ALT = 0x0;  // 0b00 ;
  const static unsigned char MISSING = 0x1;  // 0b01 ;
  
  // Basic information
  size_t M;                         // number of markers
  size_t N,N_Old;                   // number of samples, after and before filt
  size_t numArrays;                 // make arrays to store genotype data to avoid large continuous memory usage
  size_t numMarkersofEachArray;
  size_t numMarkersofLastArray;
  size_t numBytesofEachMarker;
  size_t numBytesofEachMarker_Old; 
  size_t numBytesReserve;
  
  // input file stream of .bed file
  ifstream ibedfile;  
  int bufferG1;
  unsigned char bufferG4;
  arma::vec OneMarkerG1;
  arma::vec OneMarkerG1_Old;
  std::vector<unsigned char> OneMarkerG4_Old;
  arma::vec DiagStdGeno;
  
  std::vector< std::vector<unsigned char>* > genoVecofPointers;
  
  // set environment
  void openPlink(string, arma::vec);
  void setChrVec(string, Rcpp::StringVector&, vector<string>&);
  void setArrays(float);
  void readFile(string, size_t&, string);
  float getOneMarkerPlink(vector<unsigned char>&, arma::vec);
  void setOneMarkerArray(int);
  float getinvStd(float);
  void setDiagStdGeno();
    
  arma::vec alleleFreqVec; 
  arma::vec invStdVec;
  Rcpp::StringVector chrVec;
  vector<string> chrVecParallel;
  
  // insert geno (0,1,2,3) to specific pos (0,1,2,3) of address c (1 byte)
  void setGenotype(unsigned char* c, const int pos, const int geno) {
    (*c) |= (geno << (pos << 1));
  }
  
  // extract geno (0,1,2,3) at specific pos (0,1,2,3) of address c (1 byte)  
  void getGenotype(unsigned char c, const int pos, int& geno) {
    geno = (c >> (pos << 1)) & 0b11; 
  }
  
  void setStdGenoLookUpArr(double maf, double invsd, arma::vec & stdGenoLookUpArr){
	double maf2 = 2 * maf;
	stdGenoLookUpArr(0) = (0-maf2)*invsd;
	stdGenoLookUpArr(1) = (1-maf2)*invsd;
	stdGenoLookUpArr(2) = (2-maf2)*invsd;
  }
  
public:
  void setGenoObj(string, arma::vec, float);
  int getN(){return N;}
  int getM(){return M;}
  arma::vec getalleleFreqVec(){return alleleFreqVec;}
  void getOneMarker(size_t, bool, arma::vec*);
  arma::vec* getDiagStdGeno(){return &DiagStdGeno;}
  arma::vec getAvec(){return alleleFreqVec;}
  Rcpp::StringVector getchrVec(){return chrVec;}
  vector<string> getchrVecParallel(){return chrVecParallel;}  // the same as chrVec, different data type
};


void genoClass::setGenoObj(string Plink, 
                           arma::vec posSampleInPlink, 
                           float memoryChunk)    // unit is Gb
{
  // Read in plink files
  openPlink(Plink, posSampleInPlink);
  cout << "Number of samples in Plink file:\t" << N_Old << endl;
  cout << "Number of samples:\t" << N << endl;
  cout << "Number of markers:\t" << M << endl << endl;
  
  // Allocate memory
  setArrays(memoryChunk);
  cout << "numBytesofEachMarker in Plink file:\t" << numBytesofEachMarker_Old << endl;
  cout << "numBytesofEachMarker:\t" << numBytesofEachMarker << endl;
  cout << "numBytesReserve:\t" << numBytesReserve << endl << endl;
  
  cout << "numArrays:\t" << numArrays << endl;
  cout << "numMarkersofEachArray:\t" << numMarkersofEachArray << endl;
  cout << "numMarkersofLastArray:\t" << numMarkersofLastArray << endl << endl;
  
  Rcpp::IntegerVector chrVecCounts = table(chrVec);
  Rcpp::StringVector chrVecNames = chrVecCounts.names();
  for(int i = 0; i < chrVecCounts.size(); i++){
	  cout << "Number of markers in chr " << chrVecNames(i) << ":\t" << chrVecCounts(i) << endl;
  }
  ///////////////////////////////////////// MAIN PART //////////////////////////////

  // loop for all SNPs
  float freq, invStd;
  for(size_t m=0; m<M; m++)
  {
    ibedfile.seekg(3+numBytesofEachMarker_Old*m);
    ibedfile.read((char*)(&OneMarkerG4_Old[0]), numBytesofEachMarker_Old);
    freq = getOneMarkerPlink(OneMarkerG4_Old, posSampleInPlink);
    setOneMarkerArray(m);
    invStd = getinvStd(freq);
    alleleFreqVec[m]=freq;
    invStdVec[m]=invStd;
    if((m+1) % 1000 == 0){
      cout << "Complete\t" << m+1 <<"\tSNPs!!!!" << endl;
    }
  }

  cout << endl << "Complete reading data!!" << endl;
  setDiagStdGeno();
  ibedfile.close();
}

float genoClass::getinvStd(float freq)
{
  float Std = sqrt(2*freq*(1-freq));
  if(Std == 0)
    return 0;
  else
    return 1/Std;
}

void genoClass::setOneMarkerArray(int m)
{
  // loop for all samples
  bufferG4=0;
  int pos;
  for(size_t ind=0; ind<N; ind++){
    pos = ind%4;
    bufferG1 = OneMarkerG1[ind];
    setGenotype(&bufferG4, pos, bufferG1);
    if((pos==3) | (ind==(N-1))){
      genoVecofPointers[m/numMarkersofEachArray]->push_back(bufferG4); //avoid large continuous memory usage
      bufferG4 = 0;
    }
  }
}
  
float genoClass::getOneMarkerPlink(vector<unsigned char>& OneMarkerG4_Old,
                                   arma::vec posSampleInPlink)
{
  // read in genotypes of one marker
  for(size_t ind=0; ind< N_Old; ind++){
    bufferG4 = OneMarkerG4_Old[ind/4];
    getGenotype(bufferG4, ind%4, bufferG1);
    OneMarkerG1_Old[ind] = bufferG1;
  }
  
  // select samples
  float sum=0;
  std::vector<int> indexNA;
  for(size_t index=0; index<N; index++){
    bufferG1 = OneMarkerG1_Old[posSampleInPlink[index] - 1];   // C++ start from 0
    switch(bufferG1){
      case HOM_REF: break;
      case MISSING: indexNA.push_back(index);break;
      case HET: sum+=1;break;
      case HOM_ALT: sum+=2;break;
    }
    OneMarkerG1[index] = bufferG1;
  }
  
  // genotype imputation (if needed)
  int lengthNA = indexNA.size();
  int count = N-lengthNA;
  float freq = sum/count;  // it is actually twice allele frequency
  if(lengthNA != 0){
    int bestguess = round(freq);
    unsigned char fill=0; 
    switch(bestguess){
      case 0: fill=HOM_REF;break;
      case 1: fill=HET;break;
      case 2: fill=HOM_ALT;break;
    }
    for(int i=0; i<lengthNA; i++){
      OneMarkerG1[indexNA[i]]=fill;
      // sum += fill;
      sum += bestguess;
      count += 1;
    }
    freq = sum/count;
  }
  return freq/2;
}

void genoClass::setArrays(float memoryChunk)
{
  numBytesofEachMarker = (N+3)/4;
  numBytesofEachMarker_Old = (N_Old+3)/4;
  numBytesReserve = numBytesofEachMarker*M+M*2;    // not sure about why M*2 is needed (2 bytes for each marker)
  numMarkersofEachArray = floor((memoryChunk*pow(10.0, 9.0))/numBytesofEachMarker);
  if(numMarkersofEachArray > M)
    numMarkersofEachArray=M;
  numArrays = (M-1) / numMarkersofEachArray + 1;
  numMarkersofLastArray = M - (numArrays-1)*numMarkersofEachArray;
  
  // get binray data from plink files
  OneMarkerG4_Old.reserve(numBytesofEachMarker_Old);  
  OneMarkerG4_Old.resize(numBytesofEachMarker_Old);
  OneMarkerG1_Old.zeros(N_Old);             
  OneMarkerG1.zeros(N);
  
  // set genoVecofPointers
  genoVecofPointers.resize(numArrays);
  for (size_t i = 0; i < numArrays-1 ; i++){
    genoVecofPointers[i] = new vector<unsigned char>;
    genoVecofPointers[i]->reserve(numMarkersofEachArray*numBytesofEachMarker);
  }
  genoVecofPointers[numArrays-1] = new vector<unsigned char>;
  genoVecofPointers[numArrays-1]->reserve(numMarkersofLastArray*numBytesofEachMarker);
  
  alleleFreqVec.zeros(M);
  invStdVec.zeros(M);
  
  // catch(std::bad_alloc& ba)
  // {
  //   std::cerr << "bad_alloc caught1: " << ba.what() << '\n';
  //   exit(EXIT_FAILURE);
  // }
}

void genoClass::openPlink(string Plink, arma::vec posSampleInPlink)
{
  string bimfile = Plink + ".bim";
  string famfile = Plink + ".fam";
  string bedfile = Plink + ".bed";
  readFile(famfile, N_Old, "Error! fam file not open!");
  if((max(posSampleInPlink)>N_Old) | (min(posSampleInPlink)<1)) 
    stop("posSampleInPlink should be between 1 and N");
  readFile(bimfile, M, "Error! bim file not open!");
  ibedfile.open(bedfile.c_str(), ios::binary);
  if (!ibedfile.is_open())
    stop("Error! bed file not open!");
  N = posSampleInPlink.size();
  Rcpp::StringVector chrVecTemp(M);
  chrVecParallel.resize(M);
  
  setChrVec(bimfile, chrVecTemp, chrVecParallel);
  chrVec = chrVecTemp;
}

void genoClass::setChrVec(string bimfile, Rcpp::StringVector& chrVecTemp, vector<string>& chrVecParallel){
  string temp;
  ifstream ifile;

  ifile.open(bimfile.c_str());
  int count = 0;
  while(getline(ifile, temp)){
	string chr = temp.substr(0, temp.find('\t'));
    chrVecTemp[count] = Rcpp::String(chr);
	chrVecParallel[count] = chr;
	count++;
 }
  ifile.close();
}

void genoClass::readFile(string file, size_t& n, string errInfo){
  string junk;
  ifstream ifile;
  size_t count=0;
  
  ifile.open(file.c_str());
  if(!ifile.is_open())
    stop(errInfo);
  while(getline(ifile,junk))
    count++;
  n = count;
  ifile.close();
}

void genoClass::getOneMarker(size_t m, bool ifStd, arma::vec* oneMarker){
  // avoid large continuous memory usage
  int indexOvectorPointer = m / numMarkersofEachArray;
  int markerIdxinVec = m % numMarkersofEachArray;
  // set up start byte index and end byte index
  size_t startBytesIdxinVec = numBytesofEachMarker * markerIdxinVec;
  size_t endBytesIdxinVec = startBytesIdxinVec+numBytesofEachMarker;
  size_t ind= 0;
  // use bufferG4b to replace bufferG4 in case of parallele computation
  unsigned char bufferG4b;
  size_t bufferG1b;
  size_t a,b,i,j;
  
  if(ifStd){
	 arma::vec stdGenoLookUpArr(3);
	 setStdGenoLookUpArr(alleleFreqVec[m], invStdVec[m], stdGenoLookUpArr);
	 for(i=startBytesIdxinVec; i< endBytesIdxinVec; i++){
		 bufferG4b = genoVecofPointers[indexOvectorPointer]->at(i); //avoid large continuous memory usage
		 for(j=0; (j<4)&(ind<N); j++,ind++){
			 b = bufferG4b & 1 ;
			 bufferG4b = bufferG4b >> 1;
			 a = bufferG4b & 1 ;
			 bufferG4b = bufferG4b >> 1;
			 bufferG1b = 2-(a+b);
			 oneMarker->at(ind) = stdGenoLookUpArr(bufferG1b);
		}
	} 
  }else{
	 for(i=startBytesIdxinVec; i< endBytesIdxinVec; i++){
		 bufferG4b = genoVecofPointers[indexOvectorPointer]->at(i); //avoid large continuous memory usage
		 for(j=0; (j<4)&(ind<N); j++,ind++){
			 b = bufferG4b & 1 ;
			 bufferG4b = bufferG4b >> 1;
			 a = bufferG4b & 1 ;
			 bufferG4b = bufferG4b >> 1;
			 bufferG1b = 2-(a+b);
			 oneMarker->at(ind) = bufferG1b;
		}
	} 
  }
}

void genoClass::setDiagStdGeno(){
  DiagStdGeno.zeros(N);
  arma::vec oneMarker(N);
  for(size_t m=0; m<M; m++){
    getOneMarker(m, 1, &oneMarker);  // write Standard Genotype to OneMarkerG1
    DiagStdGeno = DiagStdGeno + (oneMarker) % (oneMarker);
  }
}

float innerProduct(arma::vec& x, arma::vec& y){
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

//http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

//http://gallery.rcpp.org/articles/parallel-inner-product/
struct getKinbVecParallel : public Worker
{
  // source vectors
  arma::vec& bVec;
  unsigned int N;
  unsigned int M;
  genoClass* ptrGeno;

  // product that I have accumulated
  arma::vec KinbVec;
  int counts;

  // constructors
  getKinbVecParallel(arma::vec& bVec, genoClass* ptrGeno)
    : bVec(bVec), ptrGeno(ptrGeno) 
  {
    M = ptrGeno->getM();
    N = ptrGeno->getN();
	KinbVec.zeros(N);
	counts = 0;
  }
  getKinbVecParallel(const getKinbVecParallel& getKinbVecParallel, Split)
    : bVec(getKinbVecParallel.bVec), ptrGeno(getKinbVecParallel.ptrGeno)
  {
    N = getKinbVecParallel.N;
    M = getKinbVecParallel.M;
    KinbVec.zeros(getKinbVecParallel.N);
	counts = 0;
  }
  // process just the elements of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    arma::vec oneMarker(N);
    for(unsigned int m = begin; m < end; m++){
		ptrGeno->getOneMarker(m, 1, &oneMarker);
		// KinbVec += oneMarker * innerProduct(oneMarker, bVec);
		KinbVec += oneMarker * arma::dot(oneMarker, bVec);
		counts++;
    }
  }

  // join my value with that of another InnerProduct
  void join(const getKinbVecParallel & rhs) {
    KinbVec += rhs.KinbVec;
	counts += rhs.counts;
  }
};

arma::Mat<int> makeChrIdx(string excludeChr, vector<string> chrVecParallel, size_t M){
	if(excludeChr=="n"){
		arma::Mat<int> ChrIdx(1,2);
		ChrIdx(0,0) = 0;
		ChrIdx(0,1) = M;
		return(ChrIdx);
	}
	
	arma::Mat<int> ChrIdx;
	int idxStart = 0;
	for(int i = 0; i < M; i++){
		if(chrVecParallel[i] == excludeChr){
			if(idxStart != i){
				arma::Row<int> newChrIdx(2);
				newChrIdx(0) = idxStart;
				newChrIdx(1) = i;
				ChrIdx.insert_rows(0, newChrIdx);
			}
			idxStart = i + 1;
		}
	}
	
	if(idxStart != M){
		arma::Row<int> newChrIdx(2);
		newChrIdx(0) = idxStart;
		newChrIdx(1) = M;
		ChrIdx.insert_rows(0, newChrIdx);
	}
	
	return(ChrIdx);
}

arma::vec getKinbVec(arma::vec bVec, genoClass* ptrGeno, string excludeChr) {
  
  size_t M = ptrGeno->getM();
  size_t N = ptrGeno->getN();
  vector<string> chrVecParallel = ptrGeno->getchrVecParallel();
  
  // declare the InnerProduct instance that takes a pointer to the vector data
  getKinbVecParallel getKinbVecParallel(bVec, ptrGeno);
  
  // arma::vec KinbVec(N, arma::fill::zeros);
  // int counts = 0;
  
  arma::Mat<int> ChrIdx = makeChrIdx(excludeChr, chrVecParallel, M);
  // cout << "ChrIdx:" << endl << ChrIdx << endl;
  
  int nIdx = ChrIdx.n_rows;
  for(int i = 0; i < nIdx; i++){
	  int idxStart = ChrIdx(i,0);
	  int idxEnd = ChrIdx(i,1);
	  // call paralleReduce to start the work
      parallelReduce(idxStart, idxEnd, getKinbVecParallel);
	  // KinbVec += getKinbVecParallel.KinbVec;
	  // counts += getKinbVecParallel.counts;
  }
  
  // cout << "getKinbVecParallel.counts:\t" << getKinbVecParallel.counts << endl;
  // return the computed product
  // return KinbVec/counts;
  return getKinbVecParallel.KinbVec/getKinbVecParallel.counts;
}
