
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include <RcppArmadillo.h>

#include <vector>
#include <cmath>
#include <RcppArmadilloExtensions/sample.h> // sample

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


arma::vec nb(int n){
  return(Rcpp::rbinom(n,1,0.5));
}

double get_cpu_time(){
  return (double)clock() / CLOCKS_PER_SEC;
}

double calCV(arma::vec xVec){
  int n = xVec.size();
  double Mean = arma::mean(xVec);
  double Sd = arma::stddev(xVec);
  double CV = (Sd/Mean)/n;
  return(CV);
}

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat)
{
  double innerProd = arma::accu(x1Mat % x2Mat);
  return(innerProd);
}

std::vector<int> whichCPP(Rcpp::StringVector chrVec, std::string excludechr)
{
  std::vector<int> indexSNPs;
  for(int i = 0; i < chrVec.size(); i++){
    if(std::string(chrVec(i))==excludechr){
      indexSNPs.push_back(i);
    }
  }
  return(indexSNPs);
}

arma::mat getadjGFast(arma::vec GVec,
                      arma::mat XXR_Psi_RX,     
                      arma::mat XR_Psi_R_new, // XR_Psi_R_new ( p x n ), sum up XR_Psi_R ( p x n(J-1) ) for each subject 
                      int n, int J, int p)
{
  arma::mat adjGMat(n, J-1);
  // To increase computational efficiency when lots of GVec elements are 0
  arma::vec XR_Psi_RG1(p, arma::fill::zeros);
  for(int i = 0; i < n; i++){
    if(GVec(i)!=0){
      XR_Psi_RG1 += XR_Psi_R_new.col(i) * GVec(i);
    }
  }
  
  arma::vec TempVec = XXR_Psi_RX * XR_Psi_RG1;   // n(J-1) x 1
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      adjGMat(i,j) = GVec(i) - TempVec(index);
      index++;
    }
  }
  return(adjGMat);
}

double getStatFast(arma::vec GVec,         // n x 1
                   arma::mat RymuMat,      // n x (J-1): R %*% (y - mu)
                   int n, int J)
{
  double Stat = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      Stat += GVec(i) * RymuMat(i,j);
    }
  }
  return(Stat);
}

double getVarWFast(arma::mat adjGMat,  // n x (J-1)
                   arma::cube RPsiR,   // (J-1) x (J-1) x n
                   int n, int J)
{
  double VarW = 0;
  for(int i = 0; i < n; i++){
    for(int j1 = 0; j1 < J-1; j1++){
      VarW += RPsiR(j1,j1,i) * adjGMat(i,j1) * adjGMat(i,j1);
      for(int j2 = j1+1; j2 < J-1; j2++){
        VarW += 2 * RPsiR(j1,j2,i) * adjGMat(i,j1) * adjGMat(i,j2);
      }
    }
  }
  return(VarW);
}

// [[Rcpp::export]]
Rcpp::List outputadjGFast(arma::vec GVec,
                          Rcpp::List objP)
{
  arma::mat adjGMat = getadjGFast(GVec, objP["XXR_Psi_RX"], objP["XR_Psi_R_new"], objP["n"], objP["J"], objP["p"]);
  double Stat = getStatFast(GVec, objP["RymuMat"], objP["n"], objP["J"]);
  double VarW = getVarWFast(adjGMat, objP["RPsiR"], objP["n"], objP["J"]);
  
  Rcpp::List outList = List::create(Named("adjGMat")=adjGMat,
                                    Named("Stat")=Stat,           
                                    Named("VarW")=VarW);
  
  return(outList);
}

// [[Rcpp::export]]
arma::cube getRPsiR_v1(arma::mat muMat,
                       arma::mat iRMat,
                       int n, int J)
{
  arma::cube RPsiR(J-1, J-1, n, arma::fill::zeros);
  arma::mat muRMat = muMat.cols(0, J-2) / iRMat;
  for(int i = 0; i < n; i++){
    for(int j1 = 0; j1 < J-1; j1++){
      RPsiR(j1,j1,i) += muRMat(i,j1) / iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      for(int j2 = j1+1; j2 < J-1; j2++){
        RPsiR(j1,j2,i) -= muRMat(i,j1) * muRMat(i,j2);
        // RPsiR(j2,j1,i) -= muRMat(i,j1) * muRMat(i,j2); // v4.1: RPsiR is only used in function getVarWFast, 
      }
    }
  }
  return(RPsiR);
}

// [[Rcpp::export]]
Rcpp::List outputadjGFast_v1(arma::vec GVec,
                             Rcpp::List objP,
                             arma::cube RPsiR)
{
  arma::mat adjGMat = getadjGFast(GVec, objP["XXR_Psi_RX"], objP["XR_Psi_R_new"], objP["n"], objP["J"], objP["p"]);
  double Stat = getStatFast(GVec, objP["RymuMat"], objP["n"], objP["J"]);
  double VarW = getVarWFast(adjGMat, RPsiR, objP["n"], objP["J"]);
  
  Rcpp::List outList = List::create(Named("adjGMat")=adjGMat,
                                    Named("Stat")=Stat,           
                                    Named("VarW")=VarW);
  
  return(outList);
}

class nullModelClass
{
private:
  
  // used as dimension
  int n, J, p, M;
  genoClass* ptrGeno;
  Rcpp::List initParList, controlList;
  
  // constant vector and matrix over the process
  arma::Col<int> yVec; // should start from 1, not 0
  arma::mat yMat, Cova, CovaMat, GMatRatio;
  
  // parameters
  arma::vec beta, bVec, eps;
  arma::vec beta0, eps0;
  double tau, tau0;
  double diffBeta, difftau, diffeps;
  
  // arguments
  int maxiter, maxiterPCG, maxiterEps, tracenrun, seed, nSNPsVarRatio;
  double tolBeta, tolTau, tolPCG, tolEps, minMafVarRatio, CVcutoff;
  bool LOCO;
  Rcpp::List LOCOList;
  
  // working vectors/matrix
  arma::mat WMat, muMat, mMat, nuMat, iRMat, YMat, iSigma_CovaMat, iSigmaX_XSigmaX;
  arma::vec eta, iSigma_YVec, iSigma_VPYVec;
  
  arma::mat TraceRandMat, V_TRM, iSigma_V_TRM;
  
  // functions in setNullModel()
  void setParList(Rcpp::List);    // set up parameter list
  void setcontrolList(Rcpp::List);
  void setArray();               // set up dimension 
  arma::mat getCovaMat();
  arma::mat getyMat();
  void getTraceRandMat();
  arma::vec ZMat(arma::vec), tZMat(arma::vec);
  
  // functions in fitNullModel()
  void updateMats();
  arma::mat getiPsixMat(arma::mat);
  arma::vec getRowSums(arma::mat);
  arma::vec convert1(arma::mat xMat);
  arma::mat convert2(arma::vec xMat);
  void getPCGofSigmaAndCovaMat(arma::mat, arma::mat&, string);
  void getPCGofSigmaAndVector(arma::vec, arma::vec&, string);
  void updateParaConv(string);
  void updatePara(string);
  void updateEpsOneStep();
  void updateEps();
  arma::mat getSigmaxMat(arma::mat, string);
  arma::cube getInvBlockDiagSigma();
  arma::mat solverBlockDiagSigma(arma::cube&, arma::mat& xMat);
  void updateTau();
  
  // functions if LOCO = TRUE
  arma::mat getPsixMat(arma::mat);
  Rcpp::List getLOCO(string, std::vector<int>);
  arma::mat getVarRatio(std::vector<int>, string, Rcpp::List);
  arma::rowvec getVarOneSNP(arma::vec, string, Rcpp::List);
  arma::rowvec getVarOneSNP_v1(arma::vec, string, Rcpp::List, arma::cube);
  double getVarP(arma::mat, string);
  Rcpp::List getobjP();
  void getRPsiR(arma::cube&);
  
  // functions if LOCO = FALSE
  arma::mat getVarRatio(arma::mat, Rcpp::List);
  
public:
  
  void setNullModel(arma::mat, arma::Col<int>, genoClass*, arma::vec, arma::vec, arma::vec, double, arma::mat, Rcpp::List);
  void fitNullModel();
  Rcpp::List getNullModel();
  
};

void nullModelClass::fitNullModel()
{
  // initial vector
  updateMats();
  
  // start iteration
  cout << "Start iteration ....." << endl;
  int iter = 0;
  while(iter < maxiter)
  {
    cout << "iter " << iter << endl;
    
    // update other parameters until converge
    updateParaConv("n");
    
    // update tau
    tau0 = tau;
    updateTau();
    
    if(std::isnan(tau))
      stop("Parameter tau is NA.");
    
    cout << "beta: " << endl << beta << endl;
    cout << "tau: " << tau << endl << endl;
    iter++;
    
    double diffTau = abs(tau - tau0)/(abs(tau) + abs(tau0) + tolTau);
    if(diffTau < tolTau)
      break;
  }
  
  if(LOCO){
    Rcpp::StringVector chrVec = ptrGeno->getchrVec();
    Rcpp::StringVector uniqchr = unique(chrVec);
    cout << uniqchr << endl;
    for(int i = 0; i < uniqchr.size(); i++){
      string excludechr = string(uniqchr(i));
      std::vector<int> indexSNPs = whichCPP(chrVec, excludechr);
      indexSNPs = Rcpp::RcppArmadillo::sample(indexSNPs, indexSNPs.size(), FALSE);
      cout << endl << "Leave One Chromosome Out: Chr " << excludechr << endl;
      cout << "Number of SNPs in this chromosome:\t" << indexSNPs.size() << endl;
      Rcpp::List tempLOCO = getLOCO(excludechr, indexSNPs);
      LOCOList[excludechr] = tempLOCO;
    }
  }else{
    Rcpp::List objP = getobjP();
    // output variance matrix
    arma::mat VarRatioMat = getVarRatio(GMatRatio, objP);
    double VarRatio = arma::mean(VarRatioMat.col(4));
    
    Rcpp::List temp = List::create(Named("objP")=objP,
                                   Named("VarRatioMat")=VarRatioMat,
                                   Named("VarRatio")=VarRatio);
    LOCOList["LOCO=F"] = temp;
  }
}

// outMat = PsiMat %*% xMat, PsiMat is determined by muMat
arma::mat nullModelClass::getPsixMat(arma::mat xMat)    // matrix with dim of n x (J-1)
{
  arma::mat Psi_xMat(n, J-1);
  // loop for samples
  for(int i = 0; i < n; i++){
    arma::rowvec muVec(J-1);
    for(int j = 0; j < J-1; j++){
      Psi_xMat(i,j) = muMat(i,j) * xMat(i,j);
      muVec(j) = muMat(i,j);
    }
    double sum_mu_x = sum(Psi_xMat.row(i));
    Psi_xMat.row(i) -= muVec * sum_mu_x; 
  }
  return(Psi_xMat);
}

void nullModelClass::getRPsiR(arma::cube& RPsiR)   // (J-1) x (J-1) x n
{
  arma::mat muRMat = muMat.cols(0, J-2) / iRMat;
  for(int i = 0; i < n; i++){
    for(int j1 = 0; j1 < J-1; j1++){
      RPsiR(j1,j1,i) += muRMat(i,j1) / iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      for(int j2 = j1+1; j2 < J-1; j2++){
        RPsiR(j1,j2,i) -= muRMat(i,j1) * muRMat(i,j2);
        // RPsiR(j2,j1,i) -= muRMat(i,j1) * muRMat(i,j2); // v4.1: RPsiR is only used in function getVarWFast, 
      }
    }
  }
}

// get a list for p value calculation in step 2
Rcpp::List nullModelClass::getobjP()
{
  // output for Step 2
  arma::mat XR_Psi_R(p, n*(J-1));       // p x n(J-1)
  arma::mat xMat(n, J-1);
  arma::vec temp(n*(J-1));
  for(int k = 0; k < p; k++){
    xMat = convert2(CovaMat.col(k));
    temp = convert1(getPsixMat(xMat / iRMat) / iRMat);
    XR_Psi_R.row(k) = temp.t();
  }
  arma::mat XXR_Psi_RX = CovaMat * inv(XR_Psi_R * CovaMat);             // n(J-1) x p
  iSigmaX_XSigmaX = iSigma_CovaMat * inv(CovaMat.t() * iSigma_CovaMat);
  
  // arma::cube RPsiR(J-1, J-1, n, arma::fill::zeros);
  // getRPsiR(RPsiR);
  
  arma::mat XR_Psi_R_new(p, n, arma::fill::zeros);  // p x n
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      XR_Psi_R_new.col(i) += XR_Psi_R.col(index);
      index++;
    }
  }
  
  arma::mat ymuMat = yMat - muMat;
  arma::mat RymuMat = ymuMat.cols(0, J-2) / iRMat; // n x (J-1): R %*% (y - mu)
  
  Rcpp::List objP = List::create(Named("n")=n,
                                 Named("J")=J,
                                 Named("p")=p,
                                 Named("XXR_Psi_RX")=XXR_Psi_RX,
                                 // Named("RPsiR")=RPsiR,
                                 Named("XR_Psi_R_new")=XR_Psi_R_new,           
                                 Named("RymuMat")=RymuMat,
                                 Named("muMat")=muMat,
                                 Named("iRMat")=iRMat);
  
  return(objP);
}

double nullModelClass::getVarP(arma::mat adjGMat,
                               string excludechr)
{
  arma::vec adjGVec = convert1(adjGMat);
  arma::vec iSigmaGVec(n*(J-1), arma::fill::zeros);
  getPCGofSigmaAndVector(adjGVec, iSigmaGVec, excludechr);
  double VarP = as_scalar(adjGVec.t() * (iSigmaGVec - iSigmaX_XSigmaX * (CovaMat.t() * iSigmaGVec)));
  return(VarP);
}

arma::rowvec nullModelClass::getVarOneSNP(arma::vec GVec,
                                          string excludechr,
                                          Rcpp::List objP)
{
  arma::rowvec VarOut(5);
  double AF = sum(GVec) / GVec.size() / 2;
  if(AF > 0.5)
    AF = 1 - AF;
  
  Rcpp::List adjGList = outputadjGFast(GVec, objP);
  arma::mat adjGMat = adjGList["adjGMat"];
  double Stat = adjGList["Stat"];
  double VarW = adjGList["VarW"];
  double VarP = getVarP(adjGMat, excludechr);
  
  VarOut(0) = AF;
  VarOut(1) = Stat;
  VarOut(2) = VarW;
  VarOut(3) = VarP;
  VarOut(4) = VarP/VarW;
  return(VarOut);
}

arma::rowvec nullModelClass::getVarOneSNP_v1(arma::vec GVec,
                                             string excludechr,
                                             Rcpp::List objP,
                                             arma::cube RPsiR)
{
  arma::rowvec VarOut(5);
  double AF = sum(GVec) / GVec.size() / 2;
  if(AF > 0.5)
    AF = 1 - AF;
  
  Rcpp::List adjGList = outputadjGFast_v1(GVec, objP, RPsiR);
  arma::mat adjGMat = adjGList["adjGMat"];
  double Stat = adjGList["Stat"];
  double VarW = adjGList["VarW"];
  double VarP = getVarP(adjGMat, excludechr);
  
  VarOut(0) = AF;
  VarOut(1) = Stat;
  VarOut(2) = VarW;
  VarOut(3) = VarP;
  VarOut(4) = VarP/VarW;
  return(VarOut);
}

arma::mat nullModelClass::getVarRatio(arma::mat GMatRatio,
                                      Rcpp::List objP)
{
  arma::vec GVec(n);
  arma::rowvec VarOneSNP(5);
  
  arma::mat VarRatioMat(nSNPsVarRatio, 5);
  arma::mat newVarRatio(10, 5);
  
  // Rcpp::List adjGList = outputadjGFast(GVec, objP);
  
  int index = 0;
  int indexTot = 0;
  while(index < nSNPsVarRatio){
    indexTot++;
    GVec = GMatRatio.col(index);
    VarOneSNP = getVarOneSNP(GVec, "n", objP);
    VarRatioMat.row(index) = VarOneSNP;
    index++;
  }
  
  arma::vec VarRatio = VarRatioMat.col(4);
  double CV = calCV(VarRatio);
  cout << "nSNPs for CV: " << index << endl;
  cout << "CV: " << CV << endl;
  
  while(CV > CVcutoff && VarRatioMat.n_rows < 100){
    int indexTemp = 0;
    while(indexTemp < 10){
      indexTot++;
      GVec = GMatRatio.col(indexTot);
      VarOneSNP = getVarOneSNP(GVec, "n", objP);
      newVarRatio.row(indexTemp) = VarOneSNP;
      index++;
      indexTemp++;
    }
    VarRatioMat.insert_rows(0, newVarRatio);
    arma::vec VarRatio = VarRatioMat.col(4);
    CV = calCV(VarRatio);
    cout << "nSNPs for CV: " << index << endl;
    cout << "CV: " << CV << endl;
  }
  return(VarRatioMat);
}

arma::mat nullModelClass::getVarRatio(std::vector<int> indexSNPs,
                                      string excludechr,
                                      Rcpp::List objP)
{
  arma::vec alleleFreqVec = ptrGeno->getalleleFreqVec();
  
  arma::vec GVec(n);
  arma::rowvec VarOneSNP(5);
  
  arma::mat VarRatioMat(nSNPsVarRatio, 5);
  arma::mat newVarRatio(10, 5);
  
  // Rcpp::List adjGList = outputadjGFast(GVec, objP);
  arma::cube RPsiR = getRPsiR_v1(objP["muMat"], objP["iRMat"], objP["n"], objP["J"]);
  
  int index = 0;
  int indexTot = 0;
  while(index < nSNPsVarRatio){
    int m = indexSNPs[indexTot];
    indexTot++;
    if(alleleFreqVec(m) > minMafVarRatio && alleleFreqVec(m) < 1-minMafVarRatio){
      ptrGeno->getOneMarker(m, 0, &GVec);
      // VarOneSNP = getVarOneSNP(GVec, excludechr, objP);
      VarOneSNP = getVarOneSNP_v1(GVec, excludechr, objP, RPsiR);
      VarRatioMat.row(index) = VarOneSNP;
      index++;
    }
  }
  
  arma::vec VarRatio = VarRatioMat.col(4);
  double CV = calCV(VarRatio);
  cout << "nSNPs for CV: " << index << endl;
  cout << "CV: " << CV << endl;
  
  while(CV > CVcutoff && VarRatioMat.n_rows < 100){
    int indexTemp = 0;
    while(indexTemp < 10){
      int m = indexSNPs[indexTot]; // round double to int
      indexTot++;
      if(alleleFreqVec(m) > minMafVarRatio && alleleFreqVec(m) < 1-minMafVarRatio){
        ptrGeno->getOneMarker(m, 0, &GVec);
        // VarOneSNP = getVarOneSNP(GVec, excludechr, objP);
        VarOneSNP = getVarOneSNP_v1(GVec, excludechr, objP, RPsiR);
        newVarRatio.row(indexTemp) = VarOneSNP;
        index++;
        indexTemp++;
      }
    }
    VarRatioMat.insert_rows(0, newVarRatio);
    arma::vec VarRatio = VarRatioMat.col(4);
    CV = calCV(VarRatio);
    cout << "nSNPs for CV: " << index << endl;
    cout << "CV: " << CV << endl;
  }
  return(VarRatioMat);
}

Rcpp::List nullModelClass::getLOCO(string excludechr,
                                   std::vector<int> indexSNPs)
{
  // estimate LOCO parameter
  updateParaConv(excludechr);
  
  // get a list for p value calculation in step 2
  Rcpp::List objP = getobjP();
  
  // output variance matrix
  arma::mat VarRatioMat = getVarRatio(indexSNPs, excludechr, objP);
  double VarRatio = arma::mean(VarRatioMat.col(4));
  
  Rcpp::List outLOCO = List::create(Named("objP")=objP,
                                    Named("VarRatioMat")=VarRatioMat,
                                    Named("VarRatio")=VarRatio);
  return(outLOCO);
}

void nullModelClass::setNullModel(arma::mat CovaR,
                                  arma::Col<int> yVecR,     // should be from 1 to J
                                  genoClass* ptrGenoInput,
                                  arma::vec betaR,
                                  arma::vec bVecR,
                                  arma::vec epsR,
                                  double tauR,
                                  arma::mat GMatRatioR,
                                  Rcpp::List controlListR)
{
  n = CovaR.n_rows;
  p = CovaR.n_cols;
  J = max(yVecR);
  M = ptrGenoInput->getM();
  
  Cova = CovaR;
  CovaMat = getCovaMat();
  yVec = yVecR;
  yMat = getyMat();
  ptrGeno = ptrGenoInput;
  
  beta = betaR;
  bVec = bVecR;
  eps = epsR;
  tau = tauR;
  GMatRatio = GMatRatioR;
  
  initParList = List::create(Named("betaR")=betaR,
                             Named("bVecR")=bVecR,
                             Named("epsR")=epsR,
                             Named("tauR")=tauR);
  
  setcontrolList(controlListR);
  setArray();
  set_seed(seed);
  getTraceRandMat();
}

arma::mat nullModelClass::getCovaMat()
{
  arma::mat CovaMat(n*(J-1), p);
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      CovaMat.row(index) = Cova.row(i);
      index++;
    }
  }
  return(CovaMat);
}

// yMat: matrix with dim of n x J
arma::mat nullModelClass::getyMat()
{
  arma::mat yMat(n, J, arma::fill::zeros);
  for(int i = 0; i < n; i++){
    yMat(i, yVec(i)-1) = 1;
  }
  return(yMat);
}

void nullModelClass::setcontrolList(Rcpp::List controlListR)
{
  controlList = controlListR;
  maxiter = controlList["maxiter"];
  maxiterPCG = controlList["maxiterPCG"]; 
  maxiterEps = controlList["maxiterEps"];
  tolBeta = controlList["tolBeta"]; 
  tolTau = controlList["tolTau"]; 
  tolPCG = controlList["tolPCG"]; 
  tolEps = controlList["tolEps"];
  tracenrun = controlList["tracenrun"]; 
  seed = controlList["seed"];
  minMafVarRatio = controlList["minMafVarRatio"];
  nSNPsVarRatio = controlList["nSNPsVarRatio"];
  CVcutoff = controlList["CVcutoff"];
  LOCO = controlList["LOCO"];
}

void nullModelClass::setArray()
{
  WMat.zeros(n,J);    
  muMat.zeros(n,J);   
  mMat.zeros(n,J);
  nuMat.zeros(n,J);
  iRMat.zeros(n,J-1);
  YMat.zeros(n,J-1);
  iSigma_CovaMat.zeros(n*(J-1), p);
  //
  eta.zeros(n);
  iSigma_YVec.zeros(n*(J-1));
  iSigma_VPYVec.zeros(n*(J-1));
  //
  TraceRandMat.zeros(n*(J-1), tracenrun);
  V_TRM.zeros(n*(J-1), tracenrun);
  iSigma_V_TRM.zeros(n*(J-1), tracenrun);
}

void nullModelClass::getTraceRandMat()
{
  for(int itrace = 0; itrace < tracenrun; itrace++){
    cout << "itrace:\t" << itrace << endl;
    arma::vec uVec = nb(n*(J-1));
    uVec = uVec*2 - 1;
    TraceRandMat.col(itrace) = uVec;
    
    double cpuin  = get_cpu_time();
    V_TRM.col(itrace) = tZMat(getKinbVec(ZMat(uVec), ptrGeno, "n"));
    double cpuout  = get_cpu_time();
    cout << "CPU Time  in getKinbVec = " << cpuout - cpuin  << endl;
  }
}

// convert from n(J-1) to n by summing up each (J-1) elements
arma::vec nullModelClass::ZMat(arma::vec xVec)
{
  arma::vec y1Vec(n, arma::fill::zeros);
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      y1Vec(i) += xVec(index);
      index++;
    }
  }
  return(y1Vec);
}

// convert from n to n(J-1) by duplicating each element for (J-1) times
arma::vec nullModelClass::tZMat(arma::vec xVec)
{
  arma::vec y1Vec(n*(J-1));
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      y1Vec(index) = xVec(i);
      index++;
    }
  }
  return(y1Vec);
}


void nullModelClass::updateTau()
{
  arma::vec YVec = convert1(YMat);
  getPCGofSigmaAndCovaMat(CovaMat, iSigma_CovaMat, "n");
  getPCGofSigmaAndVector(YVec, iSigma_YVec, "n"); 
  arma::mat Cova_iSigma_CovaMat = CovaMat.t() * iSigma_CovaMat;
  arma::vec Cova_iSigma_YVec = CovaMat.t() * iSigma_YVec;
  arma::mat iSigmaXXiSigmaX = iSigma_CovaMat * inv(Cova_iSigma_CovaMat);
  arma::vec PYVec = iSigma_YVec - iSigmaXXiSigmaX * Cova_iSigma_YVec;
  arma::vec VPYVec = tZMat(getKinbVec(ZMat(PYVec), ptrGeno, "n"));
  getPCGofSigmaAndVector(VPYVec, iSigma_VPYVec, "n"); 
  arma::vec Cova_iSigma_VPYVec = CovaMat.t() * iSigma_VPYVec;
  arma::vec PVPYVec = iSigma_VPYVec - iSigmaXXiSigmaX * Cova_iSigma_VPYVec;
  double YPVPY = as_scalar(YVec.t() * PVPYVec);
  double YPVPVPY = as_scalar(VPYVec.t() * PVPYVec);
  // The below is to calculate trace
  getPCGofSigmaAndCovaMat(V_TRM, iSigma_V_TRM, "n");
  double tracePV = 0;
  int m = TraceRandMat.n_cols;
  for(int i = 0; i < m; i++){
    arma::vec iSigma_V_TRM_col = iSigma_V_TRM.col(i);
    arma::vec P_V_TRM_col = iSigma_V_TRM_col - iSigmaXXiSigmaX * (CovaMat.t() * iSigma_V_TRM_col);
    tracePV += as_scalar(TraceRandMat.col(i).t() * P_V_TRM_col);
  }
  tracePV /= m;
  // final step
  double deriv = 0.5 * YPVPY - 0.5 * tracePV;
  double AI = 0.5 * YPVPVPY;
  double dtau = deriv / AI;
  tau0 = tau;
  tau = tau0 + dtau;
  while(tau < 0){
    dtau = dtau / 2;
    tau = tau0 + dtau;
  }
  if(tau < 1e-4){
    tau = 0;
  }
}

// update (eta, WMat, muMat, mMat, nuMat, iRMat, YMat) based on (beta, bVec, eps)
void nullModelClass::updateMats()
{
  // update (eta)
  eta = Cova * beta + bVec;
  // update (WMat, muMat, mMat, nuMat)
  double tmpExp, tmpnu0, tmpnu1;
  for(int i = 0; i < n; i++){  // loop for samples
    tmpnu0 = 0;  // eps_0 = -Inf
    for(int j = 0; j < J-1; j++){ // loop from eps_1 to eps_{J-1}
      tmpExp = exp(eps(j) - eta(i));
      tmpnu1 = tmpExp / (1+tmpExp);
      muMat(i,j) = tmpnu1 - tmpnu0;
      WMat(i,j) = tmpnu1 * (1 - tmpnu1);
      mMat(i,j) = tmpnu1 + tmpnu0 - 1;
      nuMat(i,j) = tmpnu1;
      tmpnu0 = tmpnu1;
    }
    int j = J-1;      // eps_J = Inf
    tmpnu1 = 1;
    muMat(i,j) = tmpnu1 - tmpnu0;
    WMat(i,j) = tmpnu1 * (1 - tmpnu1);
    mMat(i,j) = tmpnu1 + tmpnu0 - 1;
    nuMat(i,j) = tmpnu1;
  }
  // update (iRMat)
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      iRMat(i,j) = 1 / (mMat(i,j) - mMat(i,J-1));
    }
  }
  // update (YMat)
  arma::mat xMat = yMat.cols(0,J-2) - muMat.cols(0,J-2);
  arma::mat iPsi_xMat = getiPsixMat(xMat);
  for(int i = 0; i < n; i++){  // loop for samples
    for(int j = 0; j < J-1; j++)
      YMat(i,j) = eta(i) + (iRMat(i,j) * iPsi_xMat(i,j));
  }
}

// outMat = iPsiMat %*% xMat, iPsiMat is determined by muMat
arma::mat nullModelClass::getiPsixMat(arma::mat xMat)   // matrix with dim of n x (J-1)
{
  arma::mat iPsi_xMat(n, J-1);
  for(int i = 0; i < n; i++){   // loop for samples
    double sumx = sum(xMat.row(i));
    for(int j = 0; j < J-1; j++){
      iPsi_xMat(i,j) = sumx / muMat(i,J-1) + xMat(i,j) / muMat(i,j);
    }
  }
  return(iPsi_xMat);
}

// update parameters (except tau) until converge
void nullModelClass::updateParaConv(string excludechr)
{
  int iter = 0;
  while(iter < maxiter)
  {
    beta0 = beta;
    // update beta and bVec
    updatePara(excludechr);
    updateMats();
    
    // update eps (cutpoints)
    updateEps();
    updateMats();
    
    cout << "beta: " << endl << beta << endl;
    iter++;
    
    diffBeta = max(abs(beta - beta0)/(abs(beta) + abs(beta0) + tolBeta));
    if(diffBeta < tolBeta)
      break;
  }
}

void nullModelClass::updateEps()
{
  int iter = 0;
  while(iter < maxiterEps)
  {
    arma::vec eps0 = eps;
    updateEpsOneStep();
    updateMats();
    
    diffeps = max(abs(eps - eps0));
    if(diffeps < tolEps) break;
    iter++;
  }
  cout << "UpdateEps iter: " << iter << endl;
  cout << "eps: " << endl << eps << endl;
}

// need update in case that eps(k+1) < eps(k)
void nullModelClass::updateEpsOneStep()
{
  // the first eps is fixed at 0
  arma::vec d1eps(J-2, arma::fill::zeros);
  arma::mat d2eps(J-2, J-2, arma::fill::zeros);
  double temp1, temp2, temp3;
  
  for(int k = 1; k < J-1; k++){
    for(int i = 0; i < n; i++){
      temp1 = yMat(i,k)/muMat(i,k) - yMat(i,k+1)/muMat(i,k+1);
      temp2 = - yMat(i,k)/muMat(i,k)/muMat(i,k) - yMat(i,k+1)/muMat(i,k+1)/muMat(i,k+1);
      d1eps(k-1) += WMat(i,k) * temp1;
      d2eps(k-1, k-1) += WMat(i,k) * (1-2*nuMat(i,k)) * temp1 + WMat(i,k)*WMat(i,k)*temp2;
      if(k < J-2){
        temp3 = WMat(i,k) * WMat(i,k+1) * yMat(i,k+1)/muMat(i,k+1)/muMat(i,k+1);
        d2eps(k-1, k) += temp3;
        d2eps(k, k-1) += temp3;
      }
    }
  }
  
  arma::vec deps = -1 * inv(d2eps) * d1eps;
  for(int k = 1; k < J-1; k++){
    eps(k) += deps(k-1);
  }
}

void nullModelClass::updatePara(string excludechr)
{
  getPCGofSigmaAndCovaMat(CovaMat, iSigma_CovaMat, excludechr);
  arma::vec YVec = convert1(YMat);
  getPCGofSigmaAndVector(YVec, iSigma_YVec, excludechr); 
  
  // update beta
  arma::mat Cova_iSigma_CovaMat = CovaMat.t() * iSigma_CovaMat;
  arma::vec Cova_iSigma_YVec = CovaMat.t() * iSigma_YVec;
  beta = inv(Cova_iSigma_CovaMat) * Cova_iSigma_YVec;
  
  // update bVec
  arma::vec Z_iSigma_YVec = ZMat(iSigma_YVec);
  arma::vec Z_iSigma_Xbeta = ZMat(iSigma_CovaMat * beta);
  bVec = tau * getKinbVec(Z_iSigma_YVec - Z_iSigma_Xbeta, ptrGeno, excludechr);
}

// use PCG to calculate iSigma_xMat = Sigma^-1 %*% xMat
void nullModelClass::getPCGofSigmaAndCovaMat(arma::mat xMat,             // matrix with dim of n(J-1) x p
                                             arma::mat& iSigma_xMat,      // matrix with dim of n(J-1) x p
                                             string excludechr)
{
  int p1 = xMat.n_cols;
  for(int i = 0; i < p1; i++){
    arma::vec y1Vec = xMat.col(i);
    arma::vec iSigma_y1Vec = iSigma_xMat.col(i);
    getPCGofSigmaAndVector(y1Vec, iSigma_y1Vec, excludechr);
    iSigma_xMat.col(i) = iSigma_y1Vec;
  }
}

// use PCG to calculate xVec = Sigma^-1 %*% yVec
void nullModelClass::getPCGofSigmaAndVector(arma::vec y1Vec,    // vector with length of n(J-1)
                                            arma::vec& xVec,   // vector with length of n(J-1)
                                            string excludechr)
{
  arma::mat xMat = convert2(xVec);
  arma::mat y1Mat = convert2(y1Vec);
  // r2Vec and z2Vec are for the current step; r1Vec and z1Vec are for the previous step
  int iter = 0;
  arma::mat r2Mat = y1Mat - getSigmaxMat(xMat, excludechr);  // n x (J-1): r0 = y1Mat- Sigma %*% xMat
  double meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(n*(J-1));
  if(meanL2 <= tolPCG){
    // do nothing, xMat is already close to (Sigma)^-1 %*% y1Mat
  }else{
    iter++;
    
    arma::cube InvBlockDiagSigma = getInvBlockDiagSigma();
    arma::mat z2Mat = solverBlockDiagSigma(InvBlockDiagSigma, r2Mat);
    
    //
    arma::mat z1Mat, r1Mat;
    double beta1 = 0;
    arma::mat pMat = z2Mat;
    arma::mat ApMat = getSigmaxMat(pMat, excludechr);
    double alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
    xMat = xMat + alpha * pMat;
    r1Mat = r2Mat;
    z1Mat = z2Mat;
    r2Mat = r1Mat - alpha * ApMat;
    
    meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(n*(J-1));
    
    while (meanL2 > tolPCG && iter < maxiterPCG){
      iter++;
      
      //  z2Mat = minvMat % r2Mat;
      z2Mat = solverBlockDiagSigma(InvBlockDiagSigma, r2Mat);
      //
      beta1 = getInnerProd(z2Mat, r2Mat) / getInnerProd(z1Mat, r1Mat);
      pMat = z2Mat + beta1 * pMat;
      ApMat = getSigmaxMat(pMat, excludechr);
      alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
      xMat = xMat + alpha * pMat;
      r1Mat = r2Mat;
      z1Mat = z2Mat;
      r2Mat = r1Mat - alpha * ApMat;
      meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(n*(J-1));
    }
  }
  
  xVec = convert1(xMat);
  if (iter >= maxiterPCG){
    cout << "pcg did not converge. You may increase maxiter number." << endl;
  }
  cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
}

arma::mat nullModelClass::solverBlockDiagSigma(arma::cube& InvBlockDiagSigma,      // (J-1) x (J-1) x n
                                               arma::mat& xMat)                 // n x (J-1)
{
  arma::mat outMat(n, J-1);
  for(int i = 0; i < n; i++){
    outMat.row(i) = xMat.row(i) * InvBlockDiagSigma.slice(i); // could invert matrix?? be careful!
  }
  return(outMat);
}


// used in getPCGofSigmaAndVector()
arma::cube nullModelClass::getInvBlockDiagSigma()
{
  // get diagonal elements of GRM
  arma::vec* pDiagStdGeno = ptrGeno->getDiagStdGeno();
  arma::vec DiagGRM = tau * (*pDiagStdGeno) / M;   // n x 1
  
  double temp;
  //
  arma::cube InvBlockDiagSigma(J-1, J-1, n, arma::fill::zeros);
  for(int i = 0; i < n; i++){
    for(int j2 = 0; j2 < J-1; j2++){
      for(int j1 = 0; j1 < J-1; j1++){
        temp = iRMat(i,j2) * (1 / muMat(i,J-1)) * iRMat(i,j1) + DiagGRM(i);
        if(j2 == j1){
          temp += iRMat(i,j2) * (1 / muMat(i,j2)) * iRMat(i,j1); 
        }
        InvBlockDiagSigma(j2, j1, i) = temp;
      }
    }
    InvBlockDiagSigma.slice(i) = inv(InvBlockDiagSigma.slice(i));
  }
  return(InvBlockDiagSigma);
}

// convert from n x (J-1) matrix to n(J-1) vector
arma::vec nullModelClass::convert1(arma::mat xMat) // matrix with dim of n x (J-1)
{
  arma::vec xVec(n*(J-1));
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      xVec(index) = xMat(i,j);
      index++;
    }
  }
  return(xVec);
}

// convert n(J-1) vector to n x (J-1) matrix
arma::mat nullModelClass::convert2(arma::vec xVec) 
{
  arma::mat xMat(n,(J-1));
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      xMat(i,j) = xVec(index);
      index++;
    }
  }
  return(xMat);
}

// yMat1 = Sigma %*% xMat
arma::mat nullModelClass::getSigmaxMat(arma::mat xMat,   // matrix with dim of n x (J-1)
                                       string excludechr)
{
  arma::mat iR_xMat = iRMat % xMat;
  arma::mat iPsi_iR_xMat = getiPsixMat(iR_xMat);
  arma::mat yMat1 = iRMat % iPsi_iR_xMat;
  if(tau == 0){}
  else{
    arma::vec tZ_xMat = getRowSums(xMat);  // rowSums(xMat): n x 1
    arma::vec V_tZ_xMat = getKinbVec(tZ_xMat, ptrGeno, excludechr);
    yMat1.each_col() += tau * V_tZ_xMat;
  }
  return(yMat1);
}

arma::vec nullModelClass::getRowSums(arma::mat xMat)
{
  int n1 = xMat.n_rows;
  int n2 = xMat.n_cols;
  arma::vec y1Vec(n1, arma::fill::zeros);
  for(int i = 0; i < n1; i++){
    for(int j = 0; j < n2; j++){
      y1Vec(i) += xMat(i,j);
    }
  }
  return(y1Vec);
}


Rcpp::List nullModelClass::getNullModel()
{
  Rcpp::List outList = List::create(Named("N")=n,              // number of samples
                                    Named("M")=M,              // number of SNPs in Plink file
                                    Named("controlList")=controlList,
                                    Named("eta")=eta,          // X %*% beta + bVec
                                    Named("yVec")=yVec,        // matrix with dim of n x J: observation
                                    Named("Cova")=Cova,        // matrix with dim of n(J-1) x p: covariates
                                    Named("muMat")=muMat,      // matrix with dim of n x J: probability
                                    Named("YMat")=YMat,        // matrix with dim of n x (J-1): working variables
                                    Named("beta")=beta,        // parameter for covariates
                                    Named("bVec")=bVec,        // terms of random effect 
                                    Named("tau")=tau,          // variance component
                                    Named("eps")=eps,          // cutpoints
                                    Named("LOCOList")=LOCOList);         
  
  return(outList);
}


////////////////// main function //////////////

// [[Rcpp::export]]
Rcpp::List fitNullcpp(std::string Plink,
                      arma::vec posSampleInPlink,
                      arma::mat CovaR,
                      arma::Col<int> yVecR, // should be from 1 to J
                      arma::vec betaR,
                      arma::vec bVecR,
                      arma::vec epsR,
                      double tauR,
                      arma::mat GMatRatioR,
                      Rcpp::List controlListR)
{
  // read in genotype from Plink file for GRM
  genoClass obj;
  float memoryChunk = controlListR["memoryChunk"];
  obj.setGenoObj(Plink, posSampleInPlink, memoryChunk);
  genoClass* ptrGenoInput = &obj;
  
  // Null model fitting
  nullModelClass objm;
  objm.setNullModel(CovaR, yVecR, ptrGenoInput, betaR,  bVecR,  epsR,  tauR, GMatRatioR, controlListR);
  objm.fitNullModel();
  
  Rcpp::List outList = objm.getNullModel();
  return(outList);
}


