
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include <RcppArmadillo.h>

#include <vector>
#include <cmath>
#include <sys/time.h>
#include <RcppArmadilloExtensions/sample.h> // sample

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

////////////////// sub-functions //////////////

// yMat: matrix with dim of n x J
// [[Rcpp::export]]
arma::mat getyMatR(arma::mat yVec, int n, int J)
{
  arma::mat yMat(n, J, arma::fill::zeros);
  for(int i = 0; i < n; i++)
    yMat(i, yVec(i)-1) = 1;
  return(yMat);
}

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat)
{
  double innerProd = arma::accu(x1Mat % x2Mat);
  return(innerProd);
}

arma::vec nb(int n){
  return(Rcpp::rbinom(n,1,0.5));
}

double calCV(arma::vec xVec){
  int n = xVec.size();
  double Mean = arma::mean(xVec);
  double Sd = arma::stddev(xVec);
  double CV = (Sd/Mean)/n;
  return(CV);
}

arma::vec getTime(){
  arma::vec Time(2, arma::fill::zeros);
  struct timeval time;
  if(gettimeofday(&time,NULL)){
    Time(0) = 0;
  }else{
    Time(0) = (double)time.tv_sec + (double)time.tv_usec * .000001;
  }
  Time(1) = (double)clock() / CLOCKS_PER_SEC;
  return Time; 
}

void printTime(arma::vec t1, arma::vec t2, std::string message){
  double wallTime = t2(0) - t1(0);
  double cpuTime = t2(1) - t1(1);
  if(cpuTime < 60){
    Rprintf ("It took %f seconds (%f CPU seconds) to %s.\n", 
             wallTime, cpuTime, message.c_str());
  }else if(cpuTime < 3600){
    Rprintf ("It took %f minutes (%f CPU minutes) to %s.\n", 
             wallTime/60, cpuTime/60, message.c_str());
  }else{
    Rprintf ("It took %f hours (%f CPU hours) to %s.\n", 
             wallTime/3600, cpuTime/3600, message.c_str());
  }
}

// C++ version of which(). Note: start from 0, not 1 
std::vector<int> whichCPP(Rcpp::StringVector strVec, 
                          std::string strValue)
{
  std::vector<int> indexVec;
  for(int i = 0; i < strVec.size(); i++){
    if(std::string(strVec(i))==strValue)
      indexVec.push_back(i);
  }
  return(indexVec);
}

// outMat = PsiMat %*% xMat, PsiMat is determined by muMat
arma::mat getPsixMat(arma::mat xMat,    // matrix: n x (J-1)
                     arma::mat muMat)   // matrix: n x J
{
  int n = muMat.n_rows;
  int J = muMat.n_cols;
  
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

// duplicate each row for (J-1) times: n x p -> n(J-1) x p
arma::mat getCovaMat(arma::mat Cova,   // matrix: n x p
                     int n, int J, int p)      
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

// convert: n x (J-1) -> n(J-1) x 1
arma::vec convert1(arma::mat xMat, // matrix: n x (J-1)
                   int n, int J) 
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

// convert: n(J-1) x 1 -> n x (J-1)
arma::mat convert2(arma::vec xVec, // n(J-1) x 1 
                   int n, int J)
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

// sum each (J-1) cols to 1 col: p x n(J-1) -> p x n (OR) p x (J-1) -> p x 1
arma::mat sumCols(arma::mat inMat,
                  int J)
{
  int n = inMat.n_cols / (J-1);
  int p = inMat.n_rows;
  arma::mat outMat(p, n, arma::fill::zeros);
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      outMat.col(i) += inMat.col(index);
      index++;
    }
  }
  return(outMat);
}

// get RPsiP: (J-1) x (J-1) x n 
// Only used in getVarWFast(): 
arma::vec getRPsiR(arma::mat muMat,
                   arma::mat iRMat,
                   int n, int J, int p)   
{
  // arma::cube RPsiR(J-1, J-1, n, arma::fill::zeros);
  arma::vec RPsiRVec(n, arma::fill::zeros);
  arma::mat muRMat = muMat.cols(0, J-2) / iRMat;
  for(int i = 0; i < n; i++){
    for(int j1 = 0; j1 < J-1; j1++){
      // RPsiR(j1,j1,i) += muRMat(i,j1) / iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      RPsiRVec(i) += muRMat(i,j1) / iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      for(int j2 = j1+1; j2 < J-1; j2++){
        // RPsiR(j1,j2,i) -= muRMat(i,j1) * muRMat(i,j2);
        RPsiRVec(i) -= 2* muRMat(i,j1) * muRMat(i,j2);
      }
    }
  }
  // return(RPsiR);
  return(RPsiRVec);
}

// get a list for p value calculation in step 2
// [[Rcpp::export]]
Rcpp::List getobjP(arma::mat Cova,     // matrix: n x p
                   arma::mat yMat,
                   arma::mat muMat,    // matrix: n x J
                   arma::mat iRMat)    // matrix: n x (J-1)
{
  int n = muMat.n_rows;
  int J = muMat.n_cols;
  int p = Cova.n_cols;
  
  // output for Step 2
  arma::mat XR_Psi_R(p, n*(J-1));                // p x n(J-1)
  arma::mat CovaMat = getCovaMat(Cova, n, J, p); // n(J-1) x p
  for(int k = 0; k < p; k++){
    arma::mat xMat = convert2(CovaMat.col(k), n, J);
    arma::vec temp = convert1(getPsixMat(xMat / iRMat, muMat) / iRMat, n, J);
    XR_Psi_R.row(k) = temp.t();
  }
  // arma::mat XXR_Psi_RX = CovaMat * inv(XR_Psi_R * CovaMat);             // (n(J-1) x p) * (p x p) = n(J-1) x p
  arma::mat XXR_Psi_RX_new = Cova * inv(XR_Psi_R * CovaMat);               // (n x p) * (p x p) = n x p
  
  // sum each (J-1) rows to 1 row: p x n(J-1) -> p x n
  arma::mat XR_Psi_R_new = sumCols(XR_Psi_R, J);      // p x n
  arma::mat ymuMat = yMat - muMat;                    // n x J
  arma::mat RymuMat = ymuMat.cols(0, J-2) / iRMat;    // n x (J-1): R %*% (y - mu)
  arma::mat RymuVec = sumCols(RymuMat, J);            // n x 1
  // arma::cube RPsiR = getRPsiR(muMat, iRMat, n, J, p); // (J-1) x (J-1) x n 
  arma::vec RPsiR = getRPsiR(muMat, iRMat, n, J, p); // (J-1) x (J-1) x n 
  
  Rcpp::List objP = List::create(Named("n")=n,
                                 Named("J")=J,
                                 Named("p")=p,
                                 Named("XXR_Psi_RX_new")=XXR_Psi_RX_new,
                                 Named("XR_Psi_R_new")=XR_Psi_R_new,           
                                 Named("RymuVec")=RymuVec,
                                 Named("RPsiR")=RPsiR,
                                 Named("muMat")=muMat,
                                 Named("iRMat")=iRMat);
  return(objP);
}

double getVarWFast(arma::vec adjGVec,  // n x 1
                   arma::vec RPsiRVec, // n x 1
                   // arma::mat adjGMat,  // n x (J-1)
                   // arma::cube RPsiR,   // (J-1) x (J-1) x n
                   int n, int J)
{
  double VarW = 0;
  for(int i = 0; i < n; i++){
    VarW += RPsiRVec(i) * adjGVec(i) * adjGVec(i);
    // for(int j1 = 0; j1 < J-1; j1++){
      // VarW += RPsiR(j1,j1,i) * adjGMat(i,j1) * adjGMat(i,j1);
      // for(int j2 = j1+1; j2 < J-1; j2++){
      //   VarW += 2 * RPsiR(j1,j2,i) * adjGMat(i,j1) * adjGMat(i,j2);
      // }
    // }
  }
  return(VarW);
}

double getStatFast(arma::vec GVec,         // n x 1
                   arma::vec RymuVec,      // n x 1: row sum of the n x (J-1) matrix R %*% (yMat - muMat)
                   int n)
{
  double Stat = 0;
  for(int i = 0; i < n; i++){
    if(GVec(i) != 0){
      Stat += GVec(i) * RymuVec(i);
    }
  }
  return(Stat);
}

arma::vec getadjGFast(arma::vec GVec,
                      arma::mat XXR_Psi_RX_new,   // XXR_Psi_RX_new ( n x p )
                      arma::mat XR_Psi_R_new,     // XR_Psi_R_new ( p x n ), sum up XR_Psi_R ( p x n(J-1) ) for each subject 
                      int n, int J, int p)
{
  // arma::mat adjGMat(n, J-1);
  // arma::vec adjGVec(n);
  
  // To increase computational efficiency when lots of GVec elements are 0
  arma::vec XR_Psi_RG1(p, arma::fill::zeros);
  for(int i = 0; i < n; i++){
    if(GVec(i) != 0){
      XR_Psi_RG1 += XR_Psi_R_new.col(i) * GVec(i);
    }
  }
  
  arma::vec adjGVec = GVec - XXR_Psi_RX_new * XR_Psi_RG1;
  // arma::vec TempVec = XXR_Psi_RX_new * XR_Psi_RG1;   // n x 1
  // int index = 0;
  // for(int i = 0; i < n; i++){
    // adjGVec(i) = GVec(i) - TempVec(i);
    // for(int j = 0; j < J-1; j++){
    //   adjGMat(i,j) = GVec(i) - TempVec(index);
    //   index++;
    // }
  // }
  // return(adjGMat);
  return(adjGVec);
}

// [[Rcpp::export]]
Rcpp::List outputadjGFast(arma::vec GVec,
                          Rcpp::List objP)
{
  arma::vec adjGVec = getadjGFast(GVec, objP["XXR_Psi_RX_new"], objP["XR_Psi_R_new"], objP["n"], objP["J"], objP["p"]);
  double Stat = getStatFast(adjGVec, objP["RymuVec"], objP["n"]);
  double VarW = getVarWFast(adjGVec, objP["RPsiR"], objP["n"], objP["J"]);
  Rcpp::List outList = List::create(Named("adjGVec")=adjGVec,
                                    Named("Stat")=Stat,           
                                    Named("VarW")=VarW);
  
  return(outList);
}


arma::Mat<int> makeChrIdx(string excludeChr, vector<string> chrVecParallel, int M){
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

////////////////// genoClass //////////////

class genoClass{
private:
  
  // PLINK format
  const static unsigned char HOM_REF = 0x3;  // 0b11 ;
  const static unsigned char HET = 0x2;      // 0b10 ;
  const static unsigned char HOM_ALT = 0x0;  // 0b00 ;
  const static unsigned char MISSING = 0x1;  // 0b01 ;
  
  // Basic information
  int M;                         // number of markers
  int N,N_Old;                   // number of samples, after and before filt
  int numArrays;                 // make arrays to store genotype data to avoid large continuous memory usage
  int numMarkersofEachArray;
  int numMarkersofLastArray;
  int numBytesofEachMarker;
  int numBytesofEachMarker_Old; 
  double numBytesReserve;   // unit: Gb
  
  // input file stream of .bed file
  ifstream ibedfile;  
  int bufferG1;
  unsigned char bufferG4;
  arma::vec OneMarkerG1;
  arma::vec OneMarkerG1_Old;                   // a double vector: each element (1 byte) for 4 genotypes (4 samples)
  std::vector<unsigned char> OneMarkerG4_Old;  // a char vector: each element (1 byte) for 4 genotypes (4 samples)
  arma::vec DiagStdGeno;
  
  std::vector< std::vector<unsigned char>* > genoVecofPointers;
  
  // set environment
  void openPlink(string, arma::ivec);
  void setChrVec(string, Rcpp::StringVector&, vector<string>&);
  void setArrays(float);
  int readFile(string, string);
  double getOneMarkerPlink(vector<unsigned char>&, arma::ivec);
  void setOneMarkerArray(int);
  double getinvStd(double);
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
  void getGenotype(unsigned char* c, const int pos, int& geno) {
    geno = ((*c) >> (pos << 1)) & 0b11; 
  }
  
  void setStdGenoLookUpArr(double maf, double invsd, arma::vec& stdGenoLookUpArr){
    double maf2 = 2 * maf;
    stdGenoLookUpArr(0) = (2-maf2)*invsd; // HOM_ALT = 0x0;
    stdGenoLookUpArr(1) = 0;              // MISSING = 0x1;
    stdGenoLookUpArr(2) = (1-maf2)*invsd; // HET = 0x2;
    stdGenoLookUpArr(3) = (0-maf2)*invsd; // HOM_REF = 0x3;
  }

public:
  void setGenoObj(string, arma::ivec, double);
  int getN(){return N;}
  int getM(){return M;}
  arma::vec getalleleFreqVec(){return alleleFreqVec;}
  void getOneMarkerStd(size_t, arma::vec*);
  void getOneMarker(size_t, arma::vec*);
  arma::vec* getDiagStdGeno(){return &DiagStdGeno;}
  arma::vec getAvec(){return alleleFreqVec;}
  Rcpp::StringVector getchrVec(){return chrVec;}
  vector<string> getchrVecParallel(){return chrVecParallel;}  // the same as chrVec, different data type
};


void genoClass::setGenoObj(string Plink, 
                           arma::ivec posSampleInPlink, 
                           double memoryChunk)    // unit is Gb
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
  cout << "numBytesReserve:\t" << numBytesReserve << " Gb" << endl << endl;
  
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
  double freq, invStd;
  arma::vec t1  = getTime();
  for(int m=0; m<M; m++)
  {
    ibedfile.seekg(3+numBytesofEachMarker_Old*m);
    ibedfile.read((char*)(&OneMarkerG4_Old[0]), numBytesofEachMarker_Old);
    freq = getOneMarkerPlink(OneMarkerG4_Old, posSampleInPlink);
    setOneMarkerArray(m);
    invStd = getinvStd(freq);
    alleleFreqVec[m]=freq;
    invStdVec[m]=invStd;
    if((m+1) % 10000 == 0){
      cout << "Complete\t" << m+1 <<"\tSNPs!!!!" << endl;
      cout << "Freq\t" << freq << endl;
    }
  }
  
  arma::vec t2  = getTime();
  printTime(t1, t2, "read Plink files");
  setDiagStdGeno();
  ibedfile.close();
}

double genoClass::getinvStd(double freq)
{
  double Std = sqrt(2*freq*(1-freq));
  if(Std == 0)
    return 0;
  else
    return 1/Std;
}

// oneMarkerG1 --> genoVecofPointers
void genoClass::setOneMarkerArray(int m)
{
  // loop for all samples
  int whichArray = m / numMarkersofEachArray;

  bufferG4 = 0;
  for(int ind=0; ind<N; ind++){
    int posInByte = ind % 4;
    bufferG1 = OneMarkerG1[ind];
    setGenotype(&bufferG4, posInByte, bufferG1);
    if((posInByte==3) | (ind==(N-1))){
      genoVecofPointers[whichArray]->push_back(bufferG4); // avoid large continuous memory usage
      bufferG4 = 0;
    }
  }
}

// OneMarkerG4_Old --> OneMarkerG1
double genoClass::getOneMarkerPlink(std::vector<unsigned char>& OneMarkerG4_Old,   
                                   arma::ivec posSampleInPlink)              // position of samples in model (start from 1)
{
  // read in genotypes of one marker
  int sum = 0;
  std::vector<int> indexNA;
  for(int index=0; index<N; index++){
    int ind = posSampleInPlink[index] - 1;
    bufferG4 = OneMarkerG4_Old[ind/4];
    getGenotype(&bufferG4, ind%4, bufferG1); // bufferG4 -> bufferG1
    switch(bufferG1){
    case HOM_REF: break;
    case HET: sum+=1; break;
    case HOM_ALT: sum+=2; break;
    case MISSING: indexNA.push_back(index); break;
    }
    OneMarkerG1[index] = bufferG1;
  }
  
  // calculate freq, that is, twice minor allele frequency (MAF)
  int lengthNA = indexNA.size();
  int count = N - lengthNA;
  double freq = (double)sum/(double)count;
  return freq/2;
}

void genoClass::setArrays(float memoryChunk)
{
  numBytesofEachMarker = (N+3)/4;
  numBytesofEachMarker_Old = (N_Old+3)/4;
  numBytesReserve = (double)(numBytesofEachMarker+2) * M / pow(10.0, 9.0);    // not sure about why M*2 is needed (2 bytes for each marker)
  numMarkersofEachArray = floor((memoryChunk*pow(10.0, 9.0))/numBytesofEachMarker);
  if(numMarkersofEachArray > M)
    numMarkersofEachArray = M;
  numArrays = (M-1) / numMarkersofEachArray + 1;
  numMarkersofLastArray = M - (numArrays - 1) * numMarkersofEachArray;
  
  // get binray data from plink files
  OneMarkerG4_Old.reserve(numBytesofEachMarker_Old);  
  OneMarkerG4_Old.resize(numBytesofEachMarker_Old);
  OneMarkerG1_Old.zeros(N_Old);             
  OneMarkerG1.zeros(N);
  
  // set genoVecofPointers
  genoVecofPointers.resize(numArrays);
  for (int i = 0; i < numArrays-1 ; i++){
    genoVecofPointers[i] = new vector<unsigned char>;
    genoVecofPointers[i]->reserve(numMarkersofEachArray*numBytesofEachMarker);
  }
  genoVecofPointers[numArrays-1] = new vector<unsigned char>;
  genoVecofPointers[numArrays-1]->reserve(numMarkersofLastArray*numBytesofEachMarker);
  
  alleleFreqVec.zeros(M);
  invStdVec.zeros(M);
}

void genoClass::openPlink(string Plink, arma::ivec posSampleInPlink)
{
  string bimfile = Plink + ".bim";
  string famfile = Plink + ".fam";
  string bedfile = Plink + ".bed";
  N_Old = readFile(famfile, "Error! fam file not open!");
  M = readFile(bimfile, "Error! bim file not open!");
  N = posSampleInPlink.size();
  
  if((max(posSampleInPlink)>N_Old) | (min(posSampleInPlink)<1)) 
    stop("posSampleInPlink should be between 1 and N");
  
  ibedfile.open(bedfile.c_str(), ios::binary);
  if (!ibedfile.is_open())
    stop("Error! bed file not open!");
  
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

int genoClass::readFile(string file, string errInfo){
  ifstream ifile;
  ifile.open(file.c_str());
  if(!ifile.is_open())
    stop(errInfo);
  
  string junk;
  int count=0;
  while(getline(ifile,junk))
    count++;
  
  ifile.close();
  return count;
}

void genoClass::getOneMarkerStd(size_t m, arma::vec* oneMarker){
  // avoid large continuous memory usage
  int whichArray = m / numMarkersofEachArray;
  int posMarker = m % numMarkersofEachArray;
  
  // set up start byte index and end byte index
  int startBtIdx = numBytesofEachMarker * posMarker;
  int endBtIdx = startBtIdx + numBytesofEachMarker;
  
  // use bufferG4b to replace bufferG4 in case of parallele computation
  unsigned char bufferG4b;
  // unsigned char bufferG1b;
  int bufferG1b;
  
  arma::vec stdGenoLookUpArr(4);
  setStdGenoLookUpArr(alleleFreqVec[m], invStdVec[m], stdGenoLookUpArr);
  std::vector<unsigned char>* genoPtr = genoVecofPointers[whichArray];
  int ind = 0;
  for(int BtIdx = startBtIdx; BtIdx < endBtIdx; BtIdx++){
    bufferG4b = genoPtr->at(BtIdx); // unsigned char: 4 markers
    for(unsigned char posInByte = 0; (posInByte<4) & (ind<N); posInByte++,ind++){
      bufferG1b = (bufferG4b >> (posInByte << 1)) & 0b11;
      oneMarker->at(ind) = stdGenoLookUpArr(bufferG1b);
    }
  } 
}

void genoClass::getOneMarker(size_t m, arma::vec* oneMarker){
  // avoid large continuous memory usage
  int whichArray = m / numMarkersofEachArray;
  int posMarker = m % numMarkersofEachArray;
  
  // set up start byte index and end byte index
  int startBtIdx = numBytesofEachMarker * posMarker;
  int endBtIdx = startBtIdx + numBytesofEachMarker;
  
  // use bufferG4b to replace bufferG4 in case of parallele computation
  unsigned char bufferG4b;
  int bufferG1b;
  double bufferG1c;
  
  int ind = 0;
  for(int BtIdx = startBtIdx; BtIdx < endBtIdx; BtIdx++){
    bufferG4b = genoVecofPointers[whichArray]->at(BtIdx); //avoid large continuous memory usage
    for(int posInByte = 0; (posInByte<4) & (ind<N); posInByte++,ind++){
      getGenotype(&bufferG4b, posInByte, bufferG1b);
      switch(bufferG1b){
      case HOM_REF: bufferG1c = 0; break;
      case HET: bufferG1c = 1; break;
      case HOM_ALT: bufferG1c = 2; break;
      case MISSING: bufferG1c = alleleFreqVec(m); break;
      }
      oneMarker->at(ind) = bufferG1c;
    }
  } 
}


void genoClass::setDiagStdGeno(){
  DiagStdGeno.zeros(N);
  arma::vec oneMarker(N);
  for(int m=0; m<M; m++){
    getOneMarkerStd(m, &oneMarker);  // write Standard Genotype to OneMarkerG1
    DiagStdGeno = DiagStdGeno + (oneMarker) % (oneMarker);
  }
}


//http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

genoClass objG;

//http://gallery.rcpp.org/articles/parallel-inner-product/
struct getKinbVecParallel : public Worker
{
  // source vectors
  arma::vec& bVec;
  unsigned int N;
  unsigned int M;
  
  // product that I have accumulated
  arma::vec KinbVec;
  int counts;
  
  // constructors
  getKinbVecParallel(arma::vec& bVec)
    : bVec(bVec), counts(0) 
  {
    M = objG.getM();
    N = objG.getN();
    KinbVec.zeros(N);
  }
  getKinbVecParallel(const getKinbVecParallel& getKinbVecParallel, Split)
    : bVec(getKinbVecParallel.bVec), counts(0)
  {
    N = getKinbVecParallel.N;
    M = getKinbVecParallel.M;
    KinbVec.zeros(getKinbVecParallel.N);
  }
  // process just the elements of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    arma::vec oneMarker(N);
    for(unsigned int m = begin; m < end; m++){
      objG.getOneMarkerStd(m, &oneMarker);
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


arma::vec getKinbVec(arma::vec& bVec, genoClass* ptrGeno, string excludeChr) {
  
  size_t M = ptrGeno->getM();
  vector<string> chrVecParallel = ptrGeno->getchrVecParallel();
  
  // declare the InnerProduct instance that takes a pointer to the vector data
  getKinbVecParallel getKinbVecParallel(bVec);
  
  arma::Mat<int> ChrIdx = makeChrIdx(excludeChr, chrVecParallel, M);
  // cout << "ChrIdx:" << endl << ChrIdx << endl;
  
  int nIdx = ChrIdx.n_rows;
  for(int i = 0; i < nIdx; i++){
    int idxStart = ChrIdx(i,0);
    int idxEnd = ChrIdx(i,1);
    // cout << idxStart << "\t" << idxEnd << endl;
    // call paralleReduce to start the work
    parallelReduce(idxStart, idxEnd, getKinbVecParallel);
  }
  
  return getKinbVecParallel.KinbVec/getKinbVecParallel.counts;
}

// //http://gallery.rcpp.org/articles/parallel-inner-product/
// struct getKinbVecParallel : public Worker
// {
//   // source vectors
//   arma::vec& bVec;
//   unsigned int N;
//   unsigned int M;
//   genoClass* ptrGeno;
//   
//   // product that I have accumulated
//   arma::vec KinbVec;
//   int counts;
//   
//   // constructors
//   getKinbVecParallel(arma::vec& bVec, genoClass* ptrGeno)
//     : bVec(bVec), ptrGeno(ptrGeno), counts(0) 
//   {
//     M = ptrGeno->getM();
//     N = ptrGeno->getN();
//     KinbVec.zeros(N);
//   }
//   getKinbVecParallel(const getKinbVecParallel& getKinbVecParallel, Split)
//     : bVec(getKinbVecParallel.bVec), ptrGeno(getKinbVecParallel.ptrGeno), counts(0)
//   {
//     N = getKinbVecParallel.N;
//     M = getKinbVecParallel.M;
//     KinbVec.zeros(getKinbVecParallel.N);
//   }
//   // process just the elements of the range I've been asked to
//   void operator()(std::size_t begin, std::size_t end) {
//     arma::vec oneMarker(N);
//     for(unsigned int m = begin; m < end; m++){
//       ptrGeno->getOneMarkerStd(m, &oneMarker);
//       KinbVec += oneMarker * arma::dot(oneMarker, bVec);
//       counts ++;
//     }
//   }
//   
//   // join my value with that of another InnerProduct
//   void join(const getKinbVecParallel & rhs) {
//     KinbVec += rhs.KinbVec;
//     counts += rhs.counts;
//   }
// };
// 
// 
// arma::vec getKinbVec(arma::vec bVec, genoClass* ptrGeno, string excludeChr) {
//   
//   size_t M = ptrGeno->getM();
//   vector<string> chrVecParallel = ptrGeno->getchrVecParallel();
//   
//   // declare the InnerProduct instance that takes a pointer to the vector data
//   getKinbVecParallel getKinbVecParallel(bVec, ptrGeno);
//   
//   arma::Mat<int> ChrIdx = makeChrIdx(excludeChr, chrVecParallel, M);
//   // cout << "ChrIdx:" << endl << ChrIdx << endl;
//   
//   int nIdx = ChrIdx.n_rows;
//   for(int i = 0; i < nIdx; i++){
//     int idxStart = ChrIdx(i,0);
//     int idxEnd = ChrIdx(i,1);
//     // cout << idxStart << "\t" << idxEnd << endl;
//     // call paralleReduce to start the work
//     parallelReduce(idxStart, idxEnd, getKinbVecParallel);
//   }
//   
//   return getKinbVecParallel.KinbVec/getKinbVecParallel.counts;
// }

////////////////// nullModelClass ////////////// 

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
  arma::mat getyMat();
  void getTraceRandMat();
  arma::vec ZMat(arma::vec), tZMat(arma::vec);
  
  // functions in fitNullModel()
  void updateMats();
  arma::mat getiPsixMat(arma::mat);
  arma::vec getRowSums(arma::mat);
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
  Rcpp::List getLOCO(string, std::vector<int>);
  arma::mat getVarRatio(std::vector<int>, string);
  arma::rowvec getVarOneSNP(arma::vec, string, Rcpp::List);
  double getVarP(arma::vec, string);
  
  // functions if LOCO = FALSE
  arma::mat getVarRatio(arma::mat);
  
public:
  
  void setNullModel(arma::mat, arma::Col<int>, genoClass*, arma::vec, arma::vec, arma::vec, double, arma::mat, Rcpp::List);
  void fitNullModel();
  Rcpp::List getNullModel();
  
};

void nullModelClass::fitNullModel()
{
  // initial vector
  arma::vec t1  = getTime();
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
    // output variance matrix
    arma::mat VarRatioMat = getVarRatio(GMatRatio);
    double VarRatio = arma::mean(VarRatioMat.col(4));
    
    Rcpp::List temp = List::create(Named("muMat")=muMat,
                                   Named("iRMat")=iRMat,
                                   Named("VarRatioMat")=VarRatioMat,
                                   Named("VarRatio")=VarRatio);
    LOCOList["LOCO=F"] = temp;
  }
  // 
  arma::vec t2  = getTime();
  printTime(t1, t2, "fit the null model");
}


double nullModelClass::getVarP(arma::vec adjGVec,
                               string excludechr)
{
  arma::vec adjGVecLong = tZMat(adjGVec);
  arma::vec iSigmaGVec(n*(J-1), arma::fill::zeros);
  getPCGofSigmaAndVector(adjGVecLong, iSigmaGVec, excludechr);
  double VarP = as_scalar(adjGVecLong.t() * (iSigmaGVec - iSigmaX_XSigmaX * (CovaMat.t() * iSigmaGVec)));
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
  arma::vec adjGVec = adjGList["adjGVec"];
  double Stat = adjGList["Stat"];
  double VarW = adjGList["VarW"];
  double VarP = getVarP(adjGVec, excludechr);
  
  VarOut(0) = AF;
  VarOut(1) = Stat;
  VarOut(2) = VarW;
  VarOut(3) = VarP;
  VarOut(4) = VarP/VarW;
  return(VarOut);
}

arma::mat nullModelClass::getVarRatio(arma::mat GMatRatio)
{
  Rcpp::List objP = getobjP(Cova, yMat, muMat, iRMat);
  
  arma::vec GVec(n);
  arma::rowvec VarOneSNP(5);
  
  arma::mat VarRatioMat(nSNPsVarRatio, 5);
  arma::mat newVarRatio(10, 5);
  
  int index = 0;
  int indexTot = 0;
  while(index < nSNPsVarRatio){
    GVec = GMatRatio.col(index);
    VarOneSNP = getVarOneSNP(GVec, "n", objP);
    VarRatioMat.row(index) = VarOneSNP;
    index++;
    indexTot++;
  }
  
  arma::vec VarRatio = VarRatioMat.col(4);
  double CV = calCV(VarRatio);
  cout << "nSNPs for CV: " << index << endl;
  cout << "CV: " << CV << endl;
  
  while(CV > CVcutoff && VarRatioMat.n_rows <= 100){
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
                                      string excludechr)
{
  Rcpp::List objP = getobjP(Cova, yMat, muMat, iRMat);
  arma::vec alleleFreqVec = ptrGeno->getalleleFreqVec();
  
  arma::vec GVec(n);
  arma::rowvec VarOneSNP(5);
  
  arma::mat VarRatioMat(nSNPsVarRatio, 5);
  arma::mat newVarRatio(10, 5);
  
  int index = 0;
  int indexTot = 0;
  while(index < nSNPsVarRatio){
    int m = indexSNPs[indexTot];
    indexTot++;
    if(alleleFreqVec(m) > minMafVarRatio && alleleFreqVec(m) < 1-minMafVarRatio){
      ptrGeno->getOneMarker(m, &GVec);
      VarOneSNP = getVarOneSNP(GVec, excludechr, objP);
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
        ptrGeno->getOneMarker(m, &GVec);
        VarOneSNP = getVarOneSNP(GVec, excludechr, objP);
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
  
  // output variance matrix
  arma::mat VarRatioMat = getVarRatio(indexSNPs, excludechr);
  double VarRatio = arma::mean(VarRatioMat.col(4));
  
  Rcpp::List outLOCO = List::create(Named("muMat")=muMat,
                                    Named("iRMat")=iRMat,
                                    Named("VarRatioMat")=VarRatioMat,
                                    Named("VarRatio")=VarRatio);
  return(outLOCO);
}

void nullModelClass::setNullModel(arma::mat CovaR,
                                  arma::Col<int> yVecR,     // should be from 1 to J
                                  genoClass* ptrGenoInput,
                                  arma::vec betaR,
                                  arma::vec bVecR,
                                  arma::vec epsR,           // 
                                  double tauR,
                                  arma::mat GMatRatioR,
                                  Rcpp::List controlListR)
{
  n = CovaR.n_rows;
  p = CovaR.n_cols;
  J = max(yVecR);
  M = ptrGenoInput->getM();
  
  Cova = CovaR;
  CovaMat = getCovaMat(Cova, n, J, p);
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


// yMat: matrix with dim of n x J
arma::mat nullModelClass::getyMat()
{
  arma::mat yMat(n, J, arma::fill::zeros);
  for(int i = 0; i < n; i++)
    yMat(i, yVec(i)-1) = 1;
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
  arma::vec t1  = getTime();
  for(int itrace = 0; itrace < tracenrun; itrace++){
    
    arma::vec uVec = nb(n*(J-1));
    uVec = uVec*2 - 1;
    TraceRandMat.col(itrace) = uVec;
    arma::vec ZuVec = ZMat(uVec);
    V_TRM.col(itrace) = tZMat(getKinbVec(ZuVec, ptrGeno, "n"));
  }
  
  arma::vec t2  = getTime();
  printTime(t1, t2, "calculate 30 getKinbVec");
}

// sum up each (J-1) elements: n(J-1) x 1 -> n x 1
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

// duplicate each element for (J-1) times: n x 1 -> n(J-1) x 1 
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
  cout << "Start updating tau..." << endl;
  arma::vec YVec = convert1(YMat, n, J);
  getPCGofSigmaAndCovaMat(CovaMat, iSigma_CovaMat, "n");
  getPCGofSigmaAndVector(YVec, iSigma_YVec, "n"); 
  iSigmaX_XSigmaX = iSigma_CovaMat * inv(CovaMat.t() * iSigma_CovaMat);
  arma::vec PYVec = iSigma_YVec - iSigmaX_XSigmaX * (CovaMat.t() * iSigma_YVec);
  arma::vec ZPYVec = ZMat(PYVec);
  arma::vec VPYVec = tZMat(getKinbVec(ZPYVec, ptrGeno, "n"));
  
  getPCGofSigmaAndVector(VPYVec, iSigma_VPYVec, "n");
  arma::vec PVPYVec = iSigma_VPYVec - iSigmaX_XSigmaX * (CovaMat.t() * iSigma_VPYVec);
  double YPVPY = as_scalar(YVec.t() * PVPYVec);
  double YPVPVPY = as_scalar(VPYVec.t() * PVPYVec);
  // The below is to calculate trace
  getPCGofSigmaAndCovaMat(V_TRM, iSigma_V_TRM, "n");
  double tracePV = 0;
  int m = TraceRandMat.n_cols;
  for(int i = 0; i < m; i++){
    arma::vec iSigma_V_TRM_col = iSigma_V_TRM.col(i);
    arma::vec P_V_TRM_col = iSigma_V_TRM_col - iSigmaX_XSigmaX * (CovaMat.t() * iSigma_V_TRM_col);
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
    cout << "diffBeta:\t" << diffBeta << endl << endl;
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
    
    diffeps = max(abs(eps - eps0)/(abs(eps) + abs(eps0) + tolBeta));
    iter++;
    // if(diffeps < tolEps) break;
    if(diffeps < tolEps){
      eps = eps0;
      break;
    }
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
  arma::vec YVec = convert1(YMat, n, J);
  getPCGofSigmaAndVector(YVec, iSigma_YVec, excludechr); 
  
  // update beta
  arma::mat XSigmaX = inv(CovaMat.t() * iSigma_CovaMat);
  arma::vec Cova_iSigma_YVec = CovaMat.t() * iSigma_YVec;
  beta = XSigmaX * Cova_iSigma_YVec;
  iSigmaX_XSigmaX = iSigma_CovaMat * XSigmaX;
  
  // update bVec
  arma::vec Z_iSigma_YVec = ZMat(iSigma_YVec);
  arma::vec Z_iSigma_Xbeta = ZMat(iSigma_CovaMat * beta);
  arma::vec tempVec = Z_iSigma_YVec - Z_iSigma_Xbeta;
  bVec = tau * getKinbVec(tempVec, ptrGeno, excludechr);
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
  arma::mat xMat = convert2(xVec, n, J);
  arma::mat y1Mat = convert2(y1Vec, n, J);
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
  
  xVec = convert1(xMat, n, J);
  Rcpp::checkUserInterrupt();
  if (iter >= maxiterPCG){
    cout << "pcg did not converge. You may increase maxiter number." << endl;
  }
  cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
}

arma::mat nullModelClass::solverBlockDiagSigma(arma::cube& InvBlockDiagSigma,   // (J-1) x (J-1) x n
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



// yMat1 = Sigma %*% xMat
arma::mat nullModelClass::getSigmaxMat(arma::mat xMat,   // matrix: n x (J-1) 
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

// sum for each row: n1 x n2 matrix -> n1 x 1 vector
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
                                    Named("yVec")=yVec,        // matrix with dim of n x 1: observation
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
                      arma::ivec posSampleInPlink,
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
  // genoClass objG;
  float memoryChunk = controlListR["memoryChunk"];
  objG.setGenoObj(Plink, posSampleInPlink, memoryChunk);
  genoClass* ptrGenoInput = &objG;
  
  // Null model fitting
  nullModelClass objM;
  objM.setNullModel(CovaR, yVecR, ptrGenoInput, betaR,  bVecR,  epsR,  tauR, GMatRatioR, controlListR);
  // objM.fitNullModel();
  
  Rcpp::List outList = objM.getNullModel();
  return(outList);
}


