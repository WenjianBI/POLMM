// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "util.hpp"

void imputeGeno(arma::vec& GVec, 
                double freq, 
                std::vector<uint32_t> posMissingGeno)
{
  int n = posMissingGeno.size();
  
  // should be updated later
  double imputeG = 2 * freq;
  if(freq < 0.05){
    imputeG = 0;
  }
  if(freq > 0.95){
    imputeG = 2;
  }
  
  for(int i = 0; i < n; i++){
    uint32_t posMissing = posMissingGeno.at(i);
    GVec.at(posMissing) = imputeG;
  }
}

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat)
{
  double innerProd = arma::accu(x1Mat % x2Mat);
  return(innerProd);
}

// duplicate each element for (J-1) times: n x 1 -> n(J-1) x 1 
arma::vec Vec2LongVec(arma::vec t_xVec, int n, int J)
{
  arma::vec yVec(n * (J-1));
  int index = 0;
  for(int i = 0; i < n; i ++){
    for(int j = 0; j < J-1; j ++){
      yVec(index) = t_xVec(i);
      index++;
    }
  }
  return yVec;
}

// sum up each (J-1) elements: n(J-1) x 1 -> n x 1
arma::vec LongVec2Vec(arma::vec t_xVec, int n, int J)
{
  arma::vec yVec(n, arma::fill::zeros);
  int index = 0;
  for(int i = 0; i < n; i ++){
    for(int j = 0; j < J-1; j ++){
      yVec(i) += t_xVec(index);
      index++;
    }
  }
  return yVec;
}

// convert: n(J-1) x 1 -> n x (J-1) 
arma::mat Vec2Mat(arma::vec xVec, int n, int J)
{
  arma::mat xMat(n, (J-1));
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      xMat(i, j) = xVec(index);
      index++;
    }
  }
  return(xMat);
}

// convert: n x (J-1) -> n(J-1) x 1
arma::vec Mat2Vec(arma::mat xMat, int n, int J) 
{
  arma::vec xVec(n * (J-1));
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      xVec(index) = xMat(i,j);
      index++;
    }
  }
  return(xVec);
}

// duplicate each row for (J-1) times: n x p -> n(J-1) x p
arma::mat getCovaMat(arma::mat Cova, int n, int J, int p)      
{
  arma::mat CovaMat(n * (J-1), p);
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      CovaMat.row(index) = Cova.row(i);
      index++;
    }
  }
  return CovaMat;
}

// sum each (J-1) cols to 1 col: p x n(J-1) -> p x n (OR) p x (J-1) -> p x 1
arma::mat sumCols(arma::mat t_xMat,
                  int J)
{
  int n = t_xMat.n_cols / (J-1);
  int p = t_xMat.n_rows;
  arma::mat outMat(p, n, arma::fill::zeros);
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      outMat.col(i) += t_xMat.col(index);
      index++;
    }
  }
  return(outMat);
}

// get RPsiP: (J-1) x (J-1) x n 
// Only used in getVarWFast(): 
arma::vec getRPsiR(arma::mat t_muMat,
                   arma::mat t_iRMat,
                   int t_n, int t_J, int t_p)   
{
  // arma::cube RPsiR(J-1, J-1, n, arma::fill::zeros);
  arma::vec RPsiRVec(t_n, arma::fill::zeros);
  arma::mat muRMat = t_muMat.cols(0, t_J-2) / t_iRMat;
  for(int i = 0; i < t_n; i++){
    for(int j1 = 0; j1 < t_J-1; j1++){
      // RPsiR(j1,j1,i) += muRMat(i,j1) / iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      RPsiRVec(i) += muRMat(i,j1) / t_iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      for(int j2 = j1+1; j2 < t_J-1; j2++){
        // RPsiR(j1,j2,i) -= muRMat(i,j1) * muRMat(i,j2);
        RPsiRVec(i) -= 2* muRMat(i,j1) * muRMat(i,j2);
      }
    }
  }
  // return(RPsiR);
  return(RPsiRVec);
}
