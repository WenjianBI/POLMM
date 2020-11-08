// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// setPLINKobjInR
void setPLINKobjInR(std::string t_bimFile, std::string t_famFile, std::string t_bedFile, std::vector<std::string> t_SampleInModel);
RcppExport SEXP _POLMM_setPLINKobjInR(SEXP t_bimFileSEXP, SEXP t_famFileSEXP, SEXP t_bedFileSEXP, SEXP t_SampleInModelSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_bimFile(t_bimFileSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_famFile(t_famFileSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_bedFile(t_bedFileSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type t_SampleInModel(t_SampleInModelSEXP);
    setPLINKobjInR(t_bimFile, t_famFile, t_bedFile, t_SampleInModel);
    return R_NilValue;
END_RCPP
}
// setPOLMMobjInR
void setPOLMMobjInR(arma::mat t_muMat, arma::mat t_iRMat, arma::mat t_Cova, arma::vec t_yVec, Rcpp::List t_SPmatR, double t_tau, bool t_printPCGInfo, double t_tolPCG, int t_maxiterPCG);
RcppExport SEXP _POLMM_setPOLMMobjInR(SEXP t_muMatSEXP, SEXP t_iRMatSEXP, SEXP t_CovaSEXP, SEXP t_yVecSEXP, SEXP t_SPmatRSEXP, SEXP t_tauSEXP, SEXP t_printPCGInfoSEXP, SEXP t_tolPCGSEXP, SEXP t_maxiterPCGSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type t_muMat(t_muMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_iRMat(t_iRMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_Cova(t_CovaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_yVec(t_yVecSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type t_SPmatR(t_SPmatRSEXP);
    Rcpp::traits::input_parameter< double >::type t_tau(t_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type t_printPCGInfo(t_printPCGInfoSEXP);
    Rcpp::traits::input_parameter< double >::type t_tolPCG(t_tolPCGSEXP);
    Rcpp::traits::input_parameter< int >::type t_maxiterPCG(t_maxiterPCGSEXP);
    setPOLMMobjInR(t_muMat, t_iRMat, t_Cova, t_yVec, t_SPmatR, t_tau, t_printPCGInfo, t_tolPCG, t_maxiterPCG);
    return R_NilValue;
END_RCPP
}
// MAIN_REGION
Rcpp::List MAIN_REGION(std::vector<std::string> t_MarkerReqstd, double t_NonZero_cutoff, double t_StdStat_cutoff, int t_maxMarkers, std::string t_outputFile, double t_missingRate_cutoff, double t_maxMAF_cutoff);
RcppExport SEXP _POLMM_MAIN_REGION(SEXP t_MarkerReqstdSEXP, SEXP t_NonZero_cutoffSEXP, SEXP t_StdStat_cutoffSEXP, SEXP t_maxMarkersSEXP, SEXP t_outputFileSEXP, SEXP t_missingRate_cutoffSEXP, SEXP t_maxMAF_cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type t_MarkerReqstd(t_MarkerReqstdSEXP);
    Rcpp::traits::input_parameter< double >::type t_NonZero_cutoff(t_NonZero_cutoffSEXP);
    Rcpp::traits::input_parameter< double >::type t_StdStat_cutoff(t_StdStat_cutoffSEXP);
    Rcpp::traits::input_parameter< int >::type t_maxMarkers(t_maxMarkersSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_outputFile(t_outputFileSEXP);
    Rcpp::traits::input_parameter< double >::type t_missingRate_cutoff(t_missingRate_cutoffSEXP);
    Rcpp::traits::input_parameter< double >::type t_maxMAF_cutoff(t_maxMAF_cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(MAIN_REGION(t_MarkerReqstd, t_NonZero_cutoff, t_StdStat_cutoff, t_maxMarkers, t_outputFile, t_missingRate_cutoff, t_maxMAF_cutoff));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _POLMM_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _POLMM_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _POLMM_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _POLMM_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_POLMM_setPLINKobjInR", (DL_FUNC) &_POLMM_setPLINKobjInR, 4},
    {"_POLMM_setPOLMMobjInR", (DL_FUNC) &_POLMM_setPOLMMobjInR, 9},
    {"_POLMM_MAIN_REGION", (DL_FUNC) &_POLMM_MAIN_REGION, 7},
    {"_POLMM_rcpparma_hello_world", (DL_FUNC) &_POLMM_rcpparma_hello_world, 0},
    {"_POLMM_rcpparma_outerproduct", (DL_FUNC) &_POLMM_rcpparma_outerproduct, 1},
    {"_POLMM_rcpparma_innerproduct", (DL_FUNC) &_POLMM_rcpparma_innerproduct, 1},
    {"_POLMM_rcpparma_bothproducts", (DL_FUNC) &_POLMM_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_POLMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
