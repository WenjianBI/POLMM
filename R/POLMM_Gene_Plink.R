
#' Test for association between genetic variants and an ordinal categorical variable via Proportional Odds Logistic Mixed Model (POLMM-Gene)
#' 
#' Test for association between genetic variants and an ordinal categorical variable via Proportional Odds Logistic Mixed Model (POLMM-Gene)
#' 
#' @param objNull the output of the POLMM_Null_Model() function 
#' @param GMat a numeric genotype matrix with each row as a subject and each column as a marker in a region or gene.
#'             Name of marker set, column names of marker IDs, and row names of individual IDs are all required.
#'             Missng genotype should be coded as 9 or NA. Both hard-called and imputed genotype are supported.
#' @param chrom a character to specify chromosome of the markers in GMat. Must be specified unless LOCO = F when fitting the null model.
#' @param SparseGRM an object of class "SparseGRM", check help(getSparseGRM) for more details.
#' @param SKAT.control a list of parameters for controlling the POLMM.Gene(). Check 'Details' Section for more details.
#' @return an R matrix with the following elements
#' \item{ID}{Marker IDs from colnames(GMat)}
#' \item{chr}{Chromosome name from chrVec}
#' \item{MAF}{MAFs of the markers}
#' \item{missing.rate}{Missing rates of the markers}
#' \item{Stat}{Score statistics}
#' \item{VarW}{Estimated variance (VarW) from non-mixed model}
#' \item{VarP}{Estimated variance after adjusting for variance ratio r (VarP = VarW * r)}
#' \item{beta}{Estimated effect size: Stat / VarP}
#' \item{pval.norm}{p values calculated from normal approximation}
#' \item{pval.spa}{p values calculated from saddlepoint approximation}
#' \item{switch.allele}{a logical value indicating if the REF/ALT alleles were switched, if AF > 0.5, we use GVec = 2-GVec, and then give switch.allele=T. This is useful to estimate the effect direction.}
#' @details 
#' More information about the list of 'SKAT.control'
#' \itemize{
#' \item{memoryChunk: Size (Gb) for each memory chunk when reading in Plink file [default=2].}
#' \item{seed: An integer as a random seed [default=12345678].}
#' \item{tracenrun: Number of runs for trace estimator [default=30].}
#' \item{maxiter: Maximum number of iterations used to fit the null POLMM [default=100].}
#' \item{tolBeta: Positive tolerance: the iterations converge when |beta - beta_old| / (|beta| + |beta_old| + tolBeta) < tolBeta [default=0.001].}
#' \item{tolTau: Positive tolerance: the iterations converge when |tau - tau_old| / (|tau| + |tau_old| + tolTau) < tolTau [default=0.002].}
#' \item{tau: Initial value of the variance component (tau) [default=0.2].}
#' \item{maxiterPCG: Maximum number of iterations for PCG to converge [default=100].}
#' \item{tolPCG: Positive tolerance for PCG to converge [default=1e-6].}
#' \item{maxiterEps: Maximum number of iterations for cutpoints estimation [default=100].}
#' \item{tolEps: Positive tolerance for cutpoints estimation to converge [default=1e-10].}
#' \item{minMafVarRatio: Minimal value of MAF cutoff to select markers (from Plink file) to estimate variance ratio [default=0.1].}
#' \item{maxMissingVarRatio: Maximal value of missing rate cutoff to select markers (from Plink file) to estimate variance ratio [default=0.1].}
#' \item{nSNPsVarRatio: Initial number of the selected markers to estimate variance ratio [default=20], the number will be automatically added by 10 until the coefficient of variantion (CV) of the variance ratio estimate is below CVcutoff.}
#' \item{CVcutoff: Minimal cutoff of coefficient of variantion (CV) for variance ratio estimation [default=0.0025].}
#' \item{LOCO: Whether to apply the leave-one-chromosome-out (LOCO) approach [default=TRUE].}
#' \item{numThreads: Number of threads (CPUs) to use. Only valid if dense GRM is used, default is "auto", that is, RcppParallel::defaultNumThreads() [default="auto"].}
#' \item{stackSize: Stack size (in bytes) to use for worker threads. For more details, check help(RcppParallel::setThreadOptions) [default="auto"].}
#' \item{grainSize: Grain size of a parallel algorithm sets a minimum chunk size for parallelization. In other words, at what point to stop processing input on separate threads [default=1].}
#' \item{minMafGRM: Minimal value of MAF cutoff to select markers (from Plink files) to construct dense GRM [default=0.01].}
#' \item{maxMissingGRM: Maximal value of missing rate to select markers (from Plink files) to construct dense GRM [default=0.1].}
#' \item{showInfo: Whether to show more detailed information for trouble shooting [default=TRUE].}
#' \item{onlyCheckTime: Not fit the null model, only check the computation time of reading Plink files and running 30 KinbVec() functions [default=FALSE].}
#' }
#' #' @examples 
#' ## We use a Plink file with 10,000 markers and 1,000 subjects to constract GRM for demonstration. 
#' ## For real data analysis, we recommend >= 100,000 common markers (MAF > 0.05 or 0.01).
#' ## Selection of the common markers is similar as in Principle Components Analysis (PCA).
#' famFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.fam", package = "POLMM")
#' PlinkFile = gsub("-ext.fam","-ext",famFile)
#' dataFile = system.file("extdata", "POLMM_data.csv", package = "POLMM")
#' 
#' egData = data.table::fread(dataFile)
#' 
#' ## Fit the null POLMM using the Sparse GRM
#' SparseGRMFile = system.file("SparseGRM", "SparseGRM.RData", package = "POLMM")
#' load(SparseGRMFile)   ## check getSparseGRM() for more details about how to make an R object of "SparseGRM" using Plink files. 
#' objNull = POLMM_Null_Model(as.factor(outcome)~Cova1+Cova2, 
#'                            SparseGRM = SparseGRM,
#'                            data=egData, PlinkFile = PlinkFile, subjData = egData$IID)
#'                            
#' ## If control$seed is not changed, objNull$tau should be 0.8506
#' objNull$tau
#' 
#' obj.SNPSet = list(SNPSet1=c("rs1","rs2","rs3"), SNPSet2=c("rs11","rs12","rs13"))
#' obj.POLMM.Gene = POLMM.Gene.plink(objNull, PlinkFile, obj.SNPSet, chrom="1",SparseGRM)
#' ## when using function POLMM(), argument chrVec should be from
#' names(objNull$LOCOList)
#' 
#' set.seed(123)
#' GMat = matrix(rbinom(10000,2,0.3),1000,10)
#' rownames(GMat) = egData$IID
#' colnames(GMat) = paste0("rs",1:10)
#' chrVec = chrom = "1"  # equivalant to chrVec = rep("1", ncol(GMat))
#' outPOLMM = POLMM(objNull, GMat, chrVec)
#' 
#' outList = POLMM.Gene(objNull, GMat, chrom, SparseGRM)
#' 
#' outPOLMM
#' outList
#' 
#' round(as.numeric(outPOLMM$pval.spa),2)
#' ## [1] 0.89 0.46 0.82 0.71 0.34 0.30 0.20 0.82 0.25 0.71 # using dense GRM
#' ## [1] 0.82 0.46 0.76 0.68 0.36 0.23 0.21 0.80 0.24 0.71 # using sparse GRM
#' 
#' @export
#' @import SKAT
#' @import Matrix
#' @import seqminer

POLMM.Gene.plink = function(objNull,
                            PlinkFile,                    # n x q matrix, where n is number of subjects and m is number of markers
                            obj.SNPSet,
                            chrom,
                            SparseGRM,
                            SKAT.control = NULL)
{
  if(class(objNull) != "POLMM_NULL_Model")
    stop("class(objNull) should be 'POLMM_NULL_Model'.")
  
  # check the setting of SKAT.control, if not specified, the default setting will be used
  SKAT.control = check.SKAT.control(SKAT.control)
  
  # check plink files input
  out.plink = check.PlinkFile(PlinkFile)
  fam.data = out.plink$fam.data
  bim.data = out.plink$bim.data
  
  SubjID.step1 = objNull$subjIDs
  SubjID.step2 = fam.data$V2
  subjIndex_Null = match(SubjID.step1, SubjID.step2, nomatch = 0)
  
  if(any(subjIndex_Null == 0)) 
    stop("All subjects in null model fitting should be also in plink files.")
  
  # check obj.SNPSet
  if(class(obj.SNPSet) != "list")
    stop("class(obj.SNPSet) should be 'list'.")
  
  # update SparseGRM 
  SubjID.step1 = objNull$subjIDs
  SparseGRM = updateSparseGRM(SparseGRM, objNull$subjIDs)
  J = max(objNull$yVec)
  
  if(! chrom %in% names(objNull$LOCOList))
    stop(paste("'chrom' should be from the below chromosomes:",names(objNull$LOCOList)))
  
  # set an objective for Gene-based analysis
  NonZero_cutoff = floor(log(50000, J))  # for efficient resampling (ER)
  StdStat_cutoff = SKAT.control$SPAcutoff
  setPOLMMGENEobj(objNull$controlList$maxiterPCG, 
                  objNull$controlList$tolPCG,
                  objNull$Cova,
                  objNull$yVec,
                  objNull$tau,
                  SparseGRM,
                  objNull$LOCOList,
                  objNull$eta,
                  NonZero_cutoff)
  
  setPOLMMGENEchr(objNull$LOCOList, chrom)
  
  nSet = length(obj.SNPSet)
  
  out_Multi_Set = c()
  names.SNPSet = names(obj.SNPSet)
  for(i in 1:nSet){
    SNPSet = obj.SNPSet[[i]]
    
    markerIndex = match(SNPSet, bim.data$V2, nomatch = 0)
    markerIndex = markerIndex[markerIndex!=0]
    
    GMat = seqminer::readPlinkToMatrixByIndex(PlinkFile, subjIndex_Null, markerIndex)
    colnames(GMat) = bim.data$V2[markerIndex]
    rownames(GMat) = fam.data$V2[subjIndex_Null]
    SetName = names.SNPSet[i]
    
    GMat.list = Check_GMat(GMat, SetName, SubjID.step1,
                           kernel = SKAT.control$kernel, 
                           weights.beta = SKAT.control$weights.beta, 
                           impute.method = SKAT.control$impute.method,
                           impute.MAF.cohort = SKAT.control$impute.MAF.cohort,
                           missing_cutoff = SKAT.control$missing_cutoff, 
                           max_maf = SKAT.control$max_maf)
    
    out_One_Set = POLMM.Gene.Main(GMat.list,          # output of Check_GMat()
                                  NonZero_cutoff,
                                  StdStat_cutoff)
    
    out_Multi_Set = rbind(out_Multi_Set, out_One_Set)
    
  }

  colnames(out_Multi_Set) = c("SetID", "nSNP", "P.SKAT-O", "P.SKAT", "P.Burden",
                              "error.code", "SNP.Info", "SNP.MAF","SNP.AlleleFlip","SNP.beta","SNP.pvalue")
  
  return(out_Multi_Set)
}




