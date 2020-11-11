
#' Test for association between genetic variants and an ordinal categorical variable via Proportional Odds Logistic Mixed Model (POLMM-Gene)
#' 
#' Test for association between genetic variants and an ordinal categorical variable via Proportional Odds Logistic Mixed Model (POLMM-Gene)
#' 
#' @param objNull the output of the POLMM_Null_Model() function 
#' @param PlinkFile character, represents the prefix of PLINK input file.
#' @param obj.SNPSet an R list with names of variant set and elements of variants IDs in the set
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
#' @examples 
#' famFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.fam", package = "POLMM")
#' GenoFile = gsub("-ext.fam", "-ext.bed", famFile)
#' AnnoFile = system.file("extdata", "AnnoFile.txt", package = "POLMM")
#' OutputFile = gsub("AnnoFile","OutputFile",AnnoFile)
#' SparseGRMFile = system.file("SparseGRM", "SparseGRM.RData", package = "POLMM")
#' load(SparseGRMFile)
#' objNullFile = system.file("objNull.RData", package = "POLMM")
#' load(objNullFile)
#' chrom = 1
#' 
#' OUTPUT = POLMM.Region(objNull, AnnoFile, GenoFile, GenoFileIndex = NULL, OutputFile, 
#'                       SparseGRM, chrom, POLMM.control = list(max_maf_region = 0.5))
#'      
#' @export
#' @import SKAT
#' @import Matrix
#' @import seqminer

POLMM.Region = function(objNull,
                        AnnoFile,                # column 1: SNP Set ID, column 2: SNP ID, column 3: Annotation
                        GenoFile,
                        GenoFileIndex = NULL,
                        OutputFile,
                        SparseGRM,
                        chrom,
                        POLMM.control = NULL)
{
  if(class(objNull) != "POLMM_NULL_Model")
    stop("class(objNull) should be 'POLMM_NULL_Model'.")
  
  SubjID.step1 = as.character(objNull$subjIDs);
  
  # check the setting of SKAT.control, if not specified, the default setting will be used
  POLMM.control = check.POLMM.control(POLMM.control)
  print(POLMM.control)
  genoType = setGenoInput(GenoFile, GenoFileIndex, SubjID.step1)
  AnnoList = getAnnoList(AnnoFile)
  
  if(objNull$controlList$LOCO){
    if(!chrom %in% names(objNull$LOCOList))
      stop("'chrom' should be in names(objNull$LOCOList).")
    obj.CHR = objNull$LOCOList[[chrom]]
    
    if(!chrom %in% names(SparseGRM))
      stop("'chrom' should be in names(SparseGRM).")
    SparseGRM.CHR = SparseGRM[[chrom]]
  }
  
  SPmatR.CHR = makeSPmatR(SparseGRM.CHR, objNull$subjIDs)
  
  setPOLMMobjInR(obj.CHR$muMat,
                 obj.CHR$iRMat,
                 objNull$Cova,
                 objNull$yVec,
                 SPmatR.CHR,
                 objNull$tau,
                 POLMM.control$printPCGInfo,
                 POLMM.control$tolPCG,
                 POLMM.control$maxiterPCG)
  
  memory.chunk = POLMM.control$memory.chunk
  
  # to be continued
  n = length(SubjID.step1)
  J = max(objNull$yVec)
  p = ncol(objNull$Cova)
  NonZero_cutoff = floor(log(50000, J))  # for efficient resampling (ER)
  maxMarkers = getMaxMarkers(POLMM.control$memory.chunk, n, J, p);
  
  print(paste0("The current POLMM.control$memory.chunk is ", POLMM.control$memory.chunk,"(GB)."))
  print(paste0("Based on the sample size, we divide region with more than ", maxMarkers, " markers into multiple chunks to save the memory usage."))
  print("If the memory usage still exceed the memory you request, please set a smaller POLMM.control$memory.chunk.")
  
  StdStat_cutoff = POLMM.control$SPAcutoff;
  
  OUT.Region = c()
  for(i in 1:length(AnnoList)){
    Region = names(AnnoList)[i]
    print(paste0("Analyzing Region of ", Region,"......"))
    markers = AnnoList[[i]]
    OutList = MAIN_REGION(markers,
                          NonZero_cutoff,
                          POLMM.control$SPAcutoff,
                          maxMarkers,
                          OutputFile,
                          POLMM.control$missing_cutoff,
                          POLMM.control$max_maf_region)
    
    VarSVec = diag(OutList$VarSMat)
    StdStatVec = OutList$StatVec / sqrt(VarSVec)
    QVec = StdStatVec^2
    adjPVec = PVec = pchisq(QVec, lower.tail = FALSE, df = 1)
    
    weights = Get_Weights(POLMM.control$kernel, OutList$freqVec, POLMM.control$weights.beta)
    
    # r0 = pmax(r0, 1)
    r0 = 1
    wr0 = sqrt(r0) * weights
    
    wStatVec = OutList$StatVec * weights
    wadjVarSMat = t(OutList$VarSMat * wr0) * wr0
    
    out_SKAT_List = try(SKAT:::Met_SKAT_Get_Pvalue(Score = wStatVec, 
                                                   Phi = wadjVarSMat,
                                                   r.corr = POLMM.control$r.corr, 
                                                   method = "optimal.adj", 
                                                   Score.Resampling = NULL),
                        silent = TRUE)
    
    # betaVec = StatVec / adjVarSVec; 
    if(class(out_SKAT_List) == "try-error"){
      Pvalue = c(NA, NA, NA)
      error.code = 2
    }else if(!any(c(0,1) %in% out_SKAT_List$param$rho)){
      Pvalue = c(NA, NA, NA)
      error.code = 3
    }else{
      pos0 = which(out_SKAT_List$param$rho == 0)
      pos1 = which(out_SKAT_List$param$rho == 1)
      Pvalue = c(out_SKAT_List$p.value,                  # SKAT-O
                 out_SKAT_List$param$p.val.each[pos0],   # SKAT
                 out_SKAT_List$param$p.val.each[pos1])   # Burden Test
      error.code = 0
    }
    
     
    
    OUT.Region = rbind(OUT.Region,
                       c(Region, 
                         length(OutList$markerVec), 
                         Pvalue, 
                         paste(OutList$markerVec, collapse = ","),
                         paste(OutList$freqVec, collapse = ","),
                         paste(OutList$flipVec, collapse = ","),
                         paste(OutList$StatVec, collapse = ","),
                         paste(VarSVec, collapse = ","),
                         paste(adjPVec, collapse = ",")))
  }
  
  tempFiles = list.files(path = dirname(OutputFile),
                         pattern = paste0("^",basename(OutputFile),".*\\.bin$"),
                         full.names = T)
  file.remove(tempFiles)
  # out_Multi_Set = rnorm(1)
  # colnames(out_Multi_Set) = c("SetID", "nSNP", "P.SKAT-O", "P.SKAT", "P.Burden",
  #                             "error.code", "SNP.Info", "SNP.MAF","SNP.AlleleFlip","SNP.beta","SNP.pvalue")
  
  # return(out_Multi_Set)
  return(OUT.Region)
}

check.POLMM.control = function(POLMM.control)
{
  # the below is the default setting
  default.POLMM.control = list(impute.method = "fixed",  # the below is shared parameters for both single marker testing and region-based testing
                               missing_cutoff = 0.15, 
                               max_maf_region = 0.01,
                               SPAcutoff = 2,
                               memory.chunk = 4,
                               kernel = "linear.weighted",   # the below is parameters for region-based testing
                               method = "SKAT-O",
                               weights.beta = c(1,25),
                               weights = NULL,
                               r.corr = NULL,
                               printPCGInfo = FALSE,
                               tolPCG = 1e-5,
                               maxiterPCG = 100)
  
  # use the default setting or update it
  if(!is.null(POLMM.control)){
    if(class(POLMM.control) != "list")
      stop("Argument of POLMM.control should be 'list'.")
    
    # update the setting
    ctrl.nm = names(POLMM.control)
    for(nm in ctrl.nm){
      default.POLMM.control[[nm]] = POLMM.control[[nm]]
    }
  }
  
  POLMM.control = default.POLMM.control
  
  # check the parameters
  if(! POLMM.control$kernel %in% c("linear", "linear.weighted"))
    stop("'kernel' should be 'linear' or 'linear.weighted'. Check 'Details' for more details.")
  
  if(length(POLMM.control$weights.beta)!=2)
    stop("length of 'weights.beta' should be 2. Check 'Details' for more details.")
  
  if(any(POLMM.control$weights.beta < 0))
    stop("the two elements in 'weights.beta' should be non-negative. Check 'Details' for more details.")
  
  if(! POLMM.control$impute.method %in% c("fixed", "bestguess", "random"))
    stop("'impute.method' should be 'fixed','bestguess', or 'random'. Check 'Details' for more details.")
  
  if(POLMM.control$missing_cutoff > 1 | POLMM.control$missing_cutoff < 0)
    stop("'missing_cutoff' should be between 0 and 1. We recommand using 0.15.")
  
  if(POLMM.control$max_maf_region > 1 | POLMM.control$max_maf_region <= 0)
    stop("'max_maf_region' should be between 0 and 1. We recommand using 0.05 or 0.01.")
  
  method = POLMM.control$method
  
  if(! method %in% c("SKAT", "Burden", "SKAT-O"))
    stop("'method' should be 'SKAT', 'Burden', or 'SKAT-O'. Check 'Details' for more details.")
  
  if(is.null(POLMM.control$r.corr)){
    if(method == "SKAT")
      POLMM.control$r.corr = 0;
    
    if(method == "Burden")
      POLMM.control$r.corr = 1;
    
    if(method == "SKAT-O")
      POLMM.control$r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1);  # r.corr = 0 is SKAT, r.corr = 1 is Burden Test
  }
  
  if(any(POLMM.control$r.corr < 0 | POLMM.control$r.corr > 1))
    stop("'r.corr' should be between 0 and 1. Check 'Details' for more details.")
  
  return(POLMM.control)
}

setGenoInput = function(GenoFile, GenoFileIndex, SampleInModel)
{
  if(missing(GenoFile))
    stop("Argument 'GenoFile' is required.")
  
  if(!file.exists(GenoFile))
    stop(paste("Cannot find genotype file of", GenoFile))
  
  GenoFileExt = tools::file_ext(GenoFile);
  
  if(GenoFileExt != "bed")
    stop("Current version of POLMM.Region only supports plink file input.")  # support more genotype input later
  
  # if the genotype input is plink files
  if(GenoFileExt == "bed"){
    genoType = "PLINK"
    if(is.null(GenoFileIndex)){  # default setting: the same prefix for bim and fam files
      GenoFileIndex = c(gsub("bed$","bim",GenoFile),
                        gsub("bed$","fam",GenoFile))
    }
    if(length(GenoFileIndex) != 2)
      stop("If plink file is used, argument 'GenoFileIndex' should be 'NULL' or a character vector of c(bimFile, famFile).")
    
    bimFile = GenoFileIndex[1]
    famFile = GenoFileIndex[2]
    if(!file.exists(bimFile)) stop(paste("Cannot find bim file of", bimFile))
    if(!file.exists(famFile)) stop(paste("Cannot find fam file of", famFile))
    
    setPLINKobjInR(bimFile, famFile, GenoFile, SampleInModel)
  }
  
  # add more types of genotype input later
  
  # return genotype
  return(genoType)
}

getAnnoList = function(AnnoFile)
{
  if(!file.exists(AnnoFile))
    stop(paste("Cannot find AnnoFile in", AnnoFile))
  
  AnnoData = data.table::fread(AnnoFile, header = T, stringsAsFactors = F);
  if(any(colnames(AnnoData) != c("GENE","SNP","ANNO")))
    stop("header of AnnoFile should be c('GENE', 'SNP', 'ANNO').")
  
  AnnoList = list()
  uGENE = unique(AnnoData$GENE)
  for(g in uGENE){
    uSNP = unique(AnnoData$SNP[AnnoData$GENE == g])
    AnnoList[[g]] = uSNP
  }
  return(AnnoList)
}

# make a matrix that can be passed to arma::sp_mat
makeSPmatR = function(SparseGRM,   # three columns of ID1, ID2, and value
                      subjData)
{
  SparseGRM = subset(SparseGRM, ID1 %in% subjData & ID2 %in% subjData)
  
  row.diag = which(SparseGRM$ID1 == SparseGRM$ID2)
  SparseGRM.diag = SparseGRM[row.diag,]
  SparseGRM.off.d1 = SparseGRM[-1*row.diag,]
  SparseGRM.off.d2 = data.frame(ID1 = SparseGRM.off.d1$ID2,
                                ID2 = SparseGRM.off.d1$ID1,
                                value = SparseGRM.off.d1$value)
  
  SparseGRM = rbind(SparseGRM.diag, SparseGRM.off.d1, SparseGRM.off.d2)
  
  ID1 = SparseGRM$ID1;
  ID2 = SparseGRM$ID2;
  value = SparseGRM$value;
  
  if(any(!is.element(subjData, ID1)))
    stop("At least one of subjects in `subjData` is not in `SparseGRM`.")
  
  location1 = match(ID1, subjData);
  location2 = match(ID2, subjData);
  
  locations = rbind(location1 - 1,  # -1 is to convert R to C++
                    location2 - 1)
  
  SPmatR = list(locations = locations,
                values = value)
  
  return(SPmatR)
}

getMaxMarkers = function(memory.chunk,
                         n, J, p)
{
  # THE BELOW include some large matrix that takes most of the memory usage
  
  # VarSMat = mat(indexPassingQC, indexPassingQC)
  # adjGMat = fmat(n, maxMarkers): 4*n*maxMarkers
  # ZPZ_adjGMat = fmat(n, maxMarkers): 4*n*maxMarkers
  # CovaMat = mat(n(J-1) x p): 8*n*(J-1)*p
  # arma::mat m_XXR_Psi_RX;  // XXR_Psi_RX ( n x p )
  # arma::mat m_XR_Psi_R;    // XR_Psi_R ( p x n ), sum up XR_Psi_R ( p x n(J-1) ) for each subject 
  # arma::vec m_RymuVec;     // n x 1: row sum of the n x (J-1) matrix R %*% (yMat - muMat)
  # arma::mat m_iSigmaX_XSigmaX; // n(J-1) x p
  # arma::mat m_CovaMat;     // n(J-1) x p
  
  fixed.memory = 8*(3*(n*(J-1)*p)+2*(n*p)) / 1e9
  if(memory.chunk < fixed.memory)
    stop(paste0("Please give POLMM.control$memory.chunk greater than ", fixed.memory,"."))
  
  maxMarkers = (memory.chunk - fixed.memory) * 1e9 / (8*n)
  maxMarkers = maxMarkers / 2;
  maxMarkers = floor(maxMarkers)
  
  return(maxMarkers)
}

####### ---------- Get Weights from MAF ---------- #######

Get_Weights = function(kernel, freqVec, weights.beta)
{
  if(kernel == "linear"){
    weights = rep(1, length(freqVec)) 
  }
  
  if(kernel == "linear.weighted"){
    weights = dbeta(freqVec, weights.beta[1], weights.beta[2])
  }
  
  return(weights)
}
