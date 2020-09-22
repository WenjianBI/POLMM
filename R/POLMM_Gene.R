
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

POLMM.Gene = function(objNull,
                      GMat,                    # n x q matrix, where n is number of subjects and m is number of markers
                      chrom,
                      SparseGRM,
                      SKAT.control = NULL)
{
  if(class(objNull) != "POLMM_NULL_Model")
    stop("class(objNull) should be 'POLMM_NULL_Model'.")
  
  # check the setting of SKAT.control, if not specified, the default setting will be used
  SKAT.control = check.SKAT.control(SKAT.control)
  
  
  GMat.list = Check_GMat(GMat, objNull, chrom,
                         kernel = SKAT.control$kernel, 
                         weights.beta = SKAT.control$weights.beta, 
                         impute.method = SKAT.control$impute.method,
                         impute.MAF.cohort = SKAT.control$impute.MAF.cohort,
                         missing_cutoff = SKAT.control$missing_cutoff, 
                         max_maf = SKAT.control$max_maf)
  
  if(GMat.list$error){
    outList = list(p.value = NA, param = NA, p.value.resampling = NA)
    return(outList)
  }
  
  # extract information from GMat.list
  GMat = GMat.list$GMat
  weights = GMat.list$weights
  AlleleFlip.Vec = GMat.list$AlleleFlip.Vec
  MAF.Vec = GMat.list$MAF.Vec  # this might be slightly different from "true" MAF due to the genotype imputation, but should be OK since we have limited the missing_rate
  
  # update SparseGRM 
  SparseGRM = updateSparseGRM(SparseGRM, objNull$subjIDs)
  J = max(objNull$yVec)
  NonZero_cutoff = floor(log(50000, J))  # for efficient resampling (ER)
  StdStat_cutoff = SKAT.control$SPAcutoff
  
  # set an objective for Gene-based analysis
  setPOLMMGENEobj(objNull$controlList$maxiterPCG, 
                  objNull$controlList$tolPCG,
                  objNull$Cova,
                  objNull$yVec,
                  objNull$tau,
                  SparseGRM,
                  objNull$LOCOList,
                  objNull$eta,
                  NonZero_cutoff)
  
  setPOLMMGENEchr(objNull$LOCOList, 
                  chrom)
  
  ##
  OutList = getStatVarS(GMat, NonZero_cutoff, StdStat_cutoff)

  # output basic information including Stat, stdStat, and VarSMat
  StatVec = as.vector(OutList$StatVec)  # as.vector(): transform matrix to vector
  StdStatVec = as.vector(OutList$StdStatVec)
  VarSMat = OutList$VarSMat
  VarSVec = diag(VarSMat)
  
  # adjust the results based on SPA or ER
  idxERVec = as.vector(OutList$idxERVec)
  idxSPAVec = as.vector(OutList$idxSPAVec)
  muMat = OutList$muMat
  muMat1 = muMat[,-1*J]
  iRMat = OutList$iRMat
  VarWVec = as.vector(OutList$VarWVec)
  Ratio0Vec = as.vector(OutList$Ratio0Vec)
  adjGMat = OutList$adjGMat

  out_SPA_ER = adj_SPA_ER(GMat, StatVec, VarSVec, StdStatVec,
                          idxERVec, idxSPAVec,
                          muMat1, iRMat, VarWVec, Ratio0Vec, adjGMat) # first check the non-robust version, no adjustment
  
  adjPVec = out_SPA_ER$adjPVec
  PVec = out_SPA_ER$PVec
  adjVarSVec = out_SPA_ER$adjVarSVec
  
  #########
  
  wStatVec = StatVec * weights
  
  if(any(is.infinite(adjVarSVec))) 
    stop("any(is.infinite(adjVarSVec))") # I want to know when this will happen
  
  adjVarSVec = ifelse(is.infinite(adjVarSVec), 0, adjVarSVec)
  StatVec = ifelse(is.infinite(adjVarSVec), 0, StatVec)
  
  r0 = adjVarSVec / VarSVec  # adjVarSVec might be 0
  wr0 = sqrt(r0) * weights
  
  wadjVarSMat = t(VarSMat * wr0) * wr0
  
  ## use another ratio to adjust for variance matrix based on Burden Test
  GMat.BT = matrix(colSums(t(GMat) * weights), ncol=1)
  OutList = getStatVarS(GMat.BT, 0, StdStat_cutoff)
  
  # Burden test p value is more significant than the pre-given cutoff
  if(ncol(OutList$adjGMat) == 1){
    Stat = OutList$StatVec[1,1]
    VarS = OutList$VarSMat[1,1]
    VarW = OutList$VarWMat[1,1]
    Ratio0 = OutList$Ratio0Vec[1,1]
    K1roots = c(0,0)
    posG1 = which(GMat.BT[,1] != 0)
    adjGVec = OutList$adjGMat[posG1,1]
    # calculate p value of Burden test from saddlepoint approximation
    res.spa = fastSaddle_Prob(Stat, VarS, 
                              VarW, Ratio0, 
                              K1roots,
                              adjGVec, muMat1[posG1,], iRMat[posG1,])
    pval.BT = res.spa$pval
    print(paste("pval.BT:",pval.BT))
  }
  
  # p.value_burden <- Saddle_Prob(q.sum, mu = mu, g = g.sum,
  #                               Cutoff = 2, alpha = 2.5 * 10^-6)$p.value
  # v1 = rep(1, dim(G2_adj_n)[1])
  # VarQ = t(v1) %*% G2_adj_n %*% v1
  # p.m <- dim(G)[2]
  # Q_b = p.m^2 * rowMeans(zscore.all_1)^2
  # VarQ_2 = Q_b/qchisq(p.value_burden, df = 1, ncp = 0, lower.tail = FALSE,
  #                     log.p = FALSE)
  #
  # if (VarQ_2 == 0) {
  #   r = 1
  # }     else {
  #   r = VarQ/VarQ_2
  # }
  # r = min(r, 1)
  #
  # if (dim(G2_adj_n)[2]==1){
  #   Phi_temp=as.matrix(G2_adj_n *1/r)
  # } else {
  #   Phi_temp=as.matrix(G2_adj_n %*% diag(rep(1/r, dim(G2_adj_n)[2])))
  # }
  # 
  # q = ncol(mVarSMat)
  # if(q == 1){
  #   Phi_temp = as.matrix(mVarSMat * 1/r)
  # }else{
  #   Phi_temp = as.matrix(Phi_temp %*% diag(rep(1/r, q)))
  # }
  
  r.all = SKAT.control$r.corr
  r.all = ifelse(r.all >= 0.999, 0.999, r.all)
  
  outList = try(SKAT:::Met_SKAT_Get_Pvalue(Score = wStatVec, Phi = wadjVarSMat,
                                           r.corr = r.all, method = "optimal.adj", Score.Resampling = NULL),
                silent = TRUE)
  
  outList = list(outList = outList,
                 PVec = PVec,
                 adjPVec = adjPVec)
  
  # OutMat = as.data.frame(OutMat, stringsAsFactors=F)
  return(outList)
}

adj_SPA_ER = function(GMat, StatVec, VarSVec, StdStatVec,
                      idxERVec, idxSPAVec,
                      muMat1, iRMat, VarWVec, Ratio0Vec, adjGMat)
{
  # select markers to run robust testing based on ER and SPA
  pos.ER = which(idxERVec == 1)
  pos.SPA = which(idxSPAVec == 1)
  idx.SPA = 1;
  
  QVec = StdStatVec^2
  adjPVec = PVec = pchisq(QVec, lower.tail = FALSE, df = 1)
  adjVarSVec = VarSVec
  K1roots = c(0,0)
  
  for(i in pos.ER){
    adjPVec[i] = getPvalERtoR(GMat[,i])
    adjVarSVec[i] = StatVec[i]^2 / qchisq(adjPVec[i], df = 1, lower.tail = FALSE)
  }
  
  for(i in pos.SPA){
    posG1 = which(GMat[,i] != 0)
    adjGVec = adjGMat[,idx.SPA]
    res.spa = fastSaddle_Prob(StatVec[i], VarSVec[i], 
                              VarWVec[idx.SPA], Ratio0Vec[idx.SPA], 
                              K1roots,
                              adjGVec[posG1], muMat1[posG1,], iRMat[posG1,])
    adjPVec[i] = res.spa$pval
    adjVarSVec[i] = StatVec[i]^2 / qchisq(adjPVec[i], df = 1, lower.tail = FALSE)
    idx.SPA = idx.SPA + 1
  }
  
  adjVarSVec = ifelse(adjPVec == 0, StatVec^2 / 500, adjVarSVec)
  
  out_SPA_ER = list(adjVarSVec = adjVarSVec, adjPVec = adjPVec, PVec = PVec)
  
  return(out_SPA_ER)
}

####### ---------- Check Argument of SKAT ------------ #########

check.SKAT.control = function(SKAT.control)
{
  # the below is the default setting
  default.SKAT.control = list(kernel = "linear.weighted",
                              method = "SKAT-O",
                              weights.beta = c(1,25),
                              weights = NULL,
                              impute.method = "bestguess",
                              impute.MAF.cohort = "step1",
                              r.corr = NULL,
                              missing_cutoff = 0.15, 
                              max_maf = 0.05,
                              SPAcutoff = 2)
                              # to be continued)
  
  # use the default setting or update it
  if(is.null(SKAT.control)){
    SKAT.control = default.SKAT.control
  }else{
    ctrl.nm = names(SKAT.control)
    for(nm in ctrl.nm){
      default.SKAT.control[[nm]] = SKAT.control[[nm]]
    }
    SKAT.control = default.SKAT.control
  }
  
  # check the parameters
  if(! SKAT.control$kernel %in% c("linear", "linear.weighted"))
    stop("'kernel' should be 'linear' or 'linear.weighted'. Check 'Details' for more details.")
  
  if(length(SKAT.control$weights.beta)!=2)
    stop("length of 'weights.beta' should be 2. Check 'Details' for more details.")
  
  if(any(SKAT.control$weights.beta < 0))
    stop("the two elements in 'weights.beta' should be non-negative. Check 'Details' for more details.")
  
  if(! SKAT.control$impute.method %in% c("fixed", "bestguess", "random"))
    stop("'impute.method' should be 'fixed','bestguess', or 'random'. Check 'Details' for more details.")
  
  if(! SKAT.control$impute.MAF.cohort %in% c("step1", "step2"))
    stop("'impute.MAF.cohort' should be 'step1' or 'step2'. Check 'Details' for more details.")
  
  if(SKAT.control$missing_cutoff > 1 | SKAT.control$missing_cutoff < 0)
    stop("'missing_cutoff' should be between 0 and 1. We recommand using 0.15.")
  
  if(SKAT.control$max_maf > 1 | SKAT.control$max_maf <= 0)
    stop("'max_maf' should be between 0 and 1. We recommand using 0.05 or 0.01.")
  
  method = SKAT.control$method
  if(! method %in% c("SKAT", "Burden", "SKAT-O"))
    stop("'method' should be 'SKAT', 'Burden', or 'SKAT-O'. Check 'Details' for more details.")
  
  if(is.null(SKAT.control$r.corr)){
    if(method == "SKAT")
      SKAT.control$r.corr = 0;
    
    if(method == "Burden")
      SKAT.control$r.corr = 1;
    
    if(method == "SKAT-O")
      SKAT.control$r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1);
  }
  
  if(any(SKAT.control$r.corr < 0 | SKAT.control$r.corr > 1))
    stop("'r.corr' should be between 0 and 1. Check 'Details' for more details.")
  
  return(SKAT.control)
}


####### ---------- Check the GMat matrix, and do imputation ---------- #######

Check_GMat = function(GMat, objNull, chrom, kernel, weights.beta, impute.method, impute.MAF.cohort, missing_cutoff, max_maf)
{
  # Number of subjects and SNPs
  SetID = names(GMat);
  SNPsID = colnames(GMat)
  SubjID.step2 = rownames(GMat)
  SubjID.step1 = objNull$subjIDs
  
  if(! chrom %in% names(objNull$LOCOList))
    stop(paste("'chrom' should be from the below chromosomes:",names(objNull$LOCOList)))
  
  if(is.null(SetID))
    stop("names(GMat), that is, name of SNP set, is requried.")
  
  if(is.null(SNPsID))
    stop("colnames(GMat), that is, names of SNPs in the set, is requried.")
  
  if(is.null(SubjID.step2))
    stop("rownames(GMat), that is, names of subjects in the set, is requried.")
  
  pos.Step1.In.Step2 = match(SubjID.step1, SubjID.step2, nomatch = 0)
  
  if(any(pos.Step1.In.Step2 == 0))
    stop("All subjects in step 1 (fitting the null model) should be in step 2 (association testing).")
  
  if(impute.MAF.cohort == "step2"){  # calculate MAF (for imputation) based on all subjects with genotypes
    GMat[GMat > 1.8] = 2  # change GMat close to 0/2 to 0/2
    GMat[GMat < 0.2] = 0  # mainly for the future usage of Efficient Resampling (ER)
    GMat[GMat == 9] = NA  # plink format use 9 as missing genotype
    MAF.Vec = colMeans(GMat, na.rm = T) / 2   # MAF for all markers
    GMat = GMat[pos.Step1.In.Step2,]
    MAC.Vec = colSums(GMat, na.rm = T)
  }
  
  if(impute.MAF.cohort == "step1"){  # calculate MAF (for imputation) based on subjects in null model fitting
    GMat = GMat[pos.Step1.In.Step2,] 
    GMat[GMat > 1.8] = 2  # change GMat close to 0/2 to 0/2
    GMat[GMat < 0.2] = 0  # mainly for the future usage of Efficient Resampling (ER)
    GMat[GMat == 9] = NA  # plink format use 9 as missing genotype
    MAF.Vec = colMeans(GMat, na.rm = T) / 2   # MAF for all markers
    MAC.Vec = MAF.Vec * 2 * ncol(GMat)
  }
  
  # flip allele based on minor allele frequency
  AlleleFlip.Vec = (MAF.Vec > 0.5)
  MAF.Vec[AlleleFlip.Vec] = 1 - MAF.Vec[AlleleFlip.Vec]
  MAC.Vec[AlleleFlip.Vec] = 2 * ncol(GMat) - MAC.Vec[AlleleFlip.Vec]
  
  # remove SNPs with large missing rate or MAF, or MAC == 0
  MissingRate.Vec = colMeans(is.na(GMat))                    # missing rate for all markers
  pos.PassG = which(MissingRate.Vec <= missing_cutoff & MAF.Vec <= max_maf & MAC.Vec != 0)
  
  if(length(pos.PassG) == 0){
    msg = sprintf("In %s, ALL SNPs have been excluded. P value = 1", SetID)
    warning(msg, call.=FALSE)
    
    re = list(error = 1) 
    return(re)
  }
  
  GMat = GMat[,pos.PassG,drop=F]
  MAF.Vec = MAF.Vec[pos.PassG]
  AlleleFlip.Vec = AlleleFlip.Vec[pos.PassG]
  MissingRate.Vec = MissingRate.Vec[pos.PassG]
  
  # genotype imputation
  pos.imputeG = which(MissingRate.Vec > 0)
  
  if(length(pos.imputeG) > 0){
    for(i in pos.imputeG){
      MAF = MAF.Vec[i]
      IDX_NA = which(is.na(GMat[,i]))
      n_NA = length(IDX_NA)
      # different imputation method
      if(impute.method == "fixed") 
        GMat[IDX_NA, i] = 2 * MAF
      if(impute.method == "bestguess") 
        GMat[IDX_NA, i] = 0  # ifelse(MAF < 0.5, 0, 2). Since we only consider low-frequency variants and have flipped the allele
      if(impute.method == "random")
        GMat[IDX_NA, i] = rbinom(n_NA, 2, MAF)
    }
  }
  
  # get weights for rare variants analysis
  weights = Get_Weights(kernel, MAF.Vec, weights.beta)
  
  return(list(GMat = GMat, 
              weights = weights,
              AlleleFlip.Vec = AlleleFlip.Vec,
              MAF.Vec = MAF.Vec, 
              error = 0))
}

####### ---------- Get Weights from MAF ---------- #######

Get_Weights = function(kernel, MAF.Vec, weights.beta)
{
  if(kernel == "linear"){
    weights = rep(1, length(MAF.Vec)) 
  }
  
  if(kernel == "linear.weighted"){
    weights = dbeta(MAF.Vec, weights.beta[1], weights.beta[2])
  }
  
  return(weights)
}
