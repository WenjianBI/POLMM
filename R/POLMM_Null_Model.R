
#' Fit a Proportional Odds Logistic Mixed Model (POLMM) for ordinal categorical data analysis
#' 
#' Fit a POLMM with a response varaible of ordinal categorical outcome. If the outcome is binary, POLMM degenerates to a logistic mixed model and the package performs similarly as SAIGE. 
#' Two options are provided to adjust for sample relatedness. If argument SparseGRM is specified, we use a sparse Genetic Relationship Matrix (GRM, check ?getSparseGRM for more details), otherwise we use a dense GRM.
#' @param formula a formula of the form response ~ predictors to be passed to function ordinal::clm(). The response should be an ordered factor and the predictors should not include intercept term. 
#' @param data an optional data frame containing the variables in the null model. If not specified, the variables are taken from 'environment(formula)'
#' @param subjData a character vector to specify subject IDs in 'formula' and 'data'.
#' @param PlinkFile a path to Plink files for dense Genetic Relationship Matrix (GRM). To construct dense GRM or to extract markers to estimate variant ratio.
#' @param bimFile a path to Plink bim file. Default is paste0(PlinkFile,".bim") if not specified.
#' @param bedFile a path to Plink bed file. Default is paste0(PlinkFile,".bed") if not specified.
#' @param famFile a path to Plink fam file. Default is paste0(PlinkFile,".fam") if not specified.
#' @param SparseGRM an object of class "SparseGRM", check ?getSparseGRM for more details.
#' @param GMatRatio an object of class "GMatRatio", usually not used.
#' @param control a list of parameters for controlling the fitting process. Check 'Details' Section for more details.
#' @details 
#' More information about the list of 'control'
#' \itemize{
#' \item{memoryChunk: Size (Gb) for each memory chunk when reading in Plink file [default=2].}
#' \item{seed: An integer as random seed [default=12345678].}
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
#' \item{LOCO: Whether to apply the leave-one-chromosome-out (LOCO) approach, if FALSE, argument 'GMatRatio' is required [default=TRUE].}
#' \item{numThreads: Number of threads (CPUs) to use if dense GRM is used, default is "auto", that is, RcppParallel::defaultNumThreads() [default="auto"].}
#' \item{stackSize: Stack size (in bytes) to use for worker threads. For more details, check ?RcppParallel::setThreadOptions [default="auto"].}
#' \item{grainSize: Grain size of a parallel algorithm sets a minimum chunk size for parallelization. In other words, at what point to stop processing input on separate threads [default=1].}
#' \item{minMafGRM: Minimal value of MAF cutoff to select markers (from Plink files) to construct dense GRM [default=0.01].}
#' \item{maxMissingGRM: Maximal value of missing rate to select markers (from Plink files) to construct dense GRM [default=0.1].}
#' \item{showInfo: Whether to show more detailed information for trouble shooting [default=TRUE].}
#' \item{onlyCheckTime: Not fit the null model, only check the computation time of reading Plink files and running 30 KinbVec() functions [default=FALSE].}
#' }
#' @examples 
#' ## We use a Plink file with 10,000 markers and 1,000 subjects to constract GRM for demonstration. 
#' ## For real data analysis, we recommend >= 100,000 common markers (MAF > 0.05 or 0.01).
#' ## Selection of the common markers is similar as in Principle Components Analysis (PCA).
#' famFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.fam", package = "POLMM")
#' PlinkFile = gsub("-ext.fam","-ext",famFile)
#' dataFile = system.file("extdata", "POLMM_data.csv", package = "POLMM")
#' 
#' egData = data.table::fread(dataFile)
#'
#' ## Fit the null POLMM using the Dense GRM
#' objNull = POLMM_Null_Model(as.factor(outcome)~Cova1+Cova2, 
#'                            data=egData, PlinkFile = PlinkFile, subjData = egData$IID)
#' 
#' ## If control$seed is not changed, objNull$tau should be 0.7353
#' objNull$tau
#' 
#' ## Fit the null POLMM using the Sparse GRM
#' SparseGRMFile = system.file("extdata", "SparseGRM.RData", package = "POLMM")
#' load(SparseGRMFile)   ## check getSparseGRM() for more details about how to make an R object of "SparseGRM" using Plink files. 
#' objNull = POLMM_Null_Model(as.factor(outcome)~Cova1+Cova2, 
#'                            SparseGRM = SparseGRM,
#'                            data=egData, PlinkFile = PlinkFile, subjData = egData$IID)
#'                            
#' ## If control$seed is not changed, objNull$tau should be 0.8506
#' objNull$tau
#' 
#' 
#' ## If you want more accurate parameter estimation, smaller tolBeta/tolTau can be used. 
#' objNull = POLMM_Null_Model(as.factor(outcome)~Cova1+Cova2, 
#'                            data=egData, PlinkFile = PlinkFile, subjData = egData$IID,
#'                            control=list(tolTau=0.001))  # Default toTau = 0.002
#' objNull$tau  # 0.7368
#'
#' ## This is for me to check the computational efficiency
#' objNull = POLMM_Null_Model(as.factor(outcome)~Cova1+Cova2, 
#'                            data=egData, PlinkFile = PlinkFile, subjData = egData$IID,
#'                            control=list(onlyCheckTime = T, numThreads = 4))  # Default nThreads = 0, that is, RcppParallel::defaultNumThreads()
#'   
#' @return an R object of 'POLMM_NULL' with the following elements, J is number of categories.
#' \item{N}{Number of individuals in analysis.}
#' \item{M}{Number of markers in Plink file.}
#' \item{Cova}{A stacked matrix (N(J-1) x p) of covariates; p is number of covariates.}
#' \item{beta}{Fixed effect: estimated parameters for covariates.}
#' \item{tau}{Variance component. NOTE: usually underestimated due to Penalized Quasi-Likelihood (PQL).}
#' \item{bVec}{Random effect: estimated to characterize sample relatedness}
#' \item{eps}{Cutpoints for different categories}
#' \item{eta}{Linear predicators for all individuals}
#' \item{yVec}{A vector of the categorical phenotypes.}
#' \item{muMat}{A matrix (N x J), probability in each category.}
#' \item{YMat}{A matrix (N x (J-1)), working variable in each category.}
#' \item{subjIDs}{A vector of individual IDs in null model}
#' \item{LOCOList}{Objects for Step 2 calculation. Might be large if sample size is large and number of chromosomes is large. For example, when analyzing UK Biobank with xx individuals and xx chromosomes, this object is xx Gb.}
#' \item{controlList}{List of control arguments. More details can be seen in 'Details' Section.}
#' \item{sessionInfo}{Information of R session when analysis}
#' \item{time}{Time when null model fitting is finished}
#' @param formula an object of class "formula" to be passed to function ordinal::clm() or glm(). 
#' @param data an optional data frame containing the variables in the null model. If not specified, the variables are taken from 'environment(formula)'
#' @param subjData a character vector to specify subject IDs in 'formula' and 'data'.
#' @param PlinkFile a path to Plink files for dense Genetic Relationship Matrix (GRM). Required if argument SparseGRM is not specified.
#' @param bimFile a path to Plink bim file. Default is paste0(PlinkFile,".bim") if not specified.
#' @param bedFile a path to Plink bed file. Default is paste0(PlinkFile,".bed") if not specified.
#' @param famFile a path to Plink fam file. Default is paste0(PlinkFile,".fam") if not specified.
#' @param SparseGRM an object of class "SparseGRM", check ?getSparseGRM for more details.
#' @param GMatRatio an object of class "GMatRatio", check ?getGMatRatio for more details. Required if control$LOCO=F.
#' @param control a list of parameters for controlling the fitting process. Check 'Details' Section for more details.

POLMM_Null_Model = function(formula,
                            data,
                            subjData,
                            PlinkFile,
                            bimFile,
                            bedFile,
                            famFile,
                            SparseGRM,
                            GMatRatio,
                            control)
{
  ######## -------- run ordinal::clm for initial value ---------- ##########
  ## ordinal Package version: 2019.12-10
  ## Cova: an R matrix with intercept term; it can be used for predicators of class "factor"
  if(missing(data)){
    obj.clm = summary(ordinal::clm(formula))
    Cova = model.matrix(formula)
  }else{
    obj.clm = summary(ordinal::clm(formula, data = data))
    Cova = model.matrix(formula, data = data)  
  }
  
  posNA = c(obj.clm$na.action)       # position of NA input in formula/data
  n = obj.clm$n                      # sample size in formula/data after excluding NA
  J = length(obj.clm$y.levels)       # number of response levels after excluding NA
  
  if(anyNA(obj.clm$beta))
    stop("Please check collinearity between covariates! Note that intercept term should not be explicitly incorporated.")
  if(J > 10)
    stop("Number of response levels should be <= 10. Otherwise, we suggest using methods designed for continuous trait such as fastGWA or BOLT-LMM.")
  
  ## incorporate an intercept term and fix the first cutpoint (eps) as 0
  beta = c(-1 * obj.clm$alpha[1], obj.clm$beta)
  eps = c(0, obj.clm$alpha[-1] - obj.clm$alpha[1])
  
  yVec = obj.clm$y
  
  # the following situations usually do not happen
  if(colnames(Cova)[1] != "(Intercept)")
    stop("colnames(Cova)[1] should be (intercept)!!")
  if(any(colnames(Cova)[-1] != names(obj.clm$beta)))
    stop("colnames(Cova)[-1] should be the same as names(obj.clm$beta)")

  ######### -------- check input argument ---------- ##########
  
  if(missing(subjData))
    stop("Argument 'subjData' is required to match subjects in 'formula/data' and subjects in 'PlinkFile/SparseGRM' (dense/sparse GRM).")
  
  if(!is.null(posNA)){
    if(max(posNA)>length(subjData))
      stop("length(subjData) should be the same as the sample size in 'formula/data'.")
    subjData = subjData[-1*posNA]     # remove subjID with NA input in formula/data
  }
    
  if(length(subjData)!=n)
    stop("length(subjData) should be the same as the sample size in 'formula/data'.")
  
  if(any(duplicated(subjData)))
    stop("NO duplicated name in subjData is allowed.")
  
  if(missing(SparseGRM)){
    # dense GRM
    flagSparseGRM = F
    print("Argument 'SparseGRM' is NOT specified. Will use dense GRM!")
    SparseGRM = list()
  }else{                   
    # sparse GRM
    flagSparseGRM = T
    print("Argument 'SparseGRM' is specified. Will use sparse GRM!")
    if(class(SparseGRM)!="SparseGRM") stop("class(SparseGRM)!='SparseGRM'")
    SparseGRM = updateSparseGRM(SparseGRM, subjData)
  }
  
  if(missing(control)){
    print("Argument 'control' is not specified. Check 'Details' Section for more details about the default setting of 'control'.")
    control = NULL;
  }
  
  control = updateCtrl(control);
  message("The control setting is as belows.")
  print(control)
  
  ### Check ?RcppParallel::setThreadOptions for detals about parallel setting
  RcppParallel:::setThreadOptions(numThreads = control$numThreads, stackSize = control$stackSize)
  
  if(missing(GMatRatio)){
    flagGMatRatio = F
    GMat = matrix(0)
  }else{
    flagGMatRatio = T
    if(class(GMatRatio)!="GMatRatio")
      stop("class(GMatRatio) should be 'GMatRatio'.")
    subjGMat = GMatRatio$subjGMat
    GMat = GMatRatio$GMat
    posGMat = match(subjData, subjGMat, 0)
    if(any(posGMat == 0))
      stop("All subjects in 'subjData' should be in 'GMatRatio'.")
    
    GMat = GMat[posGMat,]
  }
  
  # (no GMatRatio) or (no SparseGRM)
  if(!flagGMatRatio || !flagSparseGRM){
    
    bimFile = ifelse(!missing(bimFile), bimFile,
                     ifelse(!missing(PlinkFile), paste0(PlinkFile,".bim"), stop("Either bimFile or PlinkFile should be specified.")))
    bedFile = ifelse(!missing(bedFile), bedFile,
                     ifelse(!missing(PlinkFile), paste0(PlinkFile,".bed"), stop("Either bedFile or PlinkFile should be specified.")))
    famFile = ifelse(!missing(famFile), famFile,
                     ifelse(!missing(PlinkFile), paste0(PlinkFile,".fam"), stop("Either famFile or PlinkFile should be specified.")))
    
    if(!file.exists(bimFile)) stop("Could not find bimFile or paste0(PlinkFile,'.bim')")
    if(!file.exists(bedFile)) stop("Could not find bedFile or paste0(PlinkFile,'.bed')")
    if(!file.exists(famFile)) stop("Could not find famFile or paste0(PlinkFile,'.fam')")
    
    famData = read.table(famFile)
    subjGRM = famData[,2];  # column 2 is the subject ID (IID)
    posSampleInPlink = match(subjData, subjGRM, 0)
    if(any(posSampleInPlink == 0))
      stop("All subjects in 'subjData' should be in famFile or SparseGRM")
  }
  
  ######## -------------- fit the null POLMM --------------  ###########
  
  bVec = rep(0, n)  # initiate random effect of 0
  
  obj_Null = fitPOLMMcpp(t_flagSparseGRM = flagSparseGRM,       # if 1, then use SparseGRM, otherwise, use DenseGRM
                         t_flagGMatRatio = flagGMatRatio,       # if 1, then use GMatRatio, otherwise, extract from Plink files
                         t_bimfile = bimFile,
                         t_famfile = famFile,
                         t_bedfile = bedFile,
                         t_posSampleInPlink = posSampleInPlink,
                         t_Cova = Cova,
                         t_yVec = yVec,                         # should be from 1 to J
                         t_beta = beta,
                         t_bVec = bVec,
                         t_eps = eps,
                         t_tau = control$tau,
                         t_GMatRatio = GMat,                    # only used if m_LOCO = FALSE
                         t_SparseGRM = SparseGRM,
                         t_controlList = control)
  
  obj_Null$sessionInfo = sessionInfo()
  obj_Null$time = Sys.time()
  obj_Null$subjIDs = subjData
  
  class(obj_Null) = "POLMM_NULL_Model"
  return(obj_Null)
}

############### -------- update argument of control ------- #############

updateCtrl = function(control.new){
  # default setting
  control = list(memoryChunk = 2,
                 seed = 12345678,
                 tracenrun = 30,
                 maxiter = 100,
                 tolBeta = 0.001,
                 tolTau = 0.002,
                 tau = 0.2,
                 maxiterPCG = 100,
                 tolPCG = 1e-6,
                 maxiterEps = 100,
                 tolEps = 1e-10,
                 minMafVarRatio = 0.1,
                 maxMissingVarRatio = 0.1, 
                 nSNPsVarRatio = 20,
                 CVcutoff = 0.0025,
                 LOCO = T,
                 numThreads = "auto",
                 stackSize = "auto",
                 grainSize = 1,
                 minMafGRM = 0.01,
                 maxMissingGRM = 0.1,
                 showInfo = T,
                 onlyCheckTime = F)
  
  if(is.null(control.new))
    return(control)
  
  if(!is.list(control.new))
    stop("Argument 'control' should be a list!")
  
  for(nm in names(control.new))
    control[[nm]] = control.new[[nm]]
  
  return(control)
}

updateSparseGRM = function(SparseGRM, subjData){
  names = names(SparseGRM)
  KinMatListR = list();
  
  for(excludeChr in names){
    tempGRM1 = SparseGRM[[excludeChr]]
    # updated on 05-14-2020
    tempGRM2 = data.frame(ID1=tempGRM1$ID2,
                          ID2=tempGRM1$ID1,
                          value=tempGRM1$value)
    
    tempGRM = rbind(tempGRM1, tempGRM2)
    tempGRM = tempGRM[-1*which(duplicated(tempGRM)),]
    
    ID1 = tempGRM$ID1;
    
    if(any(!is.element(subjData, ID1)))
      stop("At least one of subjects is not in SparseGRM.")
    
    ID2 = tempGRM$ID2;
    value = tempGRM$value;
    location1 = match(ID1, subjData);
    location2 = match(ID2, subjData);
    pos = which(!is.na(location1) & !is.na(location2))
    locations = rbind(location1[pos]-1,  # -1 is to convert R to C++
                      location2[pos]-1)
    
    value = value[pos];
    nSubj = length(subjData);
    KinMatListR[[excludeChr]] = list(locations = locations,
                                     values = value,
                                     nSubj = nSubj)
  }
  
  return(KinMatListR)
}

