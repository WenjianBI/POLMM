
#' Fit a Proportional Odds Logistic Mixed Model (POLMM) for ordinal categorical data analysis
#' 
#' Fit a POLMM to analyze ordinal categorical outcome. We use a random effect with a variance of Genetic Relationship Matrix (GRM) to adjust for sample relatedness.
#' @param formula an object of class "formula": a symbolic description of the null model to be fitted
#' @param data an optional data frame containing the variables in the model. If not specified, the variables are taken from 'environment(formula)'
#' @param PlinkFile Path to plink file for creating the genetic relationship matrix (GRM).
#' @param subjData an R vector to specify subject IDs in 'formula' and 'data'
#' @param subjPlink an R vector to specify subject IDs in 'PlinkFile'
#' @param GMatRatio an R matrix with raw genotype of multiple SNPs (e.g. 30 SNPs) to calculate Variance Ratio. Only required when control$LOCO=F. Note that it should be uncorrelated with SNPs in PlinkFile.
#' @param control an R list with for controlling the fitting process. More details can be seen in 'Details' Section.
#' @details More information about list of 'control'
#' \itemize{
#' \item{memoryChunk: Size (Gb) for each memory chunk when reading in Plink file [default=2].}
#' \item{seed: A single integer value as random number seeds [default=12345678].}
#' \item{tracenrun: A single integer value, initial number of runs for trace estimator [default=30].}
#' \item{maxiter: A single integer value, maximum number of iterations used to fit the null POLMM [default=100].}
#' \item{tolBeta: Tolerance for fixed effect coefficients to converge when fitting the null POLMM [default=0.01].}
#' \item{tolTau: Tolerance for variance component (tau) of random effect to converge when fitting the null POLMM [default=0.02].}
#' \item{tau: Initial values for variance component (tau) [default=0.5].}
#' \item{maxiterPCG: Maximum number of iterations for PCG [default=100].}
#' \item{tolPCG: Tolerance for PCG to converge [default=1e-5].}
#' \item{maxiterEps: Maximum number of iterations for cutpoints estimation [default=100].}
#' \item{tolEps: Tolerance for cutpoints estimation to converge [default=1e-5].}
#' \item{minMafVarRatio: MAF cutoff to select SNPs for variance ratio estimation [default=0.1].}
#' \item{nSNPsVarRatio: Initial number of markers to be randomly selected for variance ratio estimation [default=30], the number will be automatically added by 10 untial the coefficient of variantion (CV) for the variance ratio estimate is below CVcutoff.}
#' \item{CVcutoff: coefficient of variantion (CV) cutoff for variance ratio estimation [default=0.0025].}
#' \item{LOCO: Whether to apply the leave-one-chromosome-out (LOCO) approach, if FALSE, 'GMatRatio' is required [default=TRUE].}
#' }
#' @example check help(SAIGE.POLMM) for one example.
#' @return an R object of 'POLMM_NULL' with the following elements
#' \item{N}{Number of individuals in analysis.}
#' \item{M}{Number of markers in Plink file for GRM construction.}
#' \item{Cova}{Stacked matrix (N(J-1) x p) of covariates, p is number of covariates.}
#' \item{beta}{fixed effect parameters for covariates.}
#' \item{tau}{variance component. NOTE: usually underestimated due to Penalized Quasi-likelihood.}
#' \item{bVec}{random effect for sample relatedness}
#' \item{eps}{cutpoints for different categories}
#' \item{eta}{linear predicator}
#' \item{yVec}{Matrix (N x J) to indicate the ordinal categorical phenotypes, J is number of categories.}
#' \item{muMat}{Matrix (N x J), probability in each category for all individuals.}
#' \item{YMat}{Matrix (N x (J-1)), working variable in each category for all individuals.}
#' \item{LOCOList}{Objects for Step 2 calculation. Might be large if sample size is large and number of chromosomes is large. For example, when analyzing UK Biobank with xx individuals and xx chromosomes, this object is xx Gb.}
#' \item{controlList}{List of control arguments. More details about control can be seen in 'Details' Section.}
#' \item{version}{Current version of POLMM package}
#' \item{time}{Time when finish fitting null model fitting}
POLMM_Null_Model = function(formula,
                            data,
                            PlinkFile,
                            subjData,
                            subjPlink,
                            subjMatch = FALSE,
                            GMatRatio,
                            control)
{
  if(missing(PlinkFile))
    stop("Argument 'PlinkFile' is needed to specify the plink file for GRM.\n
         Note that bed, bim and fim files should be at the same directory with the same prefix.")
  if(missing(control)){
    warning("No argument 'control' is given, The default setting of 'control' will be used!")
    control = NULL;
  }
  if(missing(data)){
    obj.clm = summary(ordinal::clm(formula))
  }else{
    obj.clm = summary(ordinal::clm(formula, data = data))
  }
  
  ## run ordinal::clm for initial value (ordinal package version: 2019.4-25)
  if(anyNA(obj.clm$beta))
    stop("Please check collinearity between covariates! Note that intercept term should not be explicitly incorporated.")
  if(length(obj.clm$alpha)==1)
    stop("Number of outcome levels should be >= 3.")
  
  beta = c(-1 * obj.clm$alpha[1], obj.clm$beta)   # incorporate an intercept and then fix the first eps as 0
  eps = c(0, obj.clm$alpha[-1] - obj.clm$alpha[1])
  
  ## read in data
  mc = match.call(expand.dots = FALSE)
  Data = getData(mc, contrasts)
  yVec = as.numeric(Data$Y)
  Cova = Data$X
  n = length(yVec)
  bVec = rep(0, n)
  
  if(subjMatch){
    warning("Arguments 'subjData' and 'subjPlink' are ignored since 'subjMatch=T'.")
    posSampleInPlink = 1:n
  }else{
    if(missing(subjData) | missing(subjPlink)){
      stop("Arguments 'subjData' and 'subjPlink' are needed to match phenotype data and plink file.\n
           If these two are of the same order, you can use 'subjMatch=T' instead.")
    }else{
      posSampleInPlink = match(subjData, subjPlink, 0)
      if(any(posSampleInPlink==0))
        stop("All elements in subjData should be in subjPlink")
      if(length(posSampleInPlink)!=n)
        stop("length of subjData should be the same with input data. Note that currently, we do not support NA input.")
    }
  }
  
  # usually the following conditions will not happen
  if(colnames(Cova)[1]!="(Intercept)")
    stop("colnames(Cova)[1] should be (intercept)!!")
  if(any(colnames(Cova)[-1] != names(obj.clm$beta)))
    stop("colnames(Cova)[-1] should be the same as names(obj.clm$beta)")
  
  control = updateCtrl(control);
  if(control$LOCO == FALSE)
  {
    if(missing(GMatRatio))
      stop("Default setting is to use LOCO (Leave One Chromosome Out), if you choose to set LOCO=F, then you have to give a matrix of GMatRatio to estimate the variance ratio r.")
    if(nrow(GMatRatio)!=n)
      stop("nrow(GMatRatio) should equal to the sample size.")
    if(ncol(GMatRatio)<100)
      stop("ncol(GMatRatio) should be >= 100.")
  }else{
    GMatRatio = matrix(1,1,1)
  }

  ### fit a POLMM

  obj_Null = fitNullcpp(Plink = PlinkFile,
                        posSampleInPlink = posSampleInPlink,
                        CovaR = Cova,
                        yVecR = yVec,
                        betaR = beta,
                        bVecR = bVec,
                        epsR = eps,
                        tauR = control$tau,
                        GMatRatioR = GMatRatio,
                        controlListR = control)
  
  obj_Null$version = "v0.1.0"
  obj_Null$time = Sys.time()
  
  class(obj_Null) = "POLMM_NULL"
  return(obj_Null)
}

getData = function(mc, contrasts)
{
  if (missing(contrasts)) 
    contrasts <- NULL
  
  m <- match(c("formula", "data"), names(mc), 0L)  # position of c("formula","data") at names(mc)
  nm <- names(as.list(mc))
  if (!"formula" %in% nm) 
    stop("Argument 'formula' is required!!!")
  mf <- mc[c(1L, m)]                               # only retain formula and data
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())  
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  X <- model.matrix(mt, mf, contrasts)
  return(list(Y = Y, X = X))
}

updateCtrl = function(control.new){
  # default setting
  control = list(memoryChunk = 2,
                 seed = 12345678,
                 tracenrun = 30,
                 maxiter = 100,
                 tolBeta = 0.01,
                 tolTau = 0.02,
                 tau = 0.5,
                 maxiterPCG = 100,
                 tolPCG = 1e-5,
                 maxiterEps = 100,
                 tolEps = 1e-5,
                 minMafVarRatio = 0.1,
                 nSNPsVarRatio = 30,
                 CVcutoff = 0.0025,
                 LOCO = T)
  
  if(is.null(control.new))
    return(control)
  
  if(!is.list(control.new))
    stop("Argument 'control' should be a list!!")
  
  nmctrl = names(control.new)
  for(nm in nmctrl){
    control[[nm]] = control.new[[nm]]
  }
  return(control)
}

#' Test for association between genotype and an ordinal categorical variable via Proportional Odds Logistic Mixed Model (POLMM)
#' 
#' Test for association between genotype and an ordinal categorical variable via Proportional Odds Logistic Mixed Model (POLMM)
#' @param GMat a numeric genotype matrix with each row as an individual and each column as a marker. Column names is required
#' @param obj_Null output object of function FitNull
#' @param chrVec a character or vector to specify the markers chromosome
#' @param minMAF lower cutoff of genetic variants minor allele frequencies (MAFs) for analysis.
#' @return an R matrix with the following elements
#' \item{SNPID}{SNP IDs from colnames(GMat)}
#' \item{MAF}{Minor allele frequencies}
#' \item{Stat}{Score statistics}
#' \item{VarW}{Estimated variance (VarW) from non mixed effect model}
#' \item{VarP}{Estimated variance after adjusting for variance ratio r (VarP = VarW * r)}
#' \item{pval.norm}{p values calculated from normal approximation}
#' \item{pval.spa}{p values calculated from saddlepoint approximation}
#' @examples 
#' ## The example files use plink files with 10000 markers and 1000 subjects. 
#' ## For real data analysis, we recommend >= 100000 markers.
#' famFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.fam", package = "POLMM")
#' PlinkFile = gsub("-ext.fam","-ext",famFile)
#' dataFile = system.file("extdata", "POLMM_data.csv", package = "POLMM")
#' 
#' famData = data.table::fread(famFile)
#' egData = data.table::fread(dataFile)
#'
#' ## Fit the null POLMM 
#' objNull = POLMM_Null_Model(as.factor(outcome)~Cova1+Cova2, data=egData, PlinkFile = PlinkFile, subjData = egData$IID, subjPlink = famData$V2)
#' 
#' ## If control$seed is not changed, objNull$tau should be 0.705
#' objNull$tau
#' 
#' ## when using SAIGE.POLMM, chrVec should be from
#' names(objNull$LOCOList)
#' 
#' GMat = matrix(rbinom(10000,2,0.3),1000,10)
#' colnames(GMat) = paste0("rs",1:10)
#' outPOLMM = POLMM(GMat, objNull, "1", 0.001)
POLMM = function(GMat,            # n x m matrix
                 obj_Null,
                 chrVec,
                 minMAF)
{
  if(class(obj_Null)!="POLMM_NULL")
    stop("class(obj_Null) should be 'POLMM_NULL'")
  
  n = nrow(GMat);
  m = ncol(GMat);
  
  if(missing(chrVec)){
   if(names(obj_Null$LOCOList)!="LOCO=F")
     stop("chrVec is required unless LOCO = F.")
    chrVec = rep("LOCO=F",m)
  }
  
  
  if(length(chrVec)==1){
    chrVec = rep(chrVec, m)
  }else{
    if(length(chrVec)!=m)
      stop("length of chrVec should be 1 or equals to ncol(GMat)!!")
  }
  
  chrVecLOCO = names(obj_Null$LOCOList)
  if(any(!is.element(chrVec, chrVecLOCO)))
    stop("all elements in chrVec should be in names(obj_Null$LOCOList)!!")
  
  cat("Totally", m, "markers to analyze!!\n")
  
  OutMat = c()
  for(i in 1:m){
    GVec = GMat[,i]
    SNPID = colnames(GMat)[i]
    AF = mean(GVec)/2
    MAF = ifelse(AF < 0.5, AF, 1-AF)
    if(MAF < minMAF)
      next
    chr = chrVec[i]
    r = obj_Null$LOCOList[[chr]]$VarRatio
    objP = obj_Null$LOCOList[[chr]]$objP
    ## 
    adjG = outputadjGFast(GVec, objP)
    adjGMat = adjG$adjGMat
    Stat = adjG$Stat
    VarW = adjG$VarW
    VarP = VarW * r
    pval.norm = 2*pnorm(-1*abs(Stat)/sqrt(VarP), lower.tail=T)
    if(pval.norm < 0.0455){  # spa-2 
      pval.spa = Saddle_Prob(Stat, VarP, VarW,
                             adjGMat, objP[["muMat"]], objP[["iRMat"]])   # n x (J-1)
    }else{
      pval.spa = pval.norm
    }
    OutMat = rbind(OutMat,
                   c(SNPID, chr, MAF, Stat, VarW, VarP, pval.norm, pval.spa))
  }
  colnames(OutMat) = c("SNPID", "chr", "MAF", "Stat", "VarW", "VarP", "pval.norm", "pval.spa")
  return(OutMat)
}

Korg = function(t, 
                muMat,        # n x (J-1) 
                cMat,         # n x (J-1)
                m1)           # sum(muMat * cMat)
{
  n = length(t)
  out = rep(0, n)
  
  for(i in 1:n){
    t1 = t[i]
    temp1Mat = - muMat + muMat * exp(cMat * t1)
    temp1 = log(1 + rowSums(temp1Mat))
    out[i] = sum(temp1) - m1 * t1
  }
  return(out)
}

K1 = function(t, 
              muMat, 
              cMat, 
              m1)
{
  n = length(t)	
  out = rep(0, n)
  
  for(i in 1:n){
    t1 = t[i]
    temp1Mat = - muMat + muMat * exp(cMat * t1)
    temp2Mat = exp(cMat * t1) * muMat * cMat
    temp1 = 1 + rowSums(temp1Mat)
    temp2 = rowSums(temp2Mat)
    out[i] = sum(temp2/temp1) - m1
  }
  return(out)
}

K2 = function(t, 
              muMat, 
              cMat)
{
  n = length(t)
  out = rep(0, n)
  
  for(i in 1:n){
    t1 = t[i]
    temp1Mat = - muMat + muMat * exp(cMat * t1)
    temp2Mat = exp(cMat * t1) * muMat * cMat
    temp3Mat = temp2Mat * cMat
    
    temp1 = 1 + rowSums(temp1Mat)
    temp2 = rowSums(temp2Mat)
    temp3 = rowSums(temp3Mat)
    
    out[i] = sum((temp3*temp1-temp2^2)/temp1^2, na.rm=TRUE)
  }
  return(out)
}


getroot_K1 = function(Stat,
                      init.t,
                      muMat, 
                      cMat, 
                      m1,
                      tol = .Machine$double.eps^0.25,
                      maxiter = 100)
{
  t = init.t;
  K1_eval = 0
  diff.t = Inf
  
  for(iter in 1:maxiter){
    old.t = t
    old.diff.t = diff.t
    old.K1 = K1_eval
    
    K1_eval = K1(t, muMat, cMat, m1) - Stat
    K2_eval = K2(t, muMat, cMat)
    
    diff.t = -1 * K1_eval / K2_eval
    if(sign(K1_eval)!=sign(old.K1)){
      while(abs(diff.t) > abs(old.diff.t) - tol){
        diff.t = diff.t/2
      }
    }
    # print(t)
    if(abs(diff.t) < tol) break;
    t = old.t + diff.t
  }
  
  if(iter == maxiter) converge = F
  else converge = T
  
  return(list(root = t,
              iter = iter,
              converge = converge))
}


Get_Saddle_Prob = function(Stat,
                           zeta, 
                           muMat, 
                           cMat,
                           m1,          # sum(muMat * cMat)
                           lower.tail)
{
  k1 = Korg(zeta, muMat, cMat, m1)
  k2 = K2(zeta, muMat, cMat)

  if(is.finite(k1) && is.finite(k2))
  {
    w = sign(zeta) * sqrt(2 * (zeta * Stat - k1))
    v = zeta * sqrt(k2)
    
    Z.test = w + 1/w * log(v/w)
    pval = pnorm(Z.test, lower.tail = lower.tail)
    
  }else{
    pval = 0
  }
  
  return(pval)
}

Saddle_Prob = function(Stat,
                       VarP,
                       VarW,
                       adjGMat, # n x (J-1)
                       muMat,  # n x J
                       iRMat)   # n x (J-1)
{
  muMat1 = muMat[,-1*ncol(muMat),drop=F]  # remove the last column
  adjStat = Stat / sqrt(VarP)
  cMat = adjGMat / iRMat / sqrt(VarW);
  m1 = sum(muMat1 * cMat)
  
  out.uni1 = getroot_K1(abs(adjStat), 0, muMat1, cMat, m1)
  out.uni2 = getroot_K1(-1*abs(adjStat), 0, muMat1, cMat, m1)
  
  if(out.uni1$converge == TRUE & out.uni2$converge == TRUE){
    p1 = Get_Saddle_Prob(abs(adjStat), out.uni1$root, muMat1, cMat, m1, lower.tail=FALSE)
    p2 = Get_Saddle_Prob(-1*abs(adjStat), out.uni2$root, muMat1, cMat, m1, lower.tail=TRUE)
    pval = p1 + p2;
  }else{
    print("SPA does not converge, use normal approximation p value.")
    pval = 2*pnorm(-1*abs(Stat)/sqrt(VarP), lower.tail=T)
  }
  
  return(pval)
}





