
#' Fit a Proportional Odds Logistic Mixed Model (POLMM) for ordinal categorical data analysis
#' 
#' Fit a POLMM with a response varaible of ordinal categorical outcome. To adjust for sample relatedness, we use a random effect with a variance estiamted from Genetic Relationship Matrix (GRM).
#' @param formula an object of class "formula" to be passed to function ordinal::clm(). 
#' @param data an optional data frame containing the variables in the null model. If not specified, the variables are taken from 'environment(formula)'
#' @param PlinkFile a path to Plink files for Genetic Relationship Matrix (GRM).
#' @param subjData a character vector to specify subject IDs in 'formula' and 'data'.
#' @param subjPlink a character vector to specify subject IDs in 'PlinkFile'.
#' @param subjMatch a logical value to indicate if subjData is the same as subjPlink. If TRUE, 'subjData' and 'subjPlink' are not needed any more. Default is FALSE.
#' @param GMatRatio only required if control$LOCO=F. a numeric matrix with genotype of multiple markers (>= 100 markers) to calculate variance ratio (r=VarP/VarW, see 'Details' Section). Each row is for one individual (the same order as subjData) and each column is for one marker. We recommand that these markers should be uncorrelated with the markers in PlinkFile.
#' @param control a list of parameters for controlling the fitting process. More details can be seen in 'Details' Section.
#' @details 
#' Variance VarP should be used to calibrate p values. However, it is computational expensive to calculate VarP for all markers. Hence, we first estimate variance ratio r = VarP / VarW, and then calculate VarW for all markers and estimate VarP = VarW * r. 
#' \itemize{
#' \item{VarP: variance estimated from mixed effect model.}
#' \item{VarW: variance estimated from non-mixed effect model.}
#' }
#' More information about the list of 'control'
#' \itemize{
#' \item{memoryChunk: Size (Gb) for each memory chunk when reading in Plink file [default=2].}
#' \item{seed: An integer as random seed [default=12345678].}
#' \item{tracenrun: Number of runs for trace estimator [default=30].}
#' \item{maxiter: Maximum number of iterations used to fit the null POLMM [default=100].}
#' \item{tolBeta: Positive tolerance for fixed effect coefficients; the iterations converge when |beta - beta_old| / (|beta| + 0.1) < tolBeta [default=0.01].}
#' \item{tolTau: Positive tolerance for random effect variance component (tau); the iterations converge when |tau - tau_old| / (|tau| + 0.1) < tolTau [default=0.02].}
#' \item{tau: Initial value of variance component (tau) [default=0.5].}
#' \item{maxiterPCG: Maximum number of iterations for PCG to converge [default=100].}
#' \item{tolPCG: Positive tolerance for PCG to converge [default=1e-5].}
#' \item{maxiterEps: Maximum number of iterations for cutpoints estimation [default=100].}
#' \item{tolEps: Positive tolerance for cutpoints estimation to converge [default=1e-5].}
#' \item{minMafVarRatio: Minimal value of MAF cutoff to select markers (from Plink file) to estimate variance ratio [default=0.1].}
#' \item{nSNPsVarRatio: Initial number of selected markers to estimate variance ratio [default=30], the number will be automatically added by 10 until the coefficient of variantion (CV) of the variance ratio estimate is below CVcutoff.}
#' \item{CVcutoff: Minimal cutoff of coefficient of variantion (CV) for variance ratio estimation [default=0.0025].}
#' \item{LOCO: Whether to apply the leave-one-chromosome-out (LOCO) approach, if FALSE, 'GMatRatio' is required [default=TRUE].}
#' }
#' @examples check help(POLMM) for an example.
#' @return an R object of 'POLMM_NULL' with the following elements, J is number of categories.
#' \item{N}{Number of individuals in analysis.}
#' \item{M}{Number of markers in Plink file.}
#' \item{Cova}{A stacked matrix (N(J-1) x p) of covariates; p is number of covariates.}
#' \item{beta}{Fixed effect: estimated parameters for covariates.}
#' \item{tau}{Variance component. NOTE: usually underestimated due to Penalized Quasi-Likelihood (PQL).}
#' \item{bVec}{Random effect: estimated to characterize sample relatedness}
#' \item{eps}{Cutpoints for different categories}
#' \item{eta}{Linear predicators for all individuals}
#' \item{yVec}{A matrix (N x J) to indicate the categorical phenotypes.}
#' \item{muMat}{A matrix (N x J), probability in each category.}
#' \item{YMat}{A matrix (N x (J-1)), working variable in each category.}
#' \item{subjIDs}{A vector of individual IDs in null model}
#' \item{LOCOList}{Objects for Step 2 calculation. Might be large if sample size is large and number of chromosomes is large. For example, when analyzing UK Biobank with xx individuals and xx chromosomes, this object is xx Gb.}
#' \item{controlList}{List of control arguments. More details can be seen in 'Details' Section.}
#' \item{sessionInfo}{Information of R session when analysis}
#' \item{time}{Time when null model fitting is finished}
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
    stop("Argument 'PlinkFile' is needed to specify the path to Plink binary files.\n
         Note that bed, bim and fim files should be of the same prefix and at the same directory.")
  if(missing(control)){
    warning("The default setting of 'control' will be used.")
    control = NULL;
  }
  
  ## run ordinal::clm for initial value (author version: 2019.4-25)
  if(missing(data)){
    obj.clm = summary(ordinal::clm(formula))
  }else{
    obj.clm = summary(ordinal::clm(formula, data = data))
  }
  
  if(anyNA(obj.clm$beta))
    stop("Please check collinearity between covariates! Note that intercept term should not be explicitly incorporated.")
  if(length(obj.clm$alpha)==1)
    stop("Number of response levels should be >= 3.")
  
  # incorporate an intercept and then fix the first eps as 0
  beta = c(-1 * obj.clm$alpha[1], obj.clm$beta)   
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
    subjData = paste0("subj",1:n)
    posSampleInPlink = 1:n
  }else{
    if(missing(subjData) | missing(subjPlink)){
      stop("Arguments 'subjData' and 'subjPlink' are needed to match phenotype data and plink file.\n
           If these two are of the same order, you can use 'subjMatch=T' instead.")
    }else{
      posSampleInPlink = match(subjData, subjPlink, 0)
      if(any(posSampleInPlink == 0))
        stop("All individuals in 'subjData' should be also in 'subjPlink'")
      if(length(posSampleInPlink) != n)
        stop("length of 'subjData' should be the same with input data. Note that currently, we do not support missing input.")
    }
  }
  
  # usually the following conditions will not happen
  if(colnames(Cova)[1] != "(Intercept)")
    stop("colnames(Cova)[1] should be (intercept)!!")
  if(any(colnames(Cova)[-1] != names(obj.clm$beta)))
    stop("colnames(Cova)[-1] should be the same as names(obj.clm$beta)")
  
  control = updateCtrl(control);
  if(control$LOCO == FALSE)
  {
    if(missing(GMatRatio))
      stop("Default setting is to use Leave One Chromosome Out (control$LOCO=TRUE). If you set control$LOCO=FALSE, then you have to give a GMatRatio to estimate variance ratio r.")
    if(nrow(GMatRatio)!=n)
      stop("nrow(GMatRatio) should == sample size.")
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
  
  obj_Null$sessionInfo = sessionInfo()
  obj_Null$time = Sys.time()
  obj_Null$subjIDs = subjData
  
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
#' @param GMat a numeric genotype matrix with each row as an individual and each column as a marker. Column names of marker IDs and row names of individual IDs are required.
#' @param obj_Null an output object of the POLMM_Null_Model() function 
#' @param chrVec a character or character vector to specify chromosome(s) of the markers in GMat.
#' @param minMAF a cutoff of the minimal minor allele frequencies (MAFs). Any markers with MAF < cutoff will be excluded from the analysis.
#' @param SPAcutoff a standard deviation cutoff (default=2). If the test statistic lies within the standard deviation cutoff of the mean, p value based on traditional score test is returned. Otherwise, p value is calculated based on saddlepoint approximation.
#' @return an R matrix with the following elements
#' \item{ID}{Marker IDs from colnames(GMat)}
#' \item{MAF}{MAFs of the markers}
#' \item{Stat}{Score statistics}
#' \item{VarW}{Estimated variance (VarW) from non-mixed model}
#' \item{VarP}{Estimated variance after adjusting for variance ratio r (VarP = VarW * r)}
#' \item{pval.norm}{p values calculated from normal approximation}
#' \item{pval.spa}{p values calculated from saddlepoint approximation}
#' @examples 
#' ## We use a Plink file with 10,000 markers and 1,000 subjects to constract GRM for demonstration. 
#' ## For real data analysis, we recommend >= 100,000 common markers (MAF > 0.05 or 0.01).
#' ## Selection of the common markers is similar as in Principle Components Analysis (PCA).
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
#' ## when using function POLMM(), argument chrVec should be from
#' names(objNull$LOCOList)
#' 
#' GMat = matrix(rbinom(10000,2,0.3),1000,10)
#' rownames(GMat) = egData$IID
#' colnames(GMat) = paste0("rs",1:10)
#' outPOLMM = POLMM(GMat, objNull, "1", 0.001)
POLMM = function(GMat,            # n x m matrix
                 obj_Null,
                 chrVec,
                 minMAF,
                 SPAcutoff=2)
{
  if(class(obj_Null)!="POLMM_NULL")
    stop("class(obj_Null) should be 'POLMM_NULL'")
  
  subjIDs = rownames(GMat)
  SNPIDs = colnames(GMat)
  
  if(is.null(subjIDs) | is.null(SNPIDs))
    stop("rownames and colnames of GMat are required!")
  if(any(duplicated(subjIDs)) | any(duplicated(SNPIDs)))
    stop("all elements of either subjIDs or SNPIDs should be unique!")
  
  subjIDs_Null = obj_Null$subjIDs
  subjPos = match(subjIDs_Null, subjIDs, 0)
  if(any(subjPos==0))
    stop("all subjucts in null model should be also in 'GMat'.")
  
  GMat = GMat[subjPos,]
  
  n = nrow(GMat);    # number of individuals
  m = ncol(GMat);
  
  if(missing(chrVec)){
   if(names(obj_Null$LOCOList)!="LOCO=F")
     stop("chrVec is required unless LOCO = F.")
    chrVec = rep("LOCO=F",m)
  }
  
  if(length(chrVec)==1){
    warning(paste0("Chromosome of all markers in analysis are ",chrVec,"."))
    chrVec = rep(chrVec, m)
  }else{
    if(length(chrVec)!=m)
      stop("length of chrVec should be 1 or equal to ncol(GMat)!")
  }
  
  chrVecLOCO = names(obj_Null$LOCOList)
  if(any(!is.element(chrVec, chrVecLOCO)))
    stop("all elements in chrVec should be from names(obj_Null$LOCOList)!")
  
  cat("Totally", m, "markers to analyze!!\n")
  
  uniq_chr = unique(chrVec)
  OutMat = c()
  for(chr in uniq_chr){
    r = obj_Null$LOCOList[[chr]]$VarRatio
    objP = obj_Null$LOCOList[[chr]]$objP
    pos = which(chrVec == uniq_chr)
    for(i in pos){
      GVec = GMat[,i]
      SNPID = colnames(GMat)[i]
      AF = mean(GVec)/2
      MAF = ifelse(AF < 0.5, AF, 1-AF)
      if(MAF < minMAF)
        next
      ## 
      adjG = outputadjGFast(GVec, objP)
      adjGMat = adjG$adjGMat
      Stat = adjG$Stat
      VarW = adjG$VarW
      VarP = VarW * r
      z = abs(Stat)/sqrt(VarP)
      pval.spa = pval.norm = 2*pnorm(-1*z, lower.tail=T)
      if(z > SPAcutoff){   
        pval.spa = Saddle_Prob(Stat, VarP, VarW,
                               adjGMat, objP[["muMat"]], objP[["iRMat"]])   # n x (J-1)
      }
      OutMat = rbind(OutMat,
                     c(SNPID, chr, MAF, Stat, VarW, VarP, pval.norm, pval.spa))
    }
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





