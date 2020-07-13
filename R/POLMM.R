
#' Test for association between genotype and an ordinal categorical variable via Proportional Odds Logistic Mixed Model (POLMM)
#' 
#' Test for association between genotype and an ordinal categorical variable via Proportional Odds Logistic Mixed Model (POLMM)
#' 
#' @param objNull output object of the POLMM_Null_Model() function 
#' @param Geno.mtx a numeric genotype matrix with each row as an individual and each column as a marker. 
#'                 Column names of marker IDs and row names of individual IDs are required.
#'                 Missng genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param chrVec a character or a character vector to specify chromosome(s) of the markers in Geno.mtx. Must be specified unless LOCO = F.
#' @param SPAcutoff a standard deviation cutoff (default=2). If the standardized test statistic < SPAcutoff, normal approximation is used, otherwise, saddlepoint approximation is used.
#' @param minMAF a cutoff of the minimal minor allele frequencies (MAFs). Any markers with MAF < minMAF will be excluded from the analysis.
#' @param maxMissing a cutoff of the maximal missing rate. Any markers with missing rate > maxMissing will be excluded from the analysis.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#'                      "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e. 2p where p is MAF).
#' @param G.model a character string (default: "Add") to specify additive ("Add"), dominant ("Dom"), or recessive ("Rec") model. 
#'                If "Dom", GVec = ifelse(GVec >= 1, 1, 0), if "Rec", GVec = ifelse(GVec <= 1, 0, 1). Be very careful if the gneotyp is imputed data.
#' @return an R matrix with the following elements
#' \item{ID}{Marker IDs from colnames(Geno.mtx)}
#' \item{chr}{Chromosome name from chrVec}
#' \item{MAF}{MAFs of the markers}
#' \item{missing.rate}{Missing rates of the markers}
#' \item{Stat}{Score statistics}
#' \item{VarW}{Estimated variance (VarW) from non-mixed model}
#' \item{VarP}{Estimated variance after adjusting for variance ratio r (VarP = VarW * r)}
#' \item{beta}{Estimated effect size: Stat / VarP}
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
#' ## when using function POLMM(), argument chrVec should be from
#' names(objNull$LOCOList)
#' 
#' set.seed(123)
#' Geno.mtx = matrix(rbinom(10000,2,0.3),1000,10)
#' rownames(Geno.mtx) = egData$IID
#' colnames(Geno.mtx) = paste0("rs",1:10)
#' chrVec = "1"  # equivalant to chrVec = rep("1", ncol(Geno.mtx))
#' outPOLMM = POLMM(objNull, Geno.mtx, chrVec)
#' 
#' outPOLMM
#' round(as.numeric(outPOLMM$pval.spa),2)
#' ## [1] 0.89 0.46 0.82 0.71 0.34 0.30 0.20 0.82 0.25 0.71 # using dense GRM
#' ## [1] 0.82 0.46 0.76 0.68 0.36 0.23 0.21 0.80 0.24 0.71 # using sparse GRM
#' 
#' outPOLMM = POLMM(objNull, Geno.mtx, chrVec, SPAcutoff = 0) # smaller SPAcutoff: more SNPs use saddlepoint approximation
#' round(as.numeric(outPOLMM$pval.spa),2)
#' 

POLMM = function(objNull,
                 Geno.mtx,            # n x m matrix
                 chrVec,
                 SPAcutoff = 2,
                 minMAF = 0.0001,
                 maxMissing = 0.15,
                 impute.method = "fixed",
                 G.model = "Add")
{
  if(!is.element(impute.method,c("fixed"))) 
    stop("Argument 'impute.method' should be 'fixed'.")
  
  if(class(objNull)!="POLMM_NULL_Model")
    stop("class(objNull) should be 'POLMM_NULL_Model'")
  
  subjIDs = rownames(Geno.mtx)
  SNPIDs = colnames(Geno.mtx)
  
  if(is.null(subjIDs) | is.null(SNPIDs))
    stop("rownames and colnames of Geno.mtx are required!")
  if(any(duplicated(subjIDs)))
    stop("all elements of rownames(Geno.mtx) should be unique!")
  
  subjIDs_Null = objNull$subjIDs
  
  if(any(subjIDs_Null != subjIDs)){
    subjPos = match(subjIDs_Null, subjIDs, 0)
    if(any(subjPos==0))
      stop("all subjucts in null model should be also in 'Geno.mtx'.")
    Geno.mtx = Geno.mtx[subjPos,]
  }
  
  n = nrow(Geno.mtx);    # number of individuals
  m = ncol(Geno.mtx);    # number of markers to test
  
  if(missing(chrVec)){
   if(names(objNull$LOCOList)!="LOCO=F")
     stop("chrVec is required unless LOCO = F.")
    chrVec = rep("LOCO=F",m)
  }
  
  if(length(chrVec)==1){
    message(paste0("Chromosome of all markers in analysis are ",chrVec,"."))
    chrVec = rep(chrVec, m)
  }else{
    if(length(chrVec)!=m)
      stop("length of chrVec should be 1 or equal to ncol(Geno.mtx)!")
  }
  
  chrVecLOCO = names(objNull$LOCOList)
  if(any(!is.element(chrVec, chrVecLOCO)))
    stop("all elements in chrVec should be from names(objNull$LOCOList)!")
  
  cat("Totally", m, "markers to analyze!!\n")
  
  uniq_chr = unique(chrVec)
  J = ncol(objNull$muMat)
  yMat = getyMatR(objNull$yVec, n, J)
  
  OutMat = c()
  for(chr in uniq_chr){
    r = objNull$LOCOList[[chr]]$VarRatio
    muMat = objNull$LOCOList[[chr]]$muMat
    iRMat = objNull$LOCOList[[chr]]$iRMat
    muMat1 = muMat[,-1*J]
    
    ## the following are the same for one chromosome
    objP = getobjP(objNull$Cova, yMat, muMat, iRMat)
    XXR_Psi_RX_new = objP[["XXR_Psi_RX_new"]]
    XR_Psi_R_new = objP[["XR_Psi_R_new"]]
    RymuVec = objP[["RymuVec"]]
    RPsiR = objP[["RPsiR"]]
    ##
    
    pos = which(chrVec == uniq_chr)
    K1roots = c(0,0)
    
    for(i in pos){
      GVec = Geno.mtx[,i]
      
      MAF = AF = mean(GVec, na.rm=T)/2
      
      ### genotype imputation
      pos.na = which(is.na(GVec))
      missing.rate = length(pos.na)/n
      if(missing.rate != 0){
        if(impute.method=="fixed")
          GVec[pos.na] = 2 * AF
      }
      
      ### additive / dominant / recessive
      if(G.model=="Add"){}   # do nothing if G.Model is "Add"
      if(G.model=="Dom") GVec = ifelse(GVec >= 1, 1, 0)
      if(G.model=="Rec") GVec = ifelse(GVec <= 1, 0, 1)
      
      SNPID = colnames(Geno.mtx)[i]
      
      if(AF > 0.5){
        MAF = 1 - AF;
        GVec = 2 - GVec
      }
      
      if(MAF < minMAF || missing.rate > maxMissing){
        OutMat = rbind(OutMat, 
                       c(SNPID, chr, MAF, missing.rate, rep(NA,6)))
        next
      }

      #############
      
      posG1 = which(GVec != 0)
      posG0 = setdiff(1:n, posG1)
      
      XR_Psi_RG1 = XR_Psi_R_new[,posG1,drop=F] %*% GVec[posG1]
      adjGVec = GVec - XXR_Psi_RX_new %*% XR_Psi_RG1
      Stat = sum(adjGVec * RymuVec)
      VarWVec = RPsiR * adjGVec^2
      VarW = sum(VarWVec)
      VarW0 = sum(VarWVec[posG0])
      Ratio0 = VarW0/VarW
      VarP = VarW * r
      z = abs(Stat)/sqrt(VarP)
      beta = Stat / VarP
      pval.spa = pval.norm = 2*pnorm(-1*z, lower.tail=T)
      
      if(z > SPAcutoff){
        res.spa <- fastSaddle_Prob(Stat, VarP, VarW, Ratio0, K1roots,
                                   adjGVec[posG1], muMat1[posG1,], iRMat[posG1,])
        pval.spa = res.spa$pval
        K1roots = res.spa$K1roots;
      }
      
      OutMat = rbind(OutMat, 
                     c(SNPID, chr, MAF, missing.rate, Stat, VarW, VarP, beta, pval.norm, pval.spa))
    }
  }
  
  colnames(OutMat) = c("SNPID", "chr", "MAF", "missing.rate", "Stat", "VarW", "VarP", "beta", "pval.norm", "pval.spa")
  OutMat = as.data.frame(OutMat, stringsAsFactors=F)
  return(OutMat)
}


imputeGMat = function(GMatRatio)
{
  m = ncol(GMatRatio)
  GMatOut = c()
  for(j in 1:m){
    GVec = as.numeric(GMatRatio[,j])
    if(any(is.na(GVec)))
      GVec[is.na(GVec)] = mean(GVec, na.rm = T)
    
    GMatOut = cbind(GMatOut, GVec)
  }
  return(GMatOut)
}


## add partial normal approximation to speed up the SPA
fastSaddle_Prob = function(Stat,
                           VarP,
                           VarW,
                           Ratio0,  # Ratio of variance (G==0)
                           K1roots,
                           adjGVec, # n1 x 1, where n1 is length(G!=0)
                           muMat1,  # n1 x (J-1)
                           iRMat)   # n1 x (J-1)
{
  adjStat = Stat / sqrt(VarP)
  n1 = nrow(muMat1)
  J = ncol(muMat1) + 1
  adjGMat = matrix(adjGVec, n1, J-1)
  
  cMat <- adjGMat / iRMat / sqrt(VarW)
  m1 <- sum(muMat1 * cMat)
  
  out.uni1 <- fastgetroot_K1(abs(adjStat), min(K1roots[1],5), Ratio0, muMat1, cMat, m1)
  
  out.uni2 <- fastgetroot_K1(-1*abs(adjStat), max(K1roots[2],-5), Ratio0, muMat1, cMat, m1)
  
  converge = F;
  if(out.uni1$converge == TRUE & out.uni2$converge == TRUE){
    p1 <- fastGet_Saddle_Prob(abs(adjStat), out.uni1$root, out.uni1$K2_eval, Ratio0, muMat1, cMat, m1, lower.tail=FALSE)
    
    p2 <- fastGet_Saddle_Prob(-1*abs(adjStat), out.uni2$root, out.uni2$K2_eval, Ratio0, muMat1, cMat, m1, lower.tail=TRUE)
    
    pval = p1 + p2;
    converge = T;
    K1roots = c(out.uni1$root, out.uni2$root)
  }else{
    print("SPA does not converge, use normal approximation p value.")
    pval = 2*pnorm(-1*abs(adjStat), lower.tail=T)
  }
  
  return(list(pval=pval, converge=converge, K1roots=K1roots))
}

fastgetroot_K1 = function(Stat,
                          init.t,
                          Ratio0,
                          muMat,
                          cMat,
                          m1,
                          tol = .Machine$double.eps^0.25,
                          maxiter = 100)
{
  t = init.t;
  K1_eval = 0
  diff.t = Inf
  converge = T
  
  for(iter in 1:maxiter){
    old.t = t
    old.diff.t = diff.t
    old.K1 = K1_eval
    
    K12_eval = K12(t, muMat, cMat, m1)
    K1_eval = K12_eval[1,1] - Stat + Ratio0 * t
    K2_eval = K12_eval[1,2] + Ratio0
    
    diff.t = -1 * K1_eval / K2_eval
    if(is.na(K1_eval) | is.infinite(K1_eval)){
      # checked it on 07/05:
      # if the solution 't' tends to infinity, 'K2_eval' tends to 0, and 'K1_eval' tends to 0 very slowly.
      # then we can set the one sided p value as 0 (instead of setting converge = F)
      t = sign(Stat)*Inf
      K2_eval = 0;
      break;
    }
    
    if(sign(K1_eval) != sign(old.K1)){
      while(abs(diff.t) > abs(old.diff.t) - tol){
        diff.t = diff.t/2
      }
    }
    if(abs(diff.t) < tol) break;
    t = old.t + diff.t
  }
  
  if(iter == maxiter) converge = F
  
  return(list(root = t,
              iter = iter,
              converge = converge,
              K2_eval = K2_eval))
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

K12 = function(t,
               muMat,
               cMat,
               m1)
{
  n = length(t)
  outMat = matrix(0, n, 2)
  
  for(i in 1:n){
    t1 = t[i]
    temp0Mat = muMat * exp(cMat * t1)
    temp1Mat = - muMat + temp0Mat
    temp2Mat = temp0Mat * cMat
    temp3Mat = temp2Mat * cMat
    
    temp1 = 1 + rowSums(temp1Mat)
    temp2 = rowSums(temp2Mat)
    temp3 = rowSums(temp3Mat)
    
    outMat[i,1] = sum(temp2/temp1) - m1
    outMat[i,2] = sum((temp3*temp1-temp2^2)/temp1^2, na.rm=TRUE)
  }
  return(outMat)
}


fastGet_Saddle_Prob = function(Stat,
                               zeta,
                               K2_eval,
                               Ratio0,
                               muMat,
                               cMat,
                               m1,          # sum(muMat * cMat)
                               lower.tail)
{
  k1 = Korg(zeta, muMat, cMat, m1) + 1/2 * zeta^2 * Ratio0
  k2 = K2_eval
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






