
#' Test for association between genetic variants and an ordinal categorical variable via Proportional Odds Logistic Mixed Model (POLMM-Gene)
#' 
#' Test for association between genetic variants and an ordinal categorical variable via Proportional Odds Logistic Mixed Model (POLMM-Gene)
#' 
#' @param objNull the output of the POLMM_Null_Model() function 
#' @param Geno.mtx a numeric genotype matrix with each row as a subject and each column as a marker in a region or gene. 
#'                 Column names of marker IDs and row names of individual IDs are required.
#'                 Missng genotype should be coded as in argument 'G.missing'. Both hard-called and imputed genotype are supported.
#' @param chrom a character to specify chromosome of the markers in Geno.mtx. Must be specified unless LOCO = F.
#' @param SPAcutoff a standard deviation cutoff (default=2). If the standardized test statistic < SPAcutoff, normal approximation is used, otherwise, saddlepoint approximation is used.
#' @param maxMAF a cutoff of the minimal minor allele frequencies (MAFs). Any markers with MAF < minMAF will be excluded from the analysis.
#' @param maxMissing a cutoff of the maximal missing rate. Any markers with missing rate > maxMissing will be excluded from the analysis.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#'                      "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e. 2p where p is MAF).
#' @param G.model a character string (default: "Add") to specify additive ("Add"), dominant ("Dom"), or recessive ("Rec") model. 
#'                If "Dom", GVec = ifelse(GVec >= 1, 1, 0), if "Rec", GVec = ifelse(GVec <= 1, 0, 1). Be very careful if the gneotyp is imputed data.
#' @param G.missing the code for missing genotype (default=NA). For plink input, G.missing = -9.
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
#' \item{switch.allele}{a logical value indicating if the REF/ALT alleles were switched, if AF > 0.5, we use GVec = 2-GVec, and then give switch.allele=T. This is useful to estimate the effect direction.}
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
#' Geno.mtx = matrix(rbinom(10000,2,0.3),1000,10)
#' rownames(Geno.mtx) = egData$IID
#' colnames(Geno.mtx) = paste0("rs",1:10)
#' chrVec = chrom = "1"  # equivalant to chrVec = rep("1", ncol(Geno.mtx))
#' outPOLMM = POLMM(objNull, Geno.mtx, chrVec)
#' 
#' outList = POLMM.Gene(objNull, Geno.mtx, chrom, SparseGRM)
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
                      Geno.mtx,            # n x q matrix, where n is number of subjects and m is number of markers
                      chrom,
                      SparseGRM,
                      SKAT.control = NULL,
                      SPAcutoff = 2,
                      G.model = "Add")
{
  n = nrow(Geno.mtx)
  # check and use the default setting of SKAT.control
  SKAT.control = check.SKAT.control(SKAT.control)
  
  GMat.list = Check_GMat(Geno.mtx, SetID, 
                         kernel = SKAT.control$kernel, 
                         weights.beta = SKAT.control$weights.beta, 
                         impute.method = SKAT.control$impute.method, 
                         missing_cutoff = SKAT.control$missing_cutoff, 
                         max_maf = SKAT.control$max_maf)
  
  if(GMat.list$error){
    outList = list(p.value = NA, param = NA, p.value.resampling = NA)
    return(outList)
  }
  
  # extract information from output of Check_GMat()
  Geno.mtx = GMat.list$Geno.mtx
  weights = GMat.list$weights
  AlleleFlip.Vec = GMat.list$AlleleFlip.Vec
  MAF.Vec = GMat.list$MAF.Vec
  MAC.Vec = MAF.Vec * 2 * n   # this might be slightly different from "true" MAC due to the genotype imputation, but should be OK since we have limited the missing_rate
  
  # update SparseGRM and set an objective for Gene-based analysis
  SparseGRM = updateSparseGRM(SparseGRM, objNull$subjIDs)
  setPOLMMGENEobj(objNull$controlList$maxiterPCG, 
                  objNull$controlList$tolPCG,
                  objNull$Cova,
                  objNull$yVec,
                  objNull$tau,
                  SparseGRM,
                  objNull$LOCOList,
                  objNull$eta)
  
  setPOLMMGENEchr(objNull$LOCOList, 
                  chrom)
  
  ##
  
  OutList = getStatVarS(Geno.mtx)
  
  StatVec = as.vector(OutList$StatVec)
  VarSMat = OutList$VarSMat
  
  VarSVec = diag(VarSMat)
  
  out_SPA_ER = adj_SPA_ER(Geno.mtx, MAC.Vec, StatVec,  VarSVec, SPA.pval.cutoff = 0, ER.MAC.cutoff = 0) # first check the non-robust version, no adjustment
  adjVarSVec = out_SPA_ER$adjVarSVec
  
  adjPVec = out_SPA_ER$adjPVec
  
  ### add something about weights
  
  weights = Get_Weights(SKAT.control$kernel, MAF.Vec, SKAT.control$weights.beta)
  
  #########
  
  wStatVec = StatVec * weights
  
  if(any(is.infinite(adjVarSVec))) stop("any(is.infinite(adjVarSVec))") # I want to know when this will happen
  adjVarSVec = ifelse(is.infinite(adjVarSVec), 0, adjVarSVec)
  StatVec = ifelse(is.infinite(adjVarSVec), 0, StatVec)
  
  r0 = adjVarSVec / VarSVec  # adjVarSVec might be 0
  wr0 = sqrt(r0) * weights
  
  wadjVarSMat = t(VarSMat * wr0) * wr0
  
  ## might need another step to introduce another ratio to adjust for variance matrix based on Burden Test
  
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
  
  # OutMat = as.data.frame(OutMat, stringsAsFactors=F)
  return(outList)
}

adj_SPA_ER = function(Geno.mtx,
                      MAC.Vec,   # imputation should be finished before this step
                      StatVec,
                      VarSVec,
                      SPA.pval.cutoff,
                      ER.MAC.cutoff)
{
  q = ncol(Geno.mtx)
  
  QVec = StatVec^2 / VarSVec  # follow a standard chi-squre distribution
  adjPVec = PVec = pchisq(QVec, lower.tail = FALSE, df = 1)
  adjVarSVec = VarSVec
  
  pos.ER = which(MAC.Vec < ER.MAC.cutoff)
  pos.SPA = which(MAC.Vec >= ER.MAC.cutoff & PVec < SPA.pval.cutoff)
  
  for(i in pos.ER){
    adjPVec[i] = getPvalER(Geno.mtx[,i])
    adjVarSVec[i] = StatVec[i]^2 / qchisq(adjPVec[i], df = 1, lower.tail = FALSE)
  }
  
  for(i in pos.SPA){
    adjPVec[i] = getPvalSPA(Geno.mtx[,i])
    adjVarSVec[i] = StatVec[i]^2 / qchisq(adjPVec[i], df = 1, lower.tail = FALSE)
  }
  
  adjVarSVec = ifelse(adjPVec == 0, StatVec^2 / 500, adjVarSVec)
  
  out_SPA_ER = list(adjVarSVec = adjVarSVec, adjPVec = adjPVec)
  
  return(out_SPA_ER)
}

getPvalSPA = function()
{
  # to be continued
}

getPvalER = function()
{
  # to be continued
}

####### ---------- Check Argument ------------ #########

check.SKAT.control = function(SKAT.control)
{
  # the below is the default setting
  default.SKAT.control = list(kernel = "linear.weighted",
                              method = "SKATO",
                              weights.beta = c(1,25),
                              weights = NULL,
                              impute.method = "bestguess",
                              r.corr = NULL,
                              is_check_genotype=TRUE,
                              is_dosage = FALSE, 
                              missing_cutoff = 0.15, 
                              max_maf = 1, 
                              estimate_MAF = 1)
  
  # use the default setting or update them
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
    stop("'kernel' should be 'linear' or 'linear.weighted'. ")
  
  if(length(SKAT.control$weights.beta)!=2)
    stop("length of 'weights.beta' should be 2.")
  
  if(! SKAT.control$impute.method %in% c("fixed", "bestguess", "random"))
    stop("'impute.method' should be 'fixed','bestguess', or 'random'. ")
  
  if(SKAT.control$missing_cutoff > 1 | SKAT.control$missing_cutoff < 0)
    stop("'missing_cutoff' should be between 0 and 1.")
  
  if(SKAT.control$max_maf > 1 | SKAT.control$max_maf <= 0)
    stop("'max_maf' should be between 0 and 1.")
  
  method = SKAT.control$method
  if(! method %in% c("SKAT", "Burden", "SKATO"))
    stop("'method' should be 'SKAT', 'Burden', or 'SKATO'.")
  
  if(is.null(SKAT.control$r.corr)){
    if(method == "SKAT")
      SKAT.control$r.corr = 0;
    
    if(method == "Burden")
      SKAT.control$r.corr = 1;
    
    if(method == "SKATO")
      SKAT.control$r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1);
  }
  
  if(any(SKAT.control$r.corr < 0 | SKAT.control$r.corr > 1))
    stop("'r.corr should be between 0 and 1.'")
  
  return(SKAT.control)
}

####### ---------- Get Weights from MAF ---------- #######

Get_Weights = function(kernel, MAF, weights.beta)
{
  if(kernel == "linear"){
    weights = rep(1, length(MAF)) 
  }
  
  if(kernel == "linear.weighted"){
    weights = dbeta(MAF, weights.beta[1], weights.beta[2])
  }
  
  return(weights)
}

####### ---------- Check the Geno.mtx matrix, and do imputation ---------- #######

Check_GMat = function(Geno.mtx, SetID, kernel, weights.beta, impute.method, missing_cutoff, max_maf)
{
  # Number of subjects and SNPs
  n = nrow(Geno.mtx)
  q = ncol(Geno.mtx)
  
  Geno.mtx[Geno.mtx > 1.8] = 2
  Geno.mtx[Geno.mtx < 0.2] = 0
  Geno.mtx[Geno.mtx == 9] = NA
  
  MissingRate.Vec = colMeans(is.na(Geno.mtx))
  MAF.Vec = colMeans(Geno.mtx, na.rm = T) / 2
  Var.Vec = rowMeans((t(Geno.mtx) - 2 * MAF.Vec)^2)
  
  AlleleFlip.Vec = ifelse(MAF.Vec > 0.5, T, F)
  MAF.Vec = ifelse(MAF.Vec > 0.5, 1 - MAF.Vec, MAF.Vec)
  
  pos.PassG = which(MissingRate.Vec <= missing_cutoff & MAF.Vec <= max_maf & Var.Vec !=0)
  
  if(length(pos.PassG) == 0){
    msg = sprintf("In %s, ALL SNPs have been excluded. P value = 1",SetID)
    warning(msg, call.=FALSE)
    
    re = list(error = 1) 
    return(re)
  }
  
  MAF.Vec = MAF.Vec[pos.PassG]
  AlleleFlip.Vec = AlleleFlip.Vec[pos.PassG]
  MissingRate.Vec = MissingRate.Vec[pos.PassG]
  
  Geno.mtx = Geno.mtx[,pos.PassG,drop=F]
  
  # genotype imputation
  pos.imputeG = which(MissingRate.Vec > 0)
  
  if(length(pos.imputeG) > 0){
    for(i in pos.imputeG){
      GVec = Geno.mtx[,i]
      MAF = MAF.Vec[i]
      IDX_NA = which(is.na(GVec))
      n_NA = length(IDX_NA)
      #
      if(impute.method == "fixed") 
        GVec[IDX_NA] = 2 * MAF
      if(impute.method == "bestguess") 
        GVec[IDX_NA] = 0  # ifelse(MAF < 0.5, 0, 2). Since we only consider low-frequency variants and have flipped the allele
      if(impute.method == "random")
        GVec[IDX_NA] = rbinom(n_NA, 2, maf)
    }
  }
  
  # get weights for rare variants analysis
  weights = Get_Weights(kernel, MAF.Vec, weights.beta)
  
  return(list(Geno.mtx = Geno.mtx, 
              weights = weights,
              AlleleFlip.Vec = AlleleFlip.Vec,
              MAF.Vec = MAF.Vec, 
              error = 0))
}

