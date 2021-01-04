
check.POLMM.control = function(POLMM.control = NULL, 
                               which.step = NULL)          # "step1" or "step2"
{
  # argument of 'which.step' can be only "step1" or "step2"
  if(is.null(which.step) | !which.step %in% c("step1", "step2")){
    stop("Argument of 'which.step' should be 'step1' or 'step2'.")
  }
  
  # default setting of POLMM.control
  if(which.step == "step1"){
    # to be continued
  }
  
  if(which.step == "step2"){
    default.POLMM.control = list(impute_method = "fixed",  # the below is shared parameters for both single marker testing and region-based testing
                                 missing_cutoff = 0.15,
                                 min_maf_marker = 0.001,
                                 min_mac_marker = 20,
                                 max_maf_region = 0.01,
                                 nMarkers_output = 10000,
                                 SPAcutoff = 2,
                                 memory_chunk = 4,
                                 kernel = "linear.weighted",   # the below is parameters for region-based testing
                                 method = "SKAT-O",
                                 weights_beta = c(1,25),
                                 r_corr = NULL,
                                 printPCGInfo = FALSE,
                                 tolPCG = 1e-5,
                                 maxiterPCG = 100)
  }
  
  # use the default setting or update it
  if(!is.null(POLMM.control)){
    if(class(POLMM.control) != "list"){
      stop("Argument of 'POLMM.control' should be 'list'.")
    }
    ctrl.nm = names(POLMM.control)
    for(nm in ctrl.nm){
      default.POLMM.control[[nm]] = POLMM.control[[nm]]
    }
  }
      
  POLMM.control = default.POLMM.control
  
  # check the parameters
  
  if(which.step == "step1")
  {
    # to be continued
  }
  
  if(which.step == "step2")
  {
    # used in both single-marker analysis and region-based analysis
    if(! POLMM.control$impute_method %in% c("fixed", "bestguess", "random"))
      stop("'impute.method' should be 'fixed','bestguess', or 'random'. Check 'Details' for more details.")
    
    if(POLMM.control$missing_cutoff > 1 | POLMM.control$missing_cutoff < 0)
      stop("'missing_cutoff' should be between 0 and 1. The default setting is 0.15.")
    
    if(POLMM.control$SPA_cutoff < 0)
      stop("'SPA_cutoff' should be greater than or equal to 0. The default setting is 2. Check 'Details' for more details.")
    
    # only used in single-marker analysis
    if(POLMM.control$min_maf_marker > 1 | POLMM.control$min_maf_marker <= 0)
      stop("'min_maf_marker' should be between 0 and 1. The default setting is 0.001.")
    
    if(POLMM.control$min_mac_marker <= 5)
      stop("'min_mac_marker' should be greater than 5. The default setting is 20.")
    
    if(POLMM.control$nMarkers_output <= 999)
      stop("nMarkers_output should be greater than or equal to 1000. The default setting is 10000.")
    
    # only used in region-based analysis
    if(POLMM.control$max_maf_region > 1 | POLMM.control$max_maf_region <= 0)
      stop("'max_maf_region' should be between 0 and 1. The default setting is 0.01.")
    
    if(! POLMM.control$kernel %in% c("linear", "linear.weighted"))
      stop("'kernel' should be 'linear' or 'linear.weighted'. Check 'Details' for more details.")
    
    if(length(POLMM.control$weights_beta) != 2 | any(POLMM.control$weights_beta < 0))
      stop("length of 'weights_beta' should be 2. The two elements in 'weights_beta' should be non-negative. Check 'Details' for more details.")
    
    method_region = POLMM.control$method_region
    
    if(! method_region %in% c("SKAT", "Burden", "SKAT-O"))
      stop("'method' should be 'SKAT', 'Burden', or 'SKAT-O'. Check 'Details' for more details.")
    
    if(is.null(POLMM.control$r_corr)){
      if(method_region == "SKAT") POLMM.control$r.corr = 0;
      if(method_region == "Burden") POLMM.control$r.corr = 1;
      if(method_region == "SKAT-O") POLMM.control$r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1);  # r.corr = 0 is SKAT, r.corr = 1 is Burden Test
    }else{
      message("Since 'r_corr' is specified, the 'method_region' is ignored.")
    }
    
    if(any(POLMM.control$r.corr < 0 | POLMM.control$r.corr > 1))
      stop("'r.corr' should be a numeric vector in which each element is between 0 and 1. Check 'Details' for more details.")
  }
  
  return(POLMM.control)
}

setGenoInput = function(GenoFile, 
                        GenoFileIndex, 
                        SampleInModel, 
                        output.marker = F,
                        output.subject = F)
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
    
    marker = subject = NULL
    if(output.marker){
      markerInfo = data.table::fread(bimFile)
      marker = markerInfo$V2
    }
    
  }
  
  # add more types of genotype input later
  
  
  # return genotype
  genoList = list(genoType = genoType,
                  marker = marker,
                  subject = subject)
  
  return(genoList)
}

split.marker = function(all.markers, nMarkers_output)
{
  n = length(all.markers)
  pos.start = seq(1, n, nMarkers_output)
  pos.end = pos.start + nMarkers_output - 1
  
  nChunks = length(pos.start)
  pos.end[nChunks] = min(pos.end[nChunks], n)
  
  marker.list = list()
  for(i in 1:nChunks){
    posMarker = pos.start[i]:pos.end[i]
    marker.list[[i]] = all.markers[posMarker]
  }
  
  return(marker.list)
}

getAnnoList = function(AnnoFile, AnnoHeader)
{
  if(!file.exists(AnnoFile))
    stop(paste("Cannot find AnnoFile in", AnnoFile))
  
  AnnoData = data.table::fread(AnnoFile, header = T, stringsAsFactors = F);
  AnnoData = as.data.frame(AnnoData)
  HeaderInAnnoFile = colnames(AnnoData)
  
  if(any(HeaderInAnnoFile[1:2] != c("GENE","SNP")))
    stop("The first two elements in the header of AnnoFile should be c('GENE', 'SNP').")
  
  if(!is.null(AnnoHeader)){
    if(any(!AnnoHeader %in% HeaderInAnnoFile))
      stop("At least one element in AnnoHeader is not in the header of AnnoFile.")
    posAnno = which(HeaderInAnnoFile %in% AnnoHeader)
  }else{
    posAnno = NULL
  }
  
  AnnoList = list()
  uGENE = unique(AnnoData$GENE)
  for(g in uGENE){
    pos = which(AnnoData$GENE == g)
    SNP = AnnoData$SNP[pos]
    AnnoMat = cbind(All=1, AnnoData[pos, posAnno, drop=F])
    rownames(AnnoMat) = SNP
    if(any(duplicated(SNP)))
      stop(paste0("Please check AnnoFile: in gene ",g,", duplicated SNPs exist."))
    AnnoList[[g]] = list(SNP = SNP,
                         AnnoMat = AnnoMat)
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
