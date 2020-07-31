
#' Please run this function before running getSparseGRM(). 
#' 
#' Please run this function before running getSparseGRM(). We strongly suggest using parallel computing for different pairs of (chrParallel, partParallel).
#' @param chrParallel a character of chromosome number (e.g. 1, 2, ..., 22), better autosomes.  
#' @param partParallel a numeric value (from 1 to nPartsGRM)
#' @param nPartsGRM a numeric value (e.g. 250): GCTA software can split subjects to multiple parts. For UK-Biobank analysis, it is recommanded to use 250 parts. 
#' @param PlinkFile a path to Plink files. The current version (gcta_1.93.1beta) of gcta software does not support different prefix names for bim, bed and fam files. 
#' @param tempDir a path to store temp files from getSparseGRMParallel(). This should be consistent to the input of getSparseGRM(). Default is system.file("SparseGRM", "temp", package = "POLMM").
#' @param subjData a character vector to specify subject IDs to retain (i.e. IID). Default is NULL, i.e. all subjects are retained in sparse GRM. If the number of subjects is less than 1,000, the GRM estimation might not be accurate.
#' @param minMafGRM Minimal value of MAF cutoff to select markers (from Plink files) to construct GRM.
#' @param maxMissingGRM Maximal value of missing rate to select markers (from Plink files) to construct GRM.
#' @param threadNum Number of threads (CPUs) to use.
#' @examples 
#' ## Please check help(getSparseGRM) for example codes.
#' @export
 
getSparseGRMParallel = function(chrParallel,
                                partParallel,
                                nPartsGRM,
                                PlinkFile,
                                tempDir = NULL,
                                subjData = NULL,
                                minMafGRM = 0.01,
                                maxMissingGRM = 0.1,
                                threadNum = 8)
{
  bimFile = paste0(PlinkFile,".bim")
  bedFile = paste0(PlinkFile,".bed")
  famFile = paste0(PlinkFile,".fam")
  
  if(!file.exists(bimFile)) stop("Could not find bimFile or paste0(PlinkFile,'.bim')")
  if(!file.exists(bedFile)) stop("Could not find bedFile or paste0(PlinkFile,'.bed')")
  if(!file.exists(famFile)) stop("Could not find famFile or paste0(PlinkFile,'.fam')")
  
  gcta64File = system.file("gcta_1.93.1beta", "gcta64", package = "POLMM");
  system(paste("chmod +x",gcta64File))
  
  if(is.null(tempDir)) tempDir = system.file("SparseGRM", "temp", package = "POLMM")
  
  PlinkName = basename(PlinkFile)
  tempFile = paste0(tempDir, "/Plink-",PlinkName,"-chr-",chrParallel,"-minMaf-",minMafGRM,"-maxMissing-",maxMissingGRM)
  
  cmd = paste(gcta64File, 
              "--bfile", PlinkFile,
              "--out", tempFile,
              "--chr", chrParallel,
              "--make-grm-part", nPartsGRM, partParallel,
              "--maf", minMafGRM,
              "--geno", maxMissingGRM,
              "--thread-num",threadNum)
  
  ## only retain parts of subjects
  if(!is.null(subjData)){
    if(length(subjData) < 1000) stop("length(subjData) < 1000, the MAF estimate might be inaccurate.")
    subjFile = paste0(tempDir,"/subjData.txt")
    famData = read.table(famFile)
    posSubj = match(subjData, famData$V2, 0)
    if(any(posSubj==0)) stop("All subjects in 'subjData' should be in IID column of 'famFile'.")
    
    write.table(famData[posSubj,c(1,2)], 
                subjFile, row.names=F, col.names=F, quote=F)
    cmd = paste(cmd, 
                "--keep", subjFile)
  }
  
  system(cmd)
}


#' Get an object of 'SparseGRM' for POLMM_Null_Model(). Only valid in Linux since GCTA software can only be used in Linux.
#' 
#' If the sample size is greater than 100,000, we recommend using sparse GRM to adjust for sample relatedness. 
#' This function is to use GCTA software (gcta_1.93.1beta, https://cnsgenomics.com/software/gcta/#Overview) to get an object of 'SparseGRM' to be passed to function POLMM_Null_Model(). 
#' \cr\cr
#' Step 1: Run getSparseGRMParallel(); please check help(getSparseGRMParallel) for more details. 
#' \cr
#' Step 2: Run getSparseGRM().
#' @param chrVec a character vector to specify all chromosomes (e.g. 1:22), better autosomes. This should be consistent to the input of getSparseGRMParallel(). 
#' @param PlinkFile a path to Plink files. The current version (gcta_1.93.1beta) of gcta software does not support different prefix names for bim, bed and fam files. 
#' @param nPartsGRM a numeric value (e.g. 250): GCTA software can split subjects to multiple parts. For UK-Biobank analysis, it is recommanded to use 250 parts. 
#' @param tempDir a path to store temp files from getSparseGRMParallel(). This should be consistent to the input of getSparseGRMParallel(). Default is system.file("SparseGRM", "temp", package = "POLMM").
#' @param relatednessCutoff a cutoff for sparse GRM, only kinship coefficient greater than this cutoff will be retained in sparse GRM.
#' @param minMafGRM Minimal value of MAF cutoff to select markers (from Plink files) to construct GRM.
#' @param maxMissingGRM Maximal value of missing rate to select markers (from Plink files) to construct GRM.
#' @param rm.tempFiles a logical value indicating if the temp files generated in getSparseGRMParallel will be deleted
#' @examples 
#' ## Input data:
#' famFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.fam", package = "POLMM")
#' PlinkFile = gsub(".fam", "", famFile)   # fam/bim/bed files should have the same prefix
#' nPartsGRM = 2;   # nPartsGRM = 250 for UK Biobank data analysis
#' chrVec = 1:22    # maybe paste0("chr",1:22), depending on the plink bim file
#' 
#' ## Step 1:
#' ## We strongly suggest parallel computing in high performance clusters (HPC) for different pairs of (chrParallel, partParallel). 
#' for(chrParallel in chrVec){
#'   for(partParallel in 1:nPartsGRM){
#'     getSparseGRMParallel(chrParallel, partParallel, nPartsGRM, PlinkFile)
#'   }
#' }
#' 
#' ## After step 1, in "tempDir (default: system.file("SparseGRM", "temp", package = "POLMM"))", there will be results (might needs large amount of space) corresponding to different pairs of (chrParallel, partParallel).
#' 
#' ## Step 2:
#' ## Combine results in step 1 to calculate an object with class of SparseGRM for POLMM_Null_Model(),
#' SparseGRM = getSparseGRM(chrVec, PlinkFile, nPartsGRM)
#' 
#' ## NOTE: You can change some options such as (minMafGRM, maxMissingGRM, nPartsGRM), but keep in mind that functions getSparseGRMParallel() and getSparseGRM() should use the same options.
#' @export

getSparseGRM = function(chrVec,
                        PlinkFile,
                        nPartsGRM,
                        tempDir = NULL,
                        relatednessCutoff = 0.05,
                        minMafGRM = 0.01,
                        maxMissingGRM = 0.1,
                        rm.tempFiles = F)
{
  PlinkName = basename(PlinkFile)
  chrVec = as.character(chrVec)
  nDigits = floor(log10(nPartsGRM)) + 1 # 1-9 -> 1; 10-99 -> 2; 100:999 -> 3.
  
  AllIDs = c()
  n0 = 0
  SparseGRM = list()
  
  if(is.null(tempDir)) tempDir = system.file("SparseGRM", "temp", package = "POLMM")
  
  ## cycle for nPartsGRM
  for(i in 1:nPartsGRM){
    print(paste("Analyzing part",i,"of total",nPartsGRM,"parts."))
    tempList = list()
    for(chr in chrVec){
      # print(paste("Analyzing chromosome",chr))
      tempFile = paste0(tempDir, "/Plink-",PlinkName,"-chr-",chr,"-minMaf-",minMafGRM,"-maxMissing-",maxMissingGRM)
      ## Three files generated by GCTA
      IDFile = paste0(tempFile,".part_",nPartsGRM,"_",formatC(i, width = nDigits, flag = '0'),".grm.id")
      BinFile = paste0(tempFile,".part_",nPartsGRM,"_",formatC(i, width = nDigits, flag = '0'),".grm.bin")
      NFile = paste0(tempFile,".part_",nPartsGRM,"_",formatC(i, width = nDigits, flag = '0'),".grm.N.bin")
      ## read in the three files
      IDs = read.table(IDFile, stringsAsFactors=F)
      ID = IDs$V2
      n1 = n0 + length(ID)
      nData = (n1 - n0) * (n0 + n1 + 1) / 2
      
      grm = readBin(BinFile, n = nData, what = numeric(0), size = 4)
      nMarkers = readBin(NFile, n = nData, what = numeric(0), size = 4)
      tempList[[chr]] = list(ID = ID, grm = grm, nMarkers = nMarkers)
      
      if(rm.tempFiles){
        file.remove(IDFile, BinFile, NFile)
      }
    }
    
    ## calculate GRM based on all chromosomes
    AllIDs = c(AllIDs, ID)
    grm.all = rep(0, nData)
    nMarkers.all = rep(0, nData)
    for(chr in chrVec){
      tempChr = tempList[[chr]]
      grm.all = grm.all + tempChr$grm * tempChr$nMarkers
      nMarkers.all = nMarkers.all + tempChr$nMarkers
    }
    grm.all = grm.all / nMarkers.all
    
    pos = which(grm.all > relatednessCutoff)
    value = grm.all[pos]
    pairs = getPairs(pos, value, AllIDs, n0, n1)
    SparseGRM[["none"]] = rbind(SparseGRM[["none"]], pairs)
    
    ## Leave-One-Chromosome-Out (LOCO)
    if(length(chrVec)>1){
      for(chr in chrVec){
        tempChr = tempList[[chr]]
        grm.LOCO = (grm.all * nMarkers.all - tempChr$grm * tempChr$nMarkers) / (nMarkers.all - tempChr$nMarkers);
        pos = which(grm.LOCO > relatednessCutoff)
        value = grm.LOCO[pos]
        pairs = getPairs(pos, value, AllIDs, n0, n1)
        SparseGRM[[chr]] = rbind(SparseGRM[[chr]], pairs)
      }
    }
    
    n0 = n1 
  }
  class(SparseGRM) = "SparseGRM"
  return(SparseGRM)
}


getPairs = function(pos, value, ID, n0, n1){
  start = c()
  end = c()
  temp = 0
  for(i in (n0+1):n1){
    start = c(start, temp)
    temp = temp + i;
    end = c(end, temp);
  }
  
  pos.row = sapply(pos, FUN=function(x){min(which(x<=end))})
  pos.col = pos - start[pos.row]
  
  ID1 = ID[n0 + pos.row];
  ID2 = ID[pos.col];
  pairs = cbind.data.frame(ID1, ID2, value);
  return(pairs)
}





