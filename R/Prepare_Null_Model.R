
#' Get an object of 'SparseGRM' for POLMM_Null_Model()
#' 
#' If you prefer using a sparse GRM to adjust for sample relatedness, you need this function to get an object of 'SparseGRM'. 
#' Our simulations show that using sparse GRM is almost the same reliable as using dense GRM, with a huge improvement of computational efficiency. 
#' \cr\cr
#' Step 1: Download GCTA software from https://cnsgenomics.com/software/gcta/#Overview; 
#' \cr
#' Step 2: Run getSparseGRMParallel(); 
#' \cr
#' Step 3: Run getSparseGRM().
#' @param outPrefix a path to specify the prefix of output files  
#' @param chrParallel chromosome number (e.g. 1, 2, ..., 22), better autosomes.  
#' @param partParallel part number (from 1 to nPartsGRM)
#' @param gcta64File a path to GCTA software. Please download from https://cnsgenomics.com/software/gcta/#Overview.
#' @param PlinkFile a path to Plink files. The current version (gcta_1.93.1beta) of gcta software does not support difference prefix names for bim, bed and fam files. 
#' @param subjData a character vector to specify subject IDs (i.e. IID). If the number of subjects is less than 1,000, the GRM estimation is not accurate.
#' @param minMafGRM Minimal value of MAF cutoff to select markers (from Plink files) to construct GRM.
#' @param maxMissingGRM Maximal value of missing rate to select markers (from Plink files) to construct GRM.
#' @param nPartsGRM GCTA software can split subjects to multiple parts. For UK-Biobank analysis, it is recommanded to use 250 parts. 
#' @param threadNum Number of threads (CPUs) to use.
#' @examples 
#' ## First download GCTA software, then you can use the following commands
#' famFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.fam", package = "POLMM")
#' PlinkFile = gsub(".fam", "", famFile)   # fam/bim/bed files should have the same prefix
#' gcta64File = "/net/snowwhite/home/wenjianb/gcta64"  # should download gcta first (https://cnsgenomics.com/software/gcta/#Download, author version is gcta_1.93.1beta)
#' outPrefix = gsub("nSNPs-10000-nsubj-1000-ext.fam","SparseGRM/sp", famFile) # Dir/prefix of the output files (in step 2)
#' SparseGRMFile = gsub("nSNPs-10000-nsubj-1000-ext.fam","SparseGRM.RData", famFile) # output file (in step 3)
#' nPartsGRM = 2;   # nPartsGRM = 250 for UK Biobank data analysis
#' chrVec = 1:4    # maybe paste0("chr",1:22), depending on the plink bim file
#' 
#' ## Step 2:
#' ## We strongly suggest parallel computing for different pairs of (chrParallel, partParallel). 
#' for(chrParallel in chrVec){
#'   for(partParallel in 1:nPartsGRM){
#'     getSparseGRMParallel(outPrefix, chrParallel, partParallel, gcta64File, PlinkFile, nPartsGRM = nPartsGRM)
#'   }
#' }
#' 
#' ## After that, in "outPrefix", there will be results (needs large amount of storage) corresponding to different pairs of (chrParallel, partParallel).
#' 
#' ## Step 3:
#' ## Combine results in step 2 to calculate an object with class of SparseGRM for POLMM_Null_Model(),
#' SparseGRM = getSparseGRM(outPrefix, chrVec, PlinkFile, nPartsGRM = nPartsGRM)
#' save(SparseGRM, file=SparseGRMFile)
#' 
#' ## NOTE: You can change some options such as (minMafGRM, maxMissingGRM, nPartsGRM), but keep in mind that functions getSparseGRMParallel() and getSparseGRM() should use the same change of these options.
getSparseGRMParallel = function(outPrefix,
                                chrParallel,
                                partParallel,
                                gcta64File,
                                PlinkFile,
                                subjData = NULL,
                                minMafGRM = 0.01,
                                maxMissingGRM = 0.1,
                                nPartsGRM = 250,
                                threadNum = 8)
{
  # bimFile = ifelse(!missing(bimFile), bimFile,
  #                  ifelse(!missing(PlinkFile), paste0(PlinkFile,".bim"), stop("Either bimFile or PlinkFile should be specified.")))
  # bedFile = ifelse(!missing(bedFile), bedFile,
  #                  ifelse(!missing(PlinkFile), paste0(PlinkFile,".bed"), stop("Either bedFile or PlinkFile should be specified.")))
  # famFile = ifelse(!missing(famFile), famFile,
  #                  ifelse(!missing(PlinkFile), paste0(PlinkFile,".fam"), stop("Either famFile or PlinkFile should be specified.")))
  
  bimFile = paste0(PlinkFile,".bim")
  bedFile = paste0(PlinkFile,".bed")
  famFile = paste0(PlinkFile,".fam")
  
  if(!file.exists(bimFile)) stop("Could not find bimFile or paste0(PlinkFile,'.bim')")
  if(!file.exists(bedFile)) stop("Could not find bedFile or paste0(PlinkFile,'.bed')")
  if(!file.exists(famFile)) stop("Could not find famFile or paste0(PlinkFile,'.fam')")
  
  if(!file.exists(gcta64File)) stop("Could not find gcta64File, please download it from https://cnsgenomics.com/software/gcta/#Download. My version is gcta_1.93.1beta.")
  
  PlinkName = basename(PlinkFile)
  outFile = paste0(outPrefix, "-Plink-",PlinkName,"-chr-",chrParallel,"-minMaf-",minMafGRM,"-maxMissing-",maxMissingGRM)
  
  cmd = paste(gcta64File, 
              "--bfile", PlinkFile,
              "--out", outFile,
              "--chr", chrParallel,
              "--make-grm-part", nPartsGRM, partParallel,
              "--maf", minMafGRM,
              "--geno", maxMissingGRM,
              "--thread-num",threadNum)
  
  ## only retain parts of subjects
  if(!is.null(subjData)){
    if(length(subjData) < 1000) stop("length(subjData) < 1000, the MAF estimate might be inaccurate.")
    subjFile = paste0(outPrefix,".subj")
    famData = read.table(famFile)
    posSubj = match(subjData, famData$V2, 0)
    if(any(posSubj==0)) stop("All subjects in 'subjData' should be in IID column of 'PlinkFile'.")
    
    write.table(famData[posSubj,c(1,2)], 
                subjFile, row.names=F, col.names=F, quote=F)
    cmd = paste(cmd, 
                "--keep", subjFile)
  }
  
  system(cmd)
}


#' Get an object of 'SparseGRM' for POLMM_Null_Model()
#' 
#' If you prefer using a sparse GRM to adjust for sample relatedness, you need this function to get an object of 'SparseGRM'. 
#' Our simulations show that using sparse GRM is almost the same reliable as using dense GRM, with a huge improvement of computational efficiency. 
#' \cr\cr
#' Step 1: Download GCTA software from https://cnsgenomics.com/software/gcta/#Overview; 
#' \cr
#' Step 2: Run getSparseGRMParallel(); 
#' \cr
#' Step 3: Run getSparseGRM().
#' @param outPrefix a path to specify the prefix of output files  
#' @param chrVec a character vector to specify all chromosomes (e.g. 1:22), better autosomes.  
#' @param PlinkFile a path to Plink files. The current version (gcta_1.93.1beta) of gcta software does not support difference prefix names for bim, bed and fam files. 
#' @param relatednessCutoff a cutoff for sparse GRM, if the kinship coefficient is greater than this cutoff, it will be in sparse GRM.
#' @param minMafGRM Minimal value of MAF cutoff to select markers (from Plink files) to construct GRM.
#' @param maxMissingGRM Maximal value of missing rate to select markers (from Plink files) to construct GRM.
#' @param nPartsGRM GCTA software can split subjects to multiple parts. For UK-Biobank analysis, it is recommanded to use 250 parts. 
#' @examples 
#' ## First download GCTA software, then you can use the following commands
#' famFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.fam", package = "POLMM")
#' PlinkFile = gsub(".fam", "", famFile)   # fam/bim/bed files should have the same prefix
#' gcta64File = "/net/snowwhite/home/wenjianb/gcta64"  # should download gcta first (https://cnsgenomics.com/software/gcta/#Download, author version is gcta_1.93.1beta)
#' outPrefix = gsub("nSNPs-10000-nsubj-1000-ext.fam","SparseGRM/sp", famFile) # Dir/prefix of the output files (in step 2)
#' SparseGRMFile = gsub("nSNPs-10000-nsubj-1000-ext.fam","SparseGRM.RData", famFile) # output file (in step 3)
#' nPartsGRM = 2;   # nPartsGRM = 250 for UK Biobank data analysis
#' chrVec = 1:22    # maybe paste0("chr",1:22), depending on the plink bim file
#' 
#' ## Step 2:
#' ## We strongly suggest parallel computing for different pairs of (chrParallel, partParallel). 
#' for(chrParallel in chrVec){
#'   for(partParallel in 1:nPartsGRM){
#'     getSparseGRMParallel(outPrefix, chrParallel, partParallel, gcta64File, PlinkFile, nPartsGRM = nPartsGRM)
#'   }
#' }
#' 
#' ## After that, in "outPrefix", there will be results (needs large amount of storage) corresponding to different pairs of (chrParallel, partParallel).
#' 
#' ## Step 3:
#' ## Combine results in step 2 to calculate an object with class of SparseGRM for POLMM_Null_Model(),
#' SparseGRM = getSparseGRM(outPrefix, chrVec, PlinkFile, nPartsGRM = nPartsGRM)
#' save(SparseGRM, file=SparseGRMFile)
#' 
#' ## NOTE: You can change some options such as (minMafGRM, maxMissingGRM, nPartsGRM), but keep in mind that functions getSparseGRMParallel() and getSparseGRM() should use the same change of these options.
getSparseGRM = function(outPrefix,
                        chrVec,
                        PlinkFile,
                        relatednessCutoff = 0.05,
                        minMafGRM = 0.01,
                        maxMissingGRM = 0.1,
                        nPartsGRM = 250)
{
  PlinkName = basename(PlinkFile)
  chrVec = as.character(chrVec)
  nDigits = floor(log10(nPartsGRM)) + 1 # 1-9 -> 1; 10-99 -> 2; 100:999 -> 3.
  
  AllIDs = c()
  n0 = 0
  SparseGRM = list()
  
  ## cycle for nPartsGRM
  for(i in 1:nPartsGRM){
    print(paste("Analyzing part",i,"of total",nPartsGRM,"parts."))
    tempList = list()
    for(chr in chrVec){
      # print(paste("Analyzing chromosome",chr))
      outFile = paste0(outPrefix, "-Plink-",PlinkName,"-chr-",chr,"-minMaf-",minMafGRM,"-maxMissing-",maxMissingGRM)
      ## Three files generated by GCTA
      IDFile = paste0(outFile,".part_",nPartsGRM,"_",formatC(i, width = nDigits, flag = '0'),".grm.id")
      BinFile = paste0(outFile,".part_",nPartsGRM,"_",formatC(i, width = nDigits, flag = '0'),".grm.bin")
      NFile = paste0(outFile,".part_",nPartsGRM,"_",formatC(i, width = nDigits, flag = '0'),".grm.N.bin")
      ## read in the three files
      IDs = read.table(IDFile, stringsAsFactors=F)
      ID = IDs$V2
      n1 = n0 + length(ID)
      nData = (n1 - n0) * (n0 + n1 + 1) / 2
      
      grm = readBin(BinFile, n = nData, what = numeric(0), size = 4)
      nMarkers = readBin(NFile, n = nData, what = numeric(0), size = 4)
      tempList[[chr]] = list(ID = ID, grm = grm, nMarkers = nMarkers)
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





