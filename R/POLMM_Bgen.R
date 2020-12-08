#' Test for association between genotype and an ordinal categorical variable via Proportional Odds Logistic Mixed Model (POLMM)
#' 
#' Test for association between genotype and an ordinal categorical variable via Proportional Odds Logistic Mixed Model (POLMM)
#' 
#' @param objNull output object of the POLMM_Null_Model() function 
#' @param plink.file character, represents the prefix of PLINK input file.
#' @param output.file character, represents the prefix of output file.
#' @param memory.chunk a numeric value (default: 4, unit=Gb) to specify how much memory is used to store genotype matrix from plink files.
#' @param chrVec a character or a character vector to specify chromosome(s) of the markers in Geno.mtx. Must be specified unless LOCO = F.
#' @param SPAcutoff a standard deviation cutoff (default=2). If the standardized test statistic < SPAcutoff, normal approximation is used, otherwise, saddlepoint approximation is used.
#' @param minMAF a cutoff of the minimal minor allele frequencies (MAFs). Any markers with MAF < minMAF will be excluded from the analysis.
#' @param maxMissing a cutoff of the maximal missing rate. Any markers with missing rate > maxMissing will be excluded from the analysis.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#'                      "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e. 2p where p is MAF).
#' @param G.model a character string (default: "Add") to specify additive ("Add"), dominant ("Dom"), or recessive ("Rec") model. 
#'                If "Dom", GVec = ifelse(GVec >= 1, 1, 0), if "Rec", GVec = ifelse(GVec <= 1, 0, 1). Be very careful if the gneotyp is imputed data.
#' @return results can be found in 'output.file': a matrix with the following elements
#' \item{ID}{Marker IDs from bim file}
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
#' outFile = gsub("-ext.fam","-ext.output",famFile)
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
#' SparseGRMFile = system.file("SparseGRM", "SparseGRM.RData", package = "POLMM")
#' load(SparseGRMFile)   ## check getSparseGRM() for more details about how to make an R object of "SparseGRM" using Plink files. 
#' objNull = POLMM_Null_Model(as.factor(outcome)~Cova1+Cova2, 
#'                            SparseGRM = SparseGRM,
#'                            data=egData, PlinkFile = PlinkFile, subjData = egData$IID)
#'                            
#' ## If control$seed is not changed, objNull$tau should be 0.8506
#' objNull$tau
#' 
#' ## when using function POLMM()/POLMM.plink(), argument chrVec/chrVec.plink should be from
#' names(objNull$LOCOList)
#' 
#' POLMM.plink(objNull, PlinkFile, outFile, 
#'             chrVec.plink = c(1,2,3,4));  # if you only want to analyze variants in chromosomes 1-4 
#' outPOLMM = read.table(outFile, header = T)
#' head(outPOLMM)
#' @export
#' @import data.table
#' @import tidyr
#' @import dbplyr
#' @import RSQLite
POLMM.bgen = function(objNull,
                      bgen.file,            # plink prefix
                      bgi.file,
                      sample.file,
                      output.file,
                      chrVec.bgen = NULL,
                      memory.chunk = 4,
                      SPAcutoff = 2,
                      minMAF = 0.0001,
                      maxMissing = 0.15,
                      impute.method = "fixed",
                      G.model = "Add")
{
  ## check if SAIGE has been installed
  
  SI = sessionInfo()
  if(!is.element("SAIGE",names(SI$otherPkgs)))
    stop("We use function in SAIGE package to read in bgen files. Please install SAIGE package (0.36.3.3: only supports linux OS) first. It works for bgen-1.3 with 8 bits.")
  
  ## check bgen files input
  
  if(!file.exists(bgen.file)) stop("Could not find bgen.file")
  if(!file.exists(bgi.file)) stop("Could not find bgi.file")
  if(!file.exists(sample.file)) stop("Could not find sample.file")
  
  if(file.exists(output.file)) 
    stop("'output.file' existed. Please give a different 'output.file' or remove the existing 'output.file'.")
  
  subjIDs_Null = objNull$subjIDs
  
  samples = data.table::fread(sample.file)
  gIDs = samples$ID_2
  
  subjIndex_Null = match(subjIDs_Null, gIDs, nomatch = 0)
  if(any(subjIndex_Null == 0)) stop("All subjects in null model fitting should be also in plink files.")
  
  RangesList = split.ranges.bgen(bgi.file, length(subjIDs_Null), memory.chunk)
  
  M = RangesList$M
  RSID.set = RangesList$RSID.set
  
  print(paste0("Totally ", M, " markers in plink files."))

  if(!is.null(chrVec.bgen))
    Ranges = subset(Ranges, is.element(chr, chrVec.bgen))
  
  # n.chunk = nrow(Ranges)
  n.chunk = length(RSID.set)
  
  print(paste0("Split all markers into ", n.chunk, " chunks."))
  # print(paste0("Each chunk includes less than ", M.chunk, " markers."))
  
  for(i in 1:n.chunk){
    print(paste0("Analyzing chunk ",i,"/",n.chunk,"."))
    # chr = Ranges$chr[i]
    # pos.start = Ranges$pos.start[i]
    # pos.end = Ranges$pos.end[i]
    rsid = RSID.set[[i]]
    
    Geno.mtx = SAIGE.read.bgen.range(bgen.file, 
                                     bgi.file, 
                                     rsid,
                                     # chr,
                                     # pos.start,
                                     # pos.end,
                                     subjIDs_Null,
                                     gIDs)
    
    if(names(objNull$LOCOList)[1] == "LOCO=F"){
      chrVec = "LOCO=F"
    }else{
      chrVec = chr
    }
    
    output = POLMM(objNull, Geno.mtx, chrVec, 
                   SPAcutoff, minMAF, maxMissing, impute.method, G.model,
                   G.missing = -9)
    
    if(i == 1){
      data.table::fwrite(output, output.file, sep = "\t", append = F, row.names = F, col.names = T,
                         na = "NA")
    }else{
      data.table::fwrite(output, output.file, sep = "\t", append = T, row.names = F, col.names = F,
                         na = "NA")
    }
    
  }
  
  print("Analysis Complete.")
  print(Sys.time())
  
  # return(output)  # output is stored in outFile
  
}

split.ranges.bgen = function(bgi.file, nSubj, memory.chunk)
{
  db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgi.file)
  on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
  
  infos = dplyr::tbl(db_con, "Variant") %>% dplyr::select(chromosome, position, rsid) %>% dplyr::arrange(chromosome, position) %>% dplyr::collect(rsid)
  
  M = nrow(infos)
  
  uchr = unique(infos$chromosome)
  
  nEachRange = floor(memory.chunk * 1e9 / nSubj / 4) # "double" needs 4 bytes
  
  # Ranges = c()
  RSID.set = list()
  idx = 1
  for(chr in uchr){
    infos1 = infos %>% dplyr::filter(chromosome == chr)
    nchr = nrow(infos1)
    idx.start = seq(1, nchr, nEachRange)
    idx.end = c(idx.start[-1]-1, nchr)
    # pos.start = infos1$position[idx.start]
    # pos.end = infos1$position[idx.end]
    # range = paste0(chr, ":", pos.start, "-", pos.end)
    # range = c(chr, pos.start, pos.end)
    # Ranges = rbind(Ranges, range)
    for(j in 1:length(idx.start)){
      RSID.set[[idx]] = unique(infos1$rsid[idx.start[j]:idx.end[j]])
      idx = idx + 1
    }
  }
  # Ranges = as.data.frame(Ranges, stringsAsFactors = F)
  # colnames(Ranges) = c("chr", "pos.start", "pos.end")
  
  # return(list(Ranges = Ranges, M = M))
  return(list(RSID.set = RSID.set, M = M))
}


SAIGE.read.bgen.range = function(bgen.file, 
                                 bgi.file, 
                                 rsid,
                                 # chr,
                                 # pos.start,
                                 # pos.end
                                 subjIDs_Null,
                                 gIDs)
{
  sampleIndex = match(gIDs, subjIDs_Null, nomatch = NA)
  N = sum(!is.na(sampleIndex))
  
  ##########
  
  bgenFile = bgen.file
  bgenFileIndex = bgi.file
  
  # ranges_to_include = data.frame(chromosome = chr, 
  #                                start = as.numeric(pos.start), 
  #                                end = as.numeric(pos.end))
  ranges_to_include = ranges_to_exclude = data.frame(chromosome = NULL, start = NULL, end = NULL)
  ids_to_exclude = as.character(vector())
  ids_to_include = rsid
  
  Mtest = SAIGE:::setgenoTest_bgenDosage(bgenFile,
                                         bgenFileIndex,
                                         ranges_to_exclude = ranges_to_exclude,
                                         ranges_to_include = ranges_to_include,
                                         ids_to_exclude = ids_to_exclude,
                                         ids_to_include = ids_to_include)
  
  SAIGE:::SetSampleIdx(sampleIndex-1, N)
  
  Geno.mtx = matrix(0, N, Mtest)
  SNPID = rep(0,Mtest)
  for(i in 1:Mtest){
    Gx = SAIGE:::getDosage_bgen_withquery()
    markerInfo = SAIGE:::getMarkerInfo()
    Geno.mtx[,i] = Gx$dosage
    SNPID[i] = Gx$variants$rsid
  }
  
  rownames(Geno.mtx) = subjIDs_Null
  colnames(Geno.mtx) = SNPID
  return(Geno.mtx)
}


