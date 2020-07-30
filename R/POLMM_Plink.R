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
#' POLMM.plink(objNull, PlinkFile, outFile, 
#'             chrVec.plink = c(1,2,3,4));  # since sparseGRM only includes these 4 chromosomes
#' outPOLMM = read.table(outFile, header = T)
#' head(outPOLMM)
#' @export
#' @import seqminer

POLMM.plink = function(objNull,
                       plink.file,            # plink prefix
                       output.file,
                       chrVec.plink,
                       memory.chunk = 4,
                       SPAcutoff = 2,
                       minMAF = 0.0001,
                       maxMissing = 0.15,
                       impute.method = "fixed",
                       G.model = "Add")
{
  ## check plink files input
  bim.file = paste0(plink.file, ".bim")
  bed.file = paste0(plink.file, ".bed")
  fam.file = paste0(plink.file, ".fam")
  
  if(!file.exists(bim.file)) stop("Could not find paste0(plink.file,'.bim')")
  if(!file.exists(bed.file)) stop("Could not find paste0(plink.file,'.bed')")
  if(!file.exists(fam.file)) stop("Could not find paste0(plink.file,'.fam')")
  if(file.exists(output.file)) stop("'output.file' existed. Please give a different 'output.file' or remove the existing 'output.file'.")
  
  fam.data = read.table(fam.file, stringsAsFactors = F)
  bim.data = read.table(bim.file, stringsAsFactors = F)
  
  subjIDs_Null = objNull$subjIDs
  gIDs = fam.data$V2
  chrVecTot = bim.data$V1
  
  subjIndex_Null = match(subjIDs_Null, gIDs, nomatch = 0)
  if(any(subjIndex_Null == 0)) stop("All subjects in null model fitting should be also in plink files.")
  
  N = length(subjIndex_Null)
  M = nrow(bim.data)
  
  if(missing(chrVec.plink)){
    pos = 1:M
  }else{
    pos = which(is.element(chrVecTot, chrVec.plink))
  }
  M1 = length(pos)
  
  print(paste0("Totally ", M, " markers in plink files."))
  print(paste0("After filtering by 'chrVec.plink', ", M1, " markers are left."))
  
  M.chunk = floor(memory.chunk * 1e9 / 4 / N)
  n.chunk = ceiling(M1 / M.chunk)
  
  print(paste0("Split all markers into ", n.chunk, " chunks."))
  print(paste0("Each chunk includes less than ", M.chunk, " markers."))
  
  for(i in 1:n.chunk){
    if(i == n.chunk){
      markerIndex = ((n.chunk-1)*M.chunk+1):M1;
    }else{
      markerIndex = 1:M.chunk + (i-1) * M.chunk;
    }
    print(paste0("Analyzing chunk ",i,"/",n.chunk,"."))
    markerIndex = pos[markerIndex]
    
    Geno.mtx = seqminer::readPlinkToMatrixByIndex(plink.file, subjIndex_Null, markerIndex)
    colnames(Geno.mtx) = bim.data$V2[markerIndex]
    chrVec = chrVecTot[markerIndex]
    
    output = POLMM(objNull, Geno.mtx, chrVec, 
                   SPAcutoff, minMAF, maxMissing, impute.method, G.model)
    
    if(i == 1){
      data.table::fwrite(output, output.file, sep = "\t", append = F, row.names = F, col.names = T)
    }else{
      data.table::fwrite(output, output.file, sep = "\t", append = T, row.names = F, col.names = F)
    }
  }
  
  print("Analysis Complete.")
  print(Sys.time())
  
  # return(output)  # output is stored in outFile
  
}