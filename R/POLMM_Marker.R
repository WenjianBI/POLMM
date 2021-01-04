
#' Proportional Odds Logistic Mixed Model (POLMM)
#' 
#' Single-marker analysis: Test for association between genetic variants and an ordinal categorical phenotype via POLMM
#' 
#' @param objNull an output object of the POLMM_Null_Model() function with a class of "POLMM_NULL_Model". 
#' @param GenoFile a character of genotype file. Three types of genotype files are supported: PLINK ("prefix.bed"), BGEN ("prefix.bgen"), and VCF ("prefix.vcf" or "prefix.vcf.gz"). 
#' @param GenoFileIndex additional index files corresponding to the "GenoFile". The default is NULL, that is, to share the same prefix as GenoFile. PLINK: c("prefix.bim", "prefix.fam"), BGEN: c("prefix.bgi"), and VCF: c("prefix.vcf.tbi") or c("prefix.vcf.gz.tbi").
#' @param OutputFile a character of output file to store the analysis results
#' @param chrom a character to specify chromosome of the markers in analysis. Must be specified unless LOCO = F when fitting the null model.
#' @param POLMM.control a list of parameters for controlling the POLMM.Marker(). The default is NULL, that is, to use the default parameters. Check 'Details' and 'Examples' for more details.
#' @details 
#' More information about the list of 'POLMM.control'
#' \itemize{
#' \item{impute_method: imputation method when genotype data is missing [default="fixed"].}
#' \item{missing_cutoff: exclude markers with missing rate greater than this cutoff [default=0.15].}
#' \item{min_maf_marker: exclude markers with MAF less than this cutoff [default=0.001].}
#' \item{min_mac_marker: exclude markers with MAC less than this cutoff [default=20].}
#' \item{nMarkers_output: output results after analyzing each chunk of fixed number of markers [default=10000].}
#' \item{SPA_cutoff: a cutoff to determine "normal distribution approximation" or "saddlepoint approximation" [default=2].}
#' }
#' @return The analysis results will be written to OutputFile with the following columns
#' \item{Marker}{Marker IDs extracted from "GenoFile" or "GenoFileIndex".}
#' \item{Info}{Marker Infomation of "CHR:POS:REF:ALT". This information is from "GenoFile" or "GenoFileIndex" and does not change even if the REF/ALT alleles are flipped in analysis}
#' \item{Freq}{Minor allele frequency (always < 0.5) in analysis.}
#' \item{Flip}{a logical value indicating if the REF/ALT alleles were switched in analysis. This information is useful to estimate the effect direction.}
#' \item{Stat}{Score statistics (changed to beta later.) The direction depends on both "Info" and "Flip".}
#' \item{Var}{Estimated variance of the score statistics (changed to se.beta later)}
#' \item{Pvalue}{p-value from normal distribution approximation or saddlepoint approximation.}
#' @examples
#' famFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.fam", package = "POLMM")
#' GenoFile = gsub("-ext.fam", "-ext.bed", famFile)
#' AnnoFile = system.file("extdata", "AnnoFile.txt", package = "POLMM")
#' OutputFile = gsub("AnnoFile","OutputFile",AnnoFile)
#' SparseGRMFile = system.file("SparseGRM", "SparseGRM.RData", package = "POLMM")
#' load(SparseGRMFile)
#' objNullFile = system.file("objNull.RData", package = "POLMM")
#' load(objNullFile)
#' chrom = 1
#' 
#' OUTPUT = POLMM.Marker(objNull, AnnoFile, GenoFile, GenoFileIndex = NULL, OutputFile, 
#'                       SparseGRM, chrom, POLMM.control = list(max_maf_region = 0.5))
#'      
#' @export
#' @import data.table

POLMM.Marker = function(objNull,
                        GenoFile,
                        GenoFileIndex = NULL,
                        OutputFile,
                        chrom,
                        POLMM.control = NULL)
{
  if(class(objNull) != "POLMM_NULL_Model")
    stop("class(objNull) should be 'POLMM_NULL_Model'.")
  
  if(file.exists(OutputFile))
    stop(paste0("'OutputFile' of '",OutputFile,"' has existed. Please use another 'OutputFile' or remove the existing 'OutputFile'."))
  
  SubjID.step1 = as.character(objNull$subjIDs);
  
  # check the setting of SKAT.control, if not specified, the default setting will be used
  POLMM.control = check.POLMM.control(POLMM.control, "step2")
  print(POLMM.control)
  genoList = setGenoInput(GenoFile, GenoFileIndex, SubjID.step1, output.marker = T)
  
  if(objNull$controlList$LOCO){
    if(!chrom %in% names(objNull$LOCOList))
      stop("'chrom' should be in names(objNull$LOCOList).")
    obj.CHR = objNull$LOCOList[[chrom]]
  }else{
    # to be continued
  }
  
  # single marker analysis does not require sparse GRM any more 
  # Note: it might be not so accurate if min_mac_marker is very low
  SPmatR.CHR = list(locations = c(0,0), values = 1)
  
  setPOLMMobjInR(obj.CHR$muMat,
                 obj.CHR$iRMat,
                 objNull$Cova,
                 objNull$yVec,          # 1 to J
                 SPmatR.CHR,
                 objNull$tau,
                 POLMM.control$printPCGInfo,
                 POLMM.control$tolPCG,
                 POLMM.control$maxiterPCG)
  
  nMarkers_output = POLMM.control$nMarkers_output;
  print(paste0("The current POLMM.control$nMarkers_output is ", nMarkers_output,"."))
  
  all.markers = genoList$marker
  marker.list = split.marker(all.markers, nMarkers_output)
  print(paste0("We split all markers to ", length(marker.list), " blocks, each of which includes no greater than ", nMarkers_output, " markers."))
  
  for(i in 1:length(marker.list)){
    print(paste0("Analyzing Block ", i, "/", length(marker.list)," ......"))
    markers = marker.list[[i]]
    OutList = MAIN_Marker(markers,
                          POLMM.control$SPA_cutoff,
                          POLMM.control$missing_cutoff,
                          POLMM.control$min_maf_region,
                          POLMM.control$min_mac_region,
                          obj.CHR$varRatio)
    
    StatVec = OutList$StatVec
    VarSVec = OutList$VarSVec
    adjPVec = OutList$pvalVec;
    markerVec = OutList$markerVec
    freqVec = OutList$freqVec
    flipVec = OutList$flipVec
    infoVec = OutList$infoVec       # marker infomation: CHR:POS:REF:ALT
    
    OUT.Marker = data.frame(Marker = markerVec,
                            Info = infoVec,
                            Freq = freqVec,
                            Flip = flipVec,
                            Stat = StatVec,
                            Var = VarSVec,
                            Pvalue = adjPVec)
    
    if(i == 1){
      data.table::fwrite(OUT.Marker, OutputFile, quote=F, sep="\t", append = F, col.names = T)
    }else{
      data.table::fwrite(OUT.Marker, OutputFile, quote=F, sep="\t", append = T, col.names = F)
    }
  }
  
  message = paste0("The analysis results has been saved to '", OutputFile,"'.")
  return(message)
}


