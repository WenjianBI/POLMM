# POLMM
Proportional Odds Logistic Mixed Model (POLMM) for ordinal categorical data analysis

### Software dependencies and operating systems
The package has been tested under linux and windows systems. 

To make a sparse GRM, the package includes codes of GCTA software (version 1.93.1beta, https://cnsgenomics.com/software/gcta/#fastGWA); to read in bgen data input, the package requires installing SAIGE package (version 0.36.3, https://github.com/weizhouUMICH/SAIGE). These two functions are only supported in linux system since the software dependencies are only supported in linux system.

License: GPL (>= 3)

### How to install and load this package

```{r}      
library(devtools)  # author version: 2.3.0
install_github("WenjianBi/POLMM")
library(POLMM)
?POLMM  # manual of POLMM() function with an example code, expected output and expected run time for demo
```
Current version is 0.2.3. The package installation typically requires < 3 minutes on a normal desktop computer. 

Please do not hesitate to contact me (wenjianb@pku.edu.cn) if you meet any problem. Suggestions or comments are also welcome.

### We support Dense GRM and Sparse GRM to adjust for sample relatedness

**Using dense GRM:** Similar strategies as in SAIGE and BOLT-LMM were used when fitting the null model.  

**Using sparse GRM (Our recommendation):** Much faster than using dense GRM and also supports LOCO option. Users should pass an R object of 'SparseGRM' to the main function POLMM_Null_Model(). The below is the manual to make an R object of 'SparseGRM' which can include multiple GRMs for different chromosomes. For one dataset, only need to make 'SparseGRM' for once.  

**How to make an R object of SparseGRM:**  
1. Use function getSparseGRMParallel() to generate GRM files for each chromosome (can split all subjects into multiple parts: we use 250 parts for UK Biobank analysis)
1. Use function getSparseGRM() to combine all GRM files to generate an R object of "SparseGRM" to be passed to main function POLMM_Null_Model() 

### About the effect direction
For POLMM() function, the input is a genotype matrix. For each marker with genotype vector of GVec, we first calculate the allele frequency (AF) of mean(GVec)/2, if the AF > 0.5, then we use GVec = 2 - GVec instead (to reduce computational time) and set switch.alllel = TRUE. Hence, if switch.allele = TRUE, the effect direction is changed.

For POLMM.plink() function, we first use seqminer::readPlinkToMatrixByIndex() function to get a genotype matrix from PLINK files. In seqminer (version 8.0), this function follows "ref-first", that is, A1 is usually major(REF) and A2 is usually minor(ALT). For each marker, genotype matrix is the counts of A2 allele and beta corresponds to A2 allele. If switch.allele = TRUE, then beta corresponds to A1 allele.  


### PheWeb for UK Biobank data analysis results

We have applied this method to analyze 258 categorical phenotypes in UK Biobank data of 408,961 samples from white British participants with European ancestry. The results can be found in http://polmm.leelabsg.org/.
