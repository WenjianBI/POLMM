# POLMM
Proportional Odds Logistic Mixed Model (POLMM) for ordinal categorical data analysis

### How to install and load this package

```{r}      
library(devtools)  # author version: 2.3.0
install_github("WenjianBi/POLMM")
library(POLMM)
?POLMM  # manual of POLMM() function with an example code
```
Current version is 0.2.2. For older version and version update information, plesase refer to OldVersions/

Please do not hesitate to contact me (wenjianb@umich.edu) if you meet any problem. Suggestions or comments are also welcome.

### We support Dense GRM and Sparse GRM to adjust for sample relatedness

**Using dense GRM:** We use similar strategies as in SAIGE and BOLT-LMM. When sample size is greater than 400K, this is still a little bit slow.  

**Using sparse GRM (Our recommendation):** Much faster than using dense GRM. If users want to use sparse GRM (which supports LOCO option), should pass argument 'SparseGRM' to main function POLMM_Null_Model().

**How to get an R object of SparseGRM:**  
1. Download gcta software from https://cnsgenomics.com/software/gcta/#Overview
1. Use function getSparseGRMParallel() to generate GRM files for each chromosome (split all subjects into multiple parts: we use 250 parts for UK Biobank analysis)
1. Use function getSparseGRM() to combine all GRM files to generate an R object of "SparseGRM" to be passed to main function POLMM_Null_Model() 

