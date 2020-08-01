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

**Using dense GRM:** Similar strategies as in SAIGE and BOLT-LMM were used when fitting the null model.  

**Using sparse GRM (Our recommendation):** Much faster than using dense GRM and also supports LOCO option. Users should pass an R object of 'SparseGRM' to the main function POLMM_Null_Model(). The below is the manual to make an R object of 'SparseGRM' which can include multiple GRMs for different chromosomes. For one dataset, only need to make 'SparseGRM' for once.  

**How to make an R object of SparseGRM:**  
1. Use function getSparseGRMParallel() to generate GRM files for each chromosome (can split all subjects into multiple parts: we use 250 parts for UK Biobank analysis)
1. Use function getSparseGRM() to combine all GRM files to generate an R object of "SparseGRM" to be passed to main function POLMM_Null_Model() 

