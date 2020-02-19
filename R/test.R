
## one SNP version after removing all unnecessary sections

POLMM.one.SNP = function(GVec,            # n vector
                         objP,
                         RPsiR,
                         SPAcutoff=2)
{
  adjG = outputadjGFast_v1(GVec, objP, RPsiR)
  adjGMat = adjG$adjGMat
  Stat = adjG$Stat
  VarW = adjG$VarW
  VarP = VarW * r
  z = abs(Stat)/sqrt(VarP)
  pval.spa = pval.norm = 2*pnorm(-1*z, lower.tail=T)
  if(z > SPAcutoff){   
    pval.spa = Saddle_Prob(Stat, VarP, VarW,
                           adjGMat, objP[["muMat"]], objP[["iRMat"]])   # n x (J-1)
  }
  
  OutVec = c(Stat, VarW, VarP, pval.norm, pval.spa)
  return(OutVec)
}



