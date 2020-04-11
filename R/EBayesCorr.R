EBayesCorr=function(NLExprMod){
  NLExprMelt=melt(NLExprMod, measure.vars=c("RefBase"))
  altvals=c(NLExprMelt$DNAalt, NLExprMelt$Alt)
  refvals=c(NLExprMelt$DNAalt+NLExprMelt$DNAref, NLExprMelt$Ref+NLExprMelt$Alt)
  rmx=which(is.na(refvals)==T)
  
  if (length(rmx)>0){
    altvals=altvals[-rmx]
    refvals=refvals[-rmx]
  }
  
  ll <- function(alpha, beta) {
    -sum(dbetabinom.ab(altvals, refvals, alpha, beta, log = TRUE))
  }
  m <- stats4::mle(ll, start = list(alpha = 1, beta = 10), method = "L-BFGS-B")
  #coef(m)
  
  for (i in 1:length(NLExprMod)){  
    NLExprMod[[i]]$DNAEBProb=(NLExprMod[[i]]$DNAalt+coef(m)[1])/(coef(m)[1]+coef(m)[2]+NLExprMod[[i]]$DNAalt+NLExprMod[[i]]$DNAref)
    NLExprMod[[i]]$RNAEBProb=(NLExprMod[[i]]$Alt+coef(m)[1])/(coef(m)[1]+coef(m)[2]+NLExprMod[[i]]$Alt+NLExprMod[[i]]$Ref)
  }
  NLExprMod
}
