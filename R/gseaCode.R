
# Wrapper to run GSEA

gseaCode=function(edrtt, hits=NULL, pway, overlap=F, GSEA=T, overlapCutoff=2){
 
  if (is.matrix(edrtt)==1){ 
    t1=as.vector(edrtt[ ,1])
    names(t1)=rownames(edrtt)
    edrtt=t1
  }
  
  ListGSC=list(pway=pway)
  if (overlap==T & is.null(hits)==T){
    hits=names(edrtt)[which(abs(edrtt)>=overlapCutoff)]
    } else if (overlap==T & !is.null(hits)) {
    gsca <- HTSanalyzeR2::GSCA(listOfGeneSetCollections=ListGSC, 
                 geneList=edrtt, hits = hits)
  } else {
    gsca <- HTSanalyzeR2::GSCA(listOfGeneSetCollections=ListGSC, 
                 geneList=edrtt)
  }
  
  gsca1 <- HTSanalyzeR2::preprocess(gsca, species="Hs", initialIDs="SYMBOL",
                      keepMultipleMappings=TRUE, duplicateRemoverMethod="max",
                      orderAbsValue=FALSE)
  
  gsca2 <- HTSanalyzeR2::analyze(gsca1, para=list(pValueCutoff=0.05, pAdjustMethod="BH",
                                    nPermutations=100, minGeneSetSize=5,
                                    exponent=1), doGSOA = overlap, doGSEA=GSEA)
  #gsca2<-appendGSTerms(gsca2, pway)
  gsca2
}
