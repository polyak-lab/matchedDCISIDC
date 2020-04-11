##################### 
##AppendDNAMutMod
######################
# 2. Append DNA Mut
#     input list of predicted neoantigens and mutational data
# AppendDNAMut


AppendDNAMut=function(NLExpr, MutList){
  NLExprMod=NLExpr
  DupNeo=lapply(NLExpr, function(x) which(duplicated(x[ ,1:2])))
  
  for (i in 1:length(NLExpr)){
    if (length(DupNeo[[i]])>0){
      for (j in DupNeo[[i]]){
        Sx=NLExpr[[i]][ c(j-1, j), ]  
        mx=which.max(Sx$Alt)
        lx=Sx[-mx, ]
        pvals=sapply(1:nrow(lx), function(x) paste(lx[x, ], collapse= " "))
        PList=sapply(1:nrow(NLExprMod[[i]]), function(x) paste(NLExprMod[[i]][x, ], collapse= " "))
        midx=match(pvals, PList)
        NLExprMod[[i]]=NLExprMod[[i]][-midx, ]
      }
    }
  }
  lSample=sapply(names(NLExprMod), function(x) match(x, names(MutList)))
  for (i in 1:length(NLExprMod)){
    rnadf=MutList[[lSample[i]]]
    #[which(MutList[[lSample[i]]]$Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins")), ]
    rnadf$Prob=rnadf$t_alt_count/(rnadf$t_alt_count+rnadf$t_ref_count) ## double 
    PcritDepth=qnbinom(0.95, 2, rnadf$Prob)
    #rnadf$Depth2=PcritDepth
    rnadfM=PcritDepth[match(paste(NLExprMod[[i]]$chr, NLExprMod[[i]]$Start.bp), paste(rnadf$Chromosome, rnadf$Start_position))]
    dnaprob=rnadf$Prob[match(paste(NLExprMod[[i]]$chr, NLExprMod[[i]]$Start.bp), paste(rnadf$Chromosome, rnadf$Start_position))]
    SSum=as.numeric(as.character(NLExprMod[[i]]$Ref))+as.numeric(as.character(NLExprMod[[i]]$Alt))
    Gene=rnadf$Hugo[match(paste(NLExprMod[[i]]$chr, NLExprMod[[i]]$Start.bp), paste(rnadf$Chromosome, rnadf$Start_position))]
    DNARef=rnadf$t_ref_count[match(paste(NLExprMod[[i]]$chr, NLExprMod[[i]]$Start.bp), paste(rnadf$Chromosome, rnadf$Start_position))]
    DNAAlt=rnadf$t_alt_count[match(paste(NLExprMod[[i]]$chr, NLExprMod[[i]]$Start.bp), paste(rnadf$Chromosome, rnadf$Start_position))]
    dbSNP=rnadf$dbSNP_RS[match(paste(NLExprMod[[i]]$chr, NLExprMod[[i]]$Start.bp), paste(rnadf$Chromosome, rnadf$Start_position))]
    NLExprMod[[i]]=cbind(NLExprMod[[i]],Gene=Gene, dnaprob=dnaprob, PcritDepth=rnadfM, SSum=SSum, DNAref=DNARef, DNAalt=DNAAlt,dbSNP=dbSNP)
  }
  NLExprMod 
}

AppendDNAMutMod=function(NLExpr, MutList){
  NLExprMod=NLExpr
  DupNeo=lapply(NLExpr, function(x) which(duplicated(x[ ,1:2])))
  
  for (i in 1:length(NLExpr)){
    if (length(DupNeo[[i]])>0){
      for (j in DupNeo[[i]]){
        Sx=NLExpr[[i]][ c(j-1, j), ]  
        mx=which.max(Sx$Alt)
        lx=Sx[-mx, ]
        pvals=sapply(1:nrow(lx), function(x) paste(lx[x, ], collapse= " "))
        PList=sapply(1:nrow(NLExprMod[[i]]), function(x) paste(NLExprMod[[i]][x, ], collapse= " "))
        midx=match(pvals, PList)
        NLExprMod[[i]]=NLExprMod[[i]][-midx, ]
      }
    }
  }
  lSample=sapply(names(NLExprMod), function(x) match(x, names(MutList)))
  
  for (i in 1:length(NLExprMod)){
    rnadf=MutList[[lSample[i]]]
    #[which(MutList[[lSample[i]]]$Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins")), ]
    rnadf$Prob=rnadf$Alternate.reads/(rnadf$Total.reads) ## double 
    PcritDepth=qnbinom(0.95, 2, rnadf$Prob)
    #rnadf$Depth2=PcritDepth
    NLEidx=paste(NLExprMod[[i]]$chr, NLExprMod[[i]]$Start.bp)
    RFidx=paste(substr(rnadf$Chrom, 4, 5), rnadf$Position.2)         
    rnadfM=PcritDepth[match(NLEidx, RFidx)]
    dnaprob=rnadf$Prob[match(NLEidx, RFidx)]
    SSum=as.numeric(as.character(NLExprMod[[i]]$Ref))+as.numeric(as.character(NLExprMod[[i]]$Alt))
    Gene=rnadf$Hugo[match(NLEidx, RFidx)]
    DNARef=rnadf$Total.reads[match(NLEidx, RFidx)]-rnadf$Alternate.reads[match(NLEidx, RFidx)]
    DNAAlt=rnadf$Alternate.reads[match(NLEidx, RFidx)]
    dbSNP=rnadf$rsID[match(NLEidx, RFidx)]
    NLExprMod[[i]]=cbind(NLExprMod[[i]],Gene=Gene, dnaprob=dnaprob, PcritDepth=rnadfM, SSum=SSum, DNAref=DNARef, DNAalt=DNAAlt,dbSNP=dbSNP)
  }
  
  NLExprMod 
}
