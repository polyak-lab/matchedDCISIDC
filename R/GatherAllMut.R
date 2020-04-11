## gather all mut in a bed file for phylowgs
GatherAllMut=function(samp1, samp2, bedName="untitled"){
  ListA=samp1[ ,c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Variant_Classification")]  
  ListB=samp2[ ,c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Variant_Classification")]  
  ListAll=rbind(ListA, ListB)
  ListAll=ListAll[which(ListAll$Variant_Classification%in%c("Silent", "IGR", "In_Frame_Del", "In_Frame_Ins", "Intron", "3'UTR", "5'Flank", "De_novo_Start_InFrame")), ]
  #ListAll$ID=paste(ListAll$Hugo_Symbol, ListAll$Start_position, sep="_")
  rmx=which(duplicated(ListAll))
  
  if (length(rmx)>0){
    DupList=ListAll[-which(duplicated(ListAll)), ]
  }else{
    DupList=ListAll
  }
  # write the bed file
  DupList2=DupList[ ,c(2:4, 1)]
  DupList2$Start_position=DupList2$Start_position-1
  write.table(DupList2, file=paste(bedName, "silent.bed", sep=""), quote=F, row.names = F, col.names = F)
}

