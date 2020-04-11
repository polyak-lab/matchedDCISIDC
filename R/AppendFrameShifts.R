AppendFrameShifts=function(getSeq, AnnotFS, fname){
  if (!is.na(AnnotFS)){
  FShifts=read.delim(AnnotFS, sep="\t")
  
  lx1=which(FShifts$Region=="INTRON")
  
  if (length(lx1)>0){
   FShifts=FShifts[-lx1, ] 
  }
  
  AAchange2=strsplit(as.character(FShifts$Amino.acid.change), "->")
  AAchange2=sapply(AAchange2, function(x) x[2])
  AAchange2New=strsplit(AAchange2, "\\*|\\.\\.\\.")
  AAchange2New=sapply(AAchange2New, function(x) x[1])
  PosEffect=strsplit(as.character(FShifts$Coordinates), ",")
  Chr1=sapply(PosEffect, function(x) x[1])
  Pos1=as.numeric(sapply(PosEffect, function(x) x[2]))+1
  #Pos1[grep("-/", FShifts$Coordinates)]=Pos1[grep("-/", FShifts$Coordinates)]+1
  
  matchIDX=match(Pos1, getSeq$Start_position)
  BB=getSeq$SequenceOfInterest2[matchIDX]
  BC=strsplit(BB, "_")
  BCD=sapply(BC, function(x) x[1])
  sidxrep=sapply(1:length(BB), function(x) regexpr(substr(AAchange2New[x],1, 4), BCD[x]))
  NSeq=toupper(paste(substr(BCD, 1, sidxrep-1), AAchange2New, sep=""))
  getSeq$SequenceOfInterest2[matchIDX]=NSeq
  }
  mx1=grep("NA[A-Z]NA", getSeq$SequenceOfInterest2)
  if (length(mx1)>0){
    seqinr::write.fasta(as.list(getSeq$SequenceOfInterest2[-mx1]), paste(getSeq$Hugo_Symbol[-mx1], getSeq$Start_position[-mx1], sep="_"),
                file.out = paste(fname, ".fasta", sep=""))
    
  }else{
    seqinr::write.fasta(as.list(getSeq$SequenceOfInterest2), paste(getSeq$Hugo_Symbol, getSeq$Start_position, sep="_"),
                file.out = paste(fname, ".fasta", sep=""))
  }
  write.table(getSeq, paste(fname, "mutation_list.txt", sep=""), sep="\t", row.names = F, quote=F)
  return(getSeq)
}
