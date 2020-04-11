
# Annotate Mutations: determine the corresponding protein sequence

AnnotateMutations=function(TumSpec2,  fname="Sample1"){
## Install this if needed to annotate mutations   
 # library(EnsDb.Hsapiens.v86)
#  edb <- EnsDb.Hsapiens.v86
  testS=c( "Nonstop_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Nonsense_Mutation")
  ax=which(TumSpec2$Variant_Classification %in% testS) #Coding=="Yes" & TumSpec2$Sequence.Ontology!="synonymous")
  getSeq=TumSpec2[ax, ]
  getSeq$protID=sapply(strsplit(as.character(getSeq$HGVS_protein_change), ":"), function(x) x[1])
  fside=strsplit(as.character(getSeq$Protein_Change), "[0-9]+")
  Nocase=regmatches(getSeq$Protein_Change, regexpr("[0-9]+", getSeq$Protein_Change))
  fsideA=sapply(fside, function(x) gsub("p\\.", "", x[1]))
  fsideB=sapply(fside, function(x)  x[2])
  
  getSeq$AAno=as.numeric(Nocase)
  getSeq$OrigAA=fsideA
  getSeq$MutAA=fsideB
  n=10
  a1=sapply(strsplit(as.character(getSeq$Annotation_Transcript), "\\."), function(x) x[1])
  
  if (length(a1)==nrow(getSeq)){
    getSeq$T2=a1
  } 
  
  #txs<-transcripts(edb, filter=GenenameFilter(getSeq$Hugo), columns=c("uniprot_id", "tx_id"))
  prts<-proteins(edb, filter=ProteinIdFilter(na.omit(getSeq$protID)), return.type="AAStringSet")
  #proteins(edb, filter=UniprotDbFilter(getSeq$Hugo_Symbol), return.type="AAStringSet")
  #prts<-proteins(edb, filter=GenenameFilter(getSeq$Hugo_Symbol), return.type="AAStringSet")
  
  #getSeq$peptide=NA
  #getSeq$ensp=NA

  for (i in 1:nrow(getSeq)){
  test1=which(mcols(prts)$protein_id==getSeq$protID[i])
  test2=test1[which(width(prts[test1]) > getSeq$AAno[i])]
  tempP=prts[test2]
   
    if (length(test2)==0){
      uP=unlist(strsplit(as.character(getSeq$dbNSFP_Uniprot_acc[i]), ";"))
      if (length(uP)>0){
        tempP<-proteins(edb, filter=UniprotFilter(uP), return.type="AAStringSet")
      }
      test2=which(width(tempP) > getSeq$AAno[i] )
      tempP=tempP[test2]
    }

  mtest=subseq(tempP, start=getSeq$AAno[i], end=getSeq$AAno[i])
  
     if (getSeq$OrigAA[i]=="_"){
    getSeq$peptide[i]=as.character(tempP[1])
    getSeq$ensp[i]=names(tempP[1])
  } else {
    getV=which(as.character(mtest)==substr(getSeq$OrigAA[i],1,1))
      if (length(getV)>0){
       getSeq$peptide[i]=as.character(tempP[getV[1]])
       getSeq$ensp[i]=names(getV)[1]
      } 
  
  }
  }
  
 
   getSeq <- getSeq %>% 
    mutate(SequenceOfInterest = peptide %>% substr(AAno-n, AAno+n)) 
  
  getSeq <- getSeq %>% 
    mutate(SequenceOfInterest2 = paste(peptide %>% substr(AAno-n,AAno-1),MutAA, peptide %>% substr(AAno+1,AAno+n), sep=""))
  
   # cases with frameshifts? and NS
  getSeq$SequenceOfInterest2[getSeq$MutAA=="*"]=substr(getSeq$SequenceOfInterest2[which(getSeq$MutAA=="*")], 1, 
                                                       unlist(gregexpr(pattern ='\\*',getSeq$SequenceOfInterest2[which(getSeq$MutAA=="*")]))-1)
  
  gS2fs=getSeq[getSeq$MutAA=="fs", ]
  
  # strand inforation
  if (nrow(gS2fs)>0){
    gS2fs$ins=ifelse(gS2fs$Tumor_Seq_Allele2=="-", "/",as.character(gS2fs$Tumor_Seq_Allele2))
    #Pos2=ifelse(gS2fs$ins=="/", as.numeric(as.character(gS2fs$Position)),as.numeric(as.character(gS2fs$Position))-1)
    
    # need to change the site names for insertions and or deletions
    ## prepare information for website: http://sift.bii.a-star.edu.sg
    coord=paste(gS2fs$Chromosome,as.numeric(as.character(gS2fs$Start_position))-1, as.numeric(as.character(gS2fs$End_position)), 1, gS2fs$ins, sep="," )
    write.table(coord, file = paste(fname, "frameshift.txt", sep=""), row.names = F, col.names = F, quote=F)
  }
  return(getSeq)
}
