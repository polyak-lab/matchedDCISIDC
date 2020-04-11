# script to re-run spat stat calculations, and redo the plots shown

library(spatstat)
library(divo)
library(RColorBrewer)
library(colorspace)
library(abind)
library(reshape2)
library(ggplot2)

a1=proc.time()

pal1=palette(brewer.pal(9, "Set3"))
cols_d4 <- colorspace::darken(pal1, 0.5)
pal2sub=rbind(pal1, cols_d4)
pal2sub=as.vector(pal2sub)

RdBu=brewer.pal(11,"RdBu")

load("data/RecurrenceCohort/cycIF/Final_cse8_data.RData")
sampid="8dcis"
palcols=pal1[c(4, 3:2, 5)]#cols_d4[c(4, 3:2, 5)]#pal1[c(4, 3:2, 5)]#cols_d4[c(4, 3:2, 5)]#pal1[c(4, 3:2, 5)]#cols_d4[c(4, 3:2, 5)] #edit this 

## define the groups
immlabels=c("bcell", "CD4", "cd8", "Tregs","macrophage", "pd1+")
stromaLabels=c("endothelial", "stroma", "vim+", "a-sma+")
lumLab2=c("CK+Ecad+" , "CK8 lo"  ,   "CKlow" ,  "col-hi" ,  "ER"   ,  "Ki67+"  ,  "myeop" ,"vim/sma", "basalLike")

memVals=case8_cellids$final2
refineArea=which(case8_cellids$type== "DCIS")# & case8_cellids$scene!="scene005")

## create spatstat object
msearch2=intersect(refineArea, which( memVals%in%lumLab2))
msearchB=intersect(refineArea, which( memVals=="CD4"))
msearchC=intersect(refineArea, which(  memVals=="cd8"))
msearchD=intersect(refineArea, which(  memVals=="Tregs"))
msearchE=intersect(refineArea, which(  memVals=="macrophage"))


windowOut=ripras(case8_cellids$XLocN[refineArea], case8_cellids$YLocN[refineArea])

newdfR=rbind(data.frame(case8_cellids[ msearch2,c("XLocN", "YLocN")], type="luminal"), 
            data.frame(case8_cellids[ msearchD,c("XLocN", "YLocN")], type="Treg"),
            data.frame(case8_cellids[ msearchC,c("XLocN", "YLocN")], type="CD8"),
            data.frame(case8_cellids[ msearchB,c("XLocN", "YLocN")], type="CD4"),
            data.frame(case8_cellids[ msearchE,c("XLocN", "YLocN")], type="Mac"))

xlimV=c(min(case8_cellids$XLocN[refineArea]), max(case8_cellids$XLocN[refineArea]))
ylimV=c(min(case8_cellids$YLocN[refineArea]), max(case8_cellids$YLocN[refineArea]))

ppOutR=ppp(newdfR$XLocN, newdfR$YLocN, marks=newdfR$type, poly=windowOut$bdry)

#######################
## interactinf distance analysis
######################
Ndist=30
d1=nndist(ppOutR, by=marks(ppOutR))
lx1a=which(ppOutR$marks=="Treg"& d1[ ,1]<Ndist)
lx1b=which(ppOutR$marks=="CD8"& d1[ ,1]<Ndist)
lx1c=which(ppOutR$marks=="CD4"& d1[ ,1]<Ndist)
lx1d=which(ppOutR$marks=="Mac"& d1[ ,1]<Ndist)

tcelldist=c(length(lx1a), length(lx1b),length(lx1c), length(lx1d))
amx1=table(newdfR$type)
Fracfound=tcelldist/amx1[-1]

## permute the labels
searchSet=intersect(refineArea, which( memVals%in%immlabels))
SumTab=table(memVals[refineArea])

DistMat=matrix(NA, nrow=1000, ncol=4)

a1=proc.time()

for (i in 1:1000){
  SampleSetcd4=sample(searchSet, SumTab[c("CD4")])
  SampleSetcd8=sample(setdiff(searchSet, SampleSetcd4), SumTab["cd8"])
  SampleSetFox=sample(setdiff(searchSet, c(SampleSetcd4, SampleSetcd8)), SumTab["Tregs"])
  SampleSetMac=sample(setdiff(searchSet, c(SampleSetcd4, SampleSetcd8, SampleSetFox)), SumTab["macrophage"])
  
  newdf=rbind(data.frame(case8_cellids[ msearch2, c("XLocN", "YLocN")], type="luminal"), 
              data.frame(case8_cellids[SampleSetFox,c("XLocN", "YLocN")], type="Treg"),
              data.frame(case8_cellids[ SampleSetcd8,c("XLocN", "YLocN")], type="CD8"),
              data.frame(case8_cellids[ SampleSetcd4,c("XLocN", "YLocN")], type="CD4"),
              data.frame(case8_cellids[ SampleSetMac,c("XLocN", "YLocN")], type="Mac"))
  
  ppOut=ppp(newdf$XLocN, newdf$YLocN, marks=newdf$type, poly=windowOut$bdry)
  
  d1=nndist(ppOut, by=marks(ppOut))
  lx1a=which(ppOut$marks=="Treg"& d1[ ,1]<Ndist)
  lx1b=which(ppOut$marks=="CD8"& d1[ ,1]<Ndist)
  lx1c=which(ppOut$marks=="CD4"& d1[ ,1]<Ndist)
  lx1d=which(ppOut$marks=="Mac"& d1[ ,1]<Ndist)
  lx2a=which(ppOut$marks=="luminal"& d1[ ,2]<Ndist)
  lx2b=which(ppOut$marks=="luminal"& d1[ ,3]<Ndist)
  lx2c=which(ppOut$marks=="luminal"& d1[ ,4]<Ndist)
  
  tdist=c(length(lx1a), length(lx1b),length(lx1c), length(lx1d))
  amx1=table(newdf$type)
  DistMat[i, ]=tdist/amx1[-1]
}

colnames(DistMat)=c("Treg", "CD8", "CD4", "Mac")
idc8_ci=sapply(1:4, function(x) quantile(DistMat[ ,x], c(0.025, 0.975)))
colnames(idc8_ci)=names(amx1)[-1]

Rx1=rbind(DistMat, Fracfound)
Rx2=scale(Rx1)
## calculate a p value:
pvals=2*pnorm(-abs(Rx2[1001, ]))

px1=melt(Rx2[-1001, ])
px2=melt(Rx2[1001,])
px2$Var2=rownames(px2)

pdf(sprintf("~/Desktop/%s_zscores.pdf", sampid), width=4, height=6)
ggplot(px1, aes(x=Var2, y=value, fill=Var2, col=Var2)) + geom_violin()+
  stat_summary(geom="point", shape=23, color="black", size=2)+
  annotate("text", label = round(pvals*100)/100, size = 3, x = px2$Var2, y = 4)+
  stat_summary(data=px2, geom="point", shape=5, col="red", size=2)+labs(y="z-score", x="cell type")+ 
  ggtitle(sampid)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values=palcols)+scale_color_manual(values=palcols)
dev.off()


## calculate the 1, 3, and 5 nearest neighbour distances and calculate the ks-statistic between pairs?
gridsize=250

ppInd=which(ppOutR$x>(xlimV[1]+gridsize) & ppOutR$x<(xlimV[2]-gridsize) & ppOutR$y>(ylimV[1]+gridsize) & 
              ppOutR$y<(ylimV[2]-gridsize))
#d1new=d1[ppInd, ]

d2=nndist(ppOutR, by=marks(ppOutR), k=c(1:5))

Dist1=d2[ ,grep("dist.1", colnames(d2))]
knn3=sapply(names(amx1), function(x) rowMeans(d2[ ,grep(x, colnames(d2))[1:3]]))
knn5=sapply(names(amx1), function(x) rowMeans(d2[ ,grep(x, colnames(d2))]))
knnDistMatrix=abind(Dist1, knn3, knn5, along=3)
Ctypes=levels(ppOutR$marks)[-1]

pdf(sprintf("~/Desktop/%s_knn_dist.pdf", sampid), width=9, height=5)
par(mfrow=c(1, 3))
#for (i in 1:3){
  i=1
  ptest=matrix(NA, nrow=4, ncol=4)
  namtest=matrix(NA, nrow=4, ncol=4)
  for (j in 1:(length(Ctypes)-1)){
    for (k in (j+1):length(Ctypes)){
    ptest[j,k]=ks.test(knnDistMatrix[ which(ppOutR$marks[ppInd]==Ctypes[j]),1, i], 
               knnDistMatrix[ which(ppOutR$marks[ppInd]==Ctypes[k]),1, i])$p.value
  #  ptest[j,k]=kSamples::ad.test(knnDistMatrix[ which(ppOutR$marks[ppInd]==Ctypes[j]),1, i], 
  #                     knnDistMatrix[ which(ppOutR$marks[ppInd]==Ctypes[k]),1, i])$p.value
    namtest[j,k]=paste(Ctypes[j], Ctypes[k])
    }}
  
  pvals=ptest[upper.tri(ptest)]
  names(pvals)=namtest[upper.tri(ptest)]
  
plot(ecdf(knnDistMatrix[ which(ppOutR$marks[ppInd]=="Treg"),1, i]), lwd=2, main=sprintf("k=%s", i*2-1 ), 
     xlim=c(0,300), col=palcols[1])
abline(v=30)
lines(ecdf(knnDistMatrix[ which(ppOutR$marks[ppInd]=="CD8"),1, i]), lwd=2, col=palcols[2])
lines(ecdf(knnDistMatrix[ which(ppOutR$marks[ppInd]=="CD4"),1, i]), lwd=2, col=palcols[3])
lines(ecdf(knnDistMatrix[ which(ppOutR$marks[ppInd]=="Mac"),1, i]), lwd=2, col=palcols[4])
legend("bottomright", legend=paste(names(pvals), round(pvals*1000)/1000), pch=1)
#}
dev.off()

###############
# Morisita Horn index testing
###############
#1. Test difference in grid size

gridsize=seq(100, 1000, by=100)

MHsumm=matrix(NA, nrow=length(gridsize), ncol=3*10)

for (i in 1:length(gridsize)){
  
  xcut=seq(xlimV[1], xlimV[2], length=round(diff(xlimV)/gridsize[i])+1)
  ycut=seq(ylimV[1], ylimV[2], length=round(diff(ylimV)/gridsize[i])+1)
  ppTess=tess(ppOutR, xgrid=xcut, ygrid=ycut)
  
  nt=quadratcount(ppOutR[ppOutR$marks=="luminal"], tess=ppTess)
  ni=quadratcount(ppOutR[ppOutR$marks=="Treg"], tess=ppTess)
  na=quadratcount(ppOutR[ppOutR$marks=="CD8"], tess=ppTess)
  ns=quadratcount(ppOutR[ppOutR$marks=="CD4"], tess=ppTess)
  nm=quadratcount(ppOutR[ppOutR$marks=="Mac"], tess=ppTess)
  
  allDat=data.frame(Tumor=as.vector(nt), Treg=as.vector(ni), CD8=as.vector(na), CD4=as.vector(ns), Mac=as.vector(nm))
  xres<-mh(allDat, resample=500)
  
  MHsumm[i, ]=c(xres$Mean[lower.tri(xres$Mean)], xres$Lower.Quantile[lower.tri(xres$Lower.Quantile)],xres$Upper.Quantile[lower.tri(xres$Upper.Quantile)])
}

colnames(MHsumm)=paste(rep(c("treg-tum", "cd8-tum", "cd4-tum" , "mac-tum","cd8-treg", "cd4-treg","mac-treg", "cd4-cd8", "mac-cd8", "mac-cd4"), 3), 
                       rep(c("mean", "low", "upper"), each=10 ))
rownames(MHsumm)=gridsize
lx=grep("mean", colnames(MHsumm))

MHmelt=cbind(melt(MHsumm[ , grep("mean", colnames(MHsumm))]), melt(MHsumm[ , grep("low", colnames(MHsumm))]), melt(MHsumm[ , grep("upper", colnames(MHsumm))]))

colnames(MHmelt)=c("xdist", "X1", "mean", "X2", "X3", "low", "X4", "X5", "upper")
MHmelt$cell=sapply(strsplit(as.character(MHmelt$X1), " "), function(x) x[1])


pdf(sprintf("~/Desktop/%s_MH_grid_variation.pdf", sampid), width=6, height=6)
ggplot(MHmelt, aes(x=factor(xdist), y=mean, col=cell))+geom_point()+geom_errorbar(aes(ymin=low, ymax=upper), width=.2)+scale_color_manual(values=pal2sub)
tempA=MHsumm[ 3,]
dfall=cbind(tempA[1:10], tempA[11:20], tempA[21:30 ])
dfall=data.frame(dfall)
dfall$comp=names(tempA)[1:10]
colnames(dfall)=c(  "mean", "lower", "upper", "comp")
dfall$type="gridfshift"
dev.off()

pdf(sprintf("~/Desktop/%s_MH_.pdf", sampid), width=6, height=6)
ggplot(dfall, aes(x=comp, y=mean))+geom_point()+geom_errorbar(aes(ymin=lower, ymax=upper), width=.2)+
  ggtitle("case8 idc Morisita-Horn variation of index")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# 2. shift the grid and test the distribution

gridsize=300
distx=seq(0, gridsize/2, by=5)
disty=seq(0,  gridsize/2, by=5)

distXb=seq( gridsize/2, 0, by=-5)
distYb=seq( gridsize/2, 0, by=-5)

MHsummGrid=matrix(NA, nrow=(length(distx)^2), ncol=30)
MHsummNoTum=matrix(NA, nrow=(length(distx)^2), ncol=30)
MHsummNoImm=matrix(NA, nrow=(length(distx)^2), ncol=30)
MHsummExcBoth=matrix(NA, nrow=(length(distx)^2), ncol=30)
count=1

for (i in 1:length(distx)){
  for (j in 1:length(distXb)){
    txlim=c(xlimV[1]+distx[i], xlimV[2]-distXb[i])
    tylim=c(ylimV[1]+disty[j], ylimV[2]-distYb[j]) 
    
    
    xcut=seq(txlim[1], txlim[2], length=round(diff(c(txlim))/gridsize)+1)
    ycut=seq(tylim[1], tylim[2], length=round(diff(c(tylim))/gridsize)+1)
    dfT=newdfR[which(newdfR$XLocN>txlim[1] & newdfR$XLocN <txlim[2] & newdfR$YLocN > tylim[1] & newdfR$YLocN < tylim[2] ), ]
    ppOut=ppp(dfT$XLocN, dfT$YLocN, marks=dfT$type, poly=windowOut$bdry)
    ppTess=tess(ppOut, xgrid=xcut, ygrid=ycut)
    
    nt=quadratcount(ppOut[ppOut$marks=="luminal"], tess=ppTess)
    ni=quadratcount(ppOut[ppOut$marks=="Treg"], tess=ppTess)
    na=quadratcount(ppOut[ppOut$marks=="CD8"], tess=ppTess)
    ns=quadratcount(ppOut[ppOut$marks=="CD4"], tess=ppTess)
    nm=quadratcount(ppOut[ppOut$marks=="Mac"], tess=ppTess)
    
    # use all points
    allDat=data.frame(Tumor=as.vector(nt), Treg=as.vector(ni), CD8=as.vector(na), CD4=as.vector(ns), Mac=as.vector(nm))
    xres<-mh(allDat, resample=10)
    MHsummGrid[count, ]=c(xres$Mean[lower.tri(xres$Mean)], xres$Lower.Quantile[lower.tri(xres$Lower.Quantile)],xres$Upper.Quantile[lower.tri(xres$Upper.Quantile)])
    # remove Imm-0
    xb=rowSums(allDat[ , -1])
    rmidx=which(allDat[ ,1]==0 & xb>0)
    allDat2=allDat[-c(rmidx), ]
    xres<-mh(allDat2, resample=10)
    MHsummNoImm[count, ]=c(xres$Mean[lower.tri(xres$Mean)], xres$Lower.Quantile[lower.tri(xres$Lower.Quantile)],xres$Upper.Quantile[lower.tri(xres$Upper.Quantile)])
   # rm Imm 0
    rmidx2=which(allDat[ ,1]>0 & xb==0)
    allDat2=allDat[-c(rmidx2), ]
    xres<-mh(allDat2, resample=10)
    MHsummNoTum[count, ]=c(xres$Mean[lower.tri(xres$Mean)], xres$Lower.Quantile[lower.tri(xres$Lower.Quantile)],xres$Upper.Quantile[lower.tri(xres$Upper.Quantile)])
    # rm both
    allDat2=allDat[-c(rmidx, rmidx2), ]
    xres<-mh(allDat2, resample=10)
    MHsummExcBoth[count, ]=c(xres$Mean[lower.tri(xres$Mean)], xres$Lower.Quantile[lower.tri(xres$Lower.Quantile)],xres$Upper.Quantile[lower.tri(xres$Upper.Quantile)])
     
   count=count+1
  }
}

colnames(MHsummGrid)=paste(rep(c("treg-tum", "cd8-tum", "cd4-tum" , "mac-tum","cd8-treg", "cd4-treg","mac-treg", "cd4-cd8", "mac-cd8", "mac-cd4"), 3), 
                       rep(c("mean", "low", "upper"), each=10 ))
colnames(MHsummNoTum)=colnames(MHsummGrid)
colnames(MHsummNoImm)=colnames(MHsummGrid)
colnames(MHsummExcBoth)=colnames(MHsummGrid)


MHexcBoth=sapply(1:10, function(x) quantile(MHsummExcBoth[ ,x], c(0.025, 0.975), na.rm = T))
dfallb=cbind(colMeans(MHsummExcBoth[ ,1:10], na.rm=T), MHexcBoth[1, ], MHexcBoth[2, ])
dfallb=data.frame(dfallb)
dfallb$comp=names(tempA)[1:10]
colnames(dfallb)=c(  "mean", "lower", "upper", "comp")
dfallb$type="shift_omit_both"

MHexcImm=sapply(1:10, function(x) quantile(MHsummNoImm[ ,x], c(0.025, 0.975), na.rm = T))
dfalli=cbind(colMeans(MHsummNoImm[ ,1:10], na.rm=T), MHexcImm[1, ], MHexcImm[2, ])
dfalli=data.frame(dfalli)
dfalli$comp=names(tempA)[1:10]
colnames(dfalli)=c(  "mean", "lower", "upper", "comp")
dfalli$type="shift_omit_imm"


MHexcTum=sapply(1:10, function(x) quantile(MHsummNoTum[ ,x], c(0.025, 0.975), na.rm = T))
dfallg=cbind(colMeans(MHsummNoTum[ ,1:10], na.rm=T), MHexcTum[1, ], MHexcTum[2, ])
dfallg=data.frame(dfallg)
dfallg$comp=names(tempA)[1:10]
colnames(dfallg)=c(  "mean", "lower", "upper", "comp")
dfallg$type="shift_omit_tum"

dfall3=rbind(dfall, dfallg, dfalli, dfallb)

# only use sample dfallg
dcisg=cbind(dfallg, samp="dcis")
idcg=cbind(dfallg, samp="idc")

df2=rbind(dcisg, idcg)



pdf(sprintf("~/Desktop/%s_MH_grid_exclude_samples.pdf", sampid), width=7, height=7)
ggplot(dfall3, aes(x=comp, y=mean, col=type))+geom_point()+geom_errorbar(aes(ymin=lower, ymax=upper), width=.2)+
  ggtitle(sampid)+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#idc8_null=DistMat
#idc8_knn=knnDistMatrix
save(DistMat,knnDistMatrix, idc8_ci,  MHsummGrid, MHsummNoTum, MHsummNoImm, MHsummExcBoth, dfall3, 
     file=sprintf("~/Desktop/%s_spatial.Rdata", sampid))

a2=proc.time()-a1


