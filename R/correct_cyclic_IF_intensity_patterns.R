## adjust background intensities
## create a matrix, which is shifted by 10, 10 in the x and y directions.
## use spatstat package to correct for background shading
## show an example
## use the updated values for follow up - re-cluster with this new information

# Pre-load data of interest, and adjust data to "asinh" place
#

## pick which data to work with: DCIS or IDC
mx1=which(all8$type=="DCIS")

# by default, the color palette used:
colsb=colorRampPalette(brewer.pal(9, "RdYlGn"))

## set the starting values, and adjust this matrix to have a 10 pixel shift
## settings that have previously worked
## Sample| TileX | Tile Y| xhisft | yshift | xstart loc | y start loc
## 8 DCIS |1845 |1845 |13 | 0 |35900 |59400
## 8 IDC | 1844 | 1857+5 | 18 | 5 | 21350 | max(all8$YLocN[mx1])

TileX=1845 #+10 ## 
TileY=1845 # +5

xshift=13#17
yshift=0 # 12

diff1a=seq(35900, 0, by=-TileX) # or 1840?
diff1b=seq(59400, 0, by=-TileY)
## reference values: DCIS8, start y at 59400seq(59400, 0, by=-1840), start x at 36200 May need to adjust this value.

diffmatX=matrix(NA, ncol=length(diff1a), nrow=length(diff1b))
diffmatX[1, ]=diff1a
for (i in 2:nrow(diffmatX)){
  diffmatX[i ,]=diffmatX[i-1, ]-xshift
}

diffmatX=cbind(diffmatX,diffmatX[ ,ncol(diffmatX)]-TileX)
diffmatX=rbind(diffmatX,diffmatX[ nrow(diffmatX),]-TileY)

diffmatY=matrix(rep(diff1b, length(diff1a), byrow=F), ncol=length(diff1a))
for (i in 2:ncol(diffmatY)){
  diffmatY[,i]=diffmatY[,i-1 ]-yshift
}

diffmatY=rbind(diffmatY, diffmatY[nrow(diffmatY), ]-TileX)
diffmatY=cbind(diffmatY, diffmatY[, ncol(diffmatY)]-TileY)

# Plot an example and set the initial background values

stain="pAkt_Ring"
#x1=densCols(all8sinh[,stain], all8sinh[,stain],colramp = colsb)

#ii <- cut(all8sinh[,stain], breaks = c(min(all8sinh[,stain]), 
#                                            quantile(all8sinh[,stain], seq(0.2, 0.95, by=0.05)),
#                                            max(all8sinh[,stain])),include.lowest = TRUE)

ii <- cut(all8sinh[,stain], breaks = c( min(all8sinh[, stain])-0.1,seq( quantile(all8sinh[,stain], 0.05), 
                                          quantile(all8sinh[,stain], 0.95), len = 98),max(all8sinh[,stain])+0.1 ),include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(brewer.pal(9, "RdYlGn"))(99)[ii]


# x1=cut(all8sinh[,"PD1_Ring"], c(-1, quantile(all8sinh[,"PD1_Ring"], seq(0.1, 0.9, by=0.1)), max(all8sinh[,"PD1_Ring"])),
#        brewer.pal(10, "RdYlGn"))
plot(all8$XLocN[mx1],all8$YLocN[mx1], col=as.character(colors[mx1]), pch='.') 
# #,
# #     xlim=c(30000, 36000), ylim=c(50000, 60000))
 points(as.vector(diffmatX),as.vector(diffmatY))

#boxplot(all8sinh[,stain]~colors, las=2)
#, xlim=c(10000, 15000),
#     ylim=c(30000, 40000)) 


## barplot of the color scheme
# pickCols=names(table(x1))
# Vals=all8sinh[match(pickCols, x1) ,stain]



####
# If the above points roughly match the locations of differences, proceed with analysis

library(spatstat)

# give each point a new x and y co-ordinate
all8$XLocM=NA
all8$YLocM=NA

for (i in 1:(ncol(diffmatX)-1)){
  for (j in 1:(nrow(diffmatX)-1)){
    lx1=which(all8$XLocN<diffmatX[j, i] & all8$XLocN>=diffmatX[j,i+1])
    lx2=which(all8$YLocN<diffmatY[j,i] & all8$YLocN>=diffmatY[j+1,i])
    lx3=intersect(lx1, lx2)
    if (length(lx3)>0){
      all8$XLocM[lx3]=all8$XLocN[lx3]-diffmatX[j, i+1]
      all8$YLocM[lx3]=all8$YLocN[lx3]-diffmatY[j+1,i]
    }
  }
}

# example of what this looks like:
plot(all8$XLocM[mx1],all8$YLocM[mx1], col=as.character(colors[mx1]), main=stain,pch='.')

# perform smoothing operation
backgrd1=ppp(all8$XLocM[mx1], all8$YLocM[mx1],c(0, TileX+1), c(0,TileY+1), marks=all8sinh[mx1, stain])
# lx2=Smooth(backgrd1, sigma=10)
lx2=Smooth(backgrd1, sigma = 10, at="points")
PD1subtr=all8sinh[mx1[-which(is.na(all8$XLocM[mx1]))], stain]-lx2
# 
# par(mfrow=c(2,1))
# 
 
 # PlotColsMod(tx1, all8[ ,c("XLocN", "YLocN")],c(min(tx1[which(tx1>0)]), max(tx1)), which(all8$type=="DCIS"), 
 #             colnames(all8sinh)[i], pch='.')
 # 
 # x1=densCols(all8sinh[mx1[-which(is.na(all8$XLocM[mx1]))], stain], colramp = colsb) #[-which(is.na(all8$XLocM[mx1]))]
 # plot(all8$XLocN[mx1[-which(is.na(all8$XLocM[mx1]))]],all8$YLocN[mx1[-which(is.na(all8$XLocM[mx1]))]], col=x1, pch='.')
 # #xlim=c(10000, 15000), ylim=c(30000, 40000)) 
  #    xlim=c(17000, 32000), ylim=c(33000, 50000))
 Out2=PlotColsMod(all8sinh[ ,stain], all8[ ,c("XLocN", "YLocN")],quantile(all8sinh[ ,stain], c(0.025, 0.975)), mx1[-which(is.na(all8$XLocM[mx1]))], 
                  colnames(all8sinh)[i], pch='.')
 
# 
Out2=PlotColsMod(PD1subtr, all8[ mx1[-which(is.na(all8$XLocM[mx1]))],c("XLocN", "YLocN")],
                 quantile(PD1subtr, c(0.025, 0.975)), c(1:length(PD1subtr)), stain, pch='.')
 
x1=densCols(PD1subtr, colramp = colsb)
plot(all8$XLocN[mx1[-which(is.na(all8$XLocM[mx1]))]],all8$YLocN[mx1[-which(is.na(all8$XLocM[mx1]))]], col=x1, pch='.')#,
# xlim=c(10000, 15000), ylim=c(30000, 40000))
# 
# 
# ## If all above is good!
# ## Apply this process to all columns in the matrix!
# 
SampMod2=all8sinh[mx1[-which(is.na(all8$XLocM[mx1]))], ]
SampMod=all8sinh[mx1[-which(is.na(all8$XLocM[mx1]))], ] #[-which(is.na(all8$XLocM[mx1]))]
allfeat2=all8[mx1[-which(is.na(all8$XLocM[mx1]))], c(1, 43:50) ]
#
 for (i in 1:ncol(all8sinh)){
   backgrd1=ppp(allfeat2$XLocM, allfeat2$YLocM,c(0, TileX+1), c(0,TileY+1), marks=SampMod2[,i])
   lx2=Smooth(backgrd1, sigma = 10, at="points")
   SampMod[, i]=SampMod[, i]-lx2
 }

# i=23
# Out2=PlotColsMod(SampMod[ ,i], allfeat1[ ,c("XLocN", "YLocN")],
#                  quantile(SampMod[ ,i], c(0.025, 0.975)), c(1:length(SampMod[ ,i])), colnames(SampMod)[i], pch='.')
# 
# Out2=PlotColsMod(SampMod2[ ,i], allfeat1[ ,c("XLocN", "YLocN")],
#                  quantile(SampMod2[ ,i], c(0.025, 0.975)), c(1:length(SampMod2[ ,i])), colnames(SampMod2)[i], pch='.')
# 
case8dcisMod=SampMod
 
save(case8dcisMod, allfeat2, file="case8_dcis_manual_adjustment.RData")

## all data

case8_bck_mod=rbind(case8dcisMod, case8idcMod)
case8_cellids=rbind(allfeat2, allfeat1)

save(case8_bck_mod, case8_cellids, file="case8_0108_back_mod.RData")

# featidc8=all8[mx1[-which(is.na(all8$XLocM[mx1]))], ] #[-which(is.na(all8$XLocM[mx1]))]
# save(SampMod, featidc8, file="case8_idc_asinh_modified.RData")
# 
# 
# ## Load results and check the intensity distribution
# load("case8_dcis_asinh_modified.RData")
# dcisMod=SampMod
# load("case8_idc_asinh_modified.RData")
# idcMod=SampMod
# 
# par(mfrow=c(3,4))
# 
# for (i in 1:ncol(dcisMod)){
#   plot(density(dcisMod[ ,i]), main=colnames(dcisMod)[i])
#   lines(density(idcMod[ ,i]), col="red")
# }
# 
# ## Load to check the original intensity distribution
# par(mfrow=c(3,4))
# mx1=which(all8$type=="DCIS")
# for (i in 1:ncol(all8sinh)){
#   plot(density(all8sinh[mx1 ,i]), main=colnames(dcisMod)[i])
#   lines(density(all8sinh[-mx1 ,i]), col="red")
# }
# 
# 
# ## save the files to run phenograph!
# case8_bckadj=rbind(dcisMod, idcMod)
# # also save the cell type
# lx1=rbind(featdcis8[ ,c("memLum", "memStr", "memcd45")], featidc8[ ,c("memLum", "memStr", "memcd45")])
# 
# save(case8_bckadj,lx1, file="case8_cycIF_bckadj.RData")
