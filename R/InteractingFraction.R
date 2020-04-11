
# refineArea=which(all9meta$type=="IDC" & all9meta$scene!="scene003")

InteractingFraction=function(memVals, refineArea, all9meta, labs=c("CD4", "CD8", "Foxp3", "Macrophage"), N=30){
  ## do a distance measurement between Tregs and luminal cells
  msearch2=intersect(refineArea, which( memVals%in%lumLab2))
  msearchB=intersect(refineArea, which( memVals==labs[1]))
  msearchC=intersect(refineArea, which(  memVals==labs[2]))
  msearchD=intersect(refineArea, which(  memVals==labs[3]))
  msearchE=intersect(refineArea, which(  memVals==labs[4]))
  
  windowOut=ripras(all9meta$XLocN[refineArea], all9meta$YLocN[refineArea])
  
  newdfR=rbind(data.frame(all9meta[ msearch2,c("XLocN", "YLocN")], type="luminal"), 
               data.frame(all9meta[ msearchD,c("XLocN", "YLocN")], type="Treg"),
               data.frame(all9meta[ msearchC,c("XLocN", "YLocN")], type="CD8"),
               data.frame(all9meta[ msearchB,c("XLocN", "YLocN")], type="CD4"),
               data.frame(all9meta[ msearchE,c("XLocN", "YLocN")], type="Mac"))
  
  ppOutR=ppp(newdfR$XLocN, newdfR$YLocN, marks=newdfR$type, poly=windowOut$bdry)
  
  ############
  ## method 1: calculate the nearest neighbour distance between immune and luminal, and compute the % which lie within 50pixels
  ############
  
  Ndist=30
  d1=nndist(ppOutR, by=marks(ppOutR))
  lx1a=which(ppOutR$marks=="Treg"& d1[ ,1]<Ndist)
  lx1b=which(ppOutR$marks=="CD8"& d1[ ,1]<Ndist)
  lx1c=which(ppOutR$marks=="CD4"& d1[ ,1]<Ndist)
  lx1d=which(ppOutR$marks=="Mac"& d1[ ,1]<Ndist)
  lx2a=which(ppOutR$marks=="luminal"& d1[ ,2]<Ndist)
  lx2b=which(ppOutR$marks=="luminal"& d1[ ,3]<Ndist)
  lx2c=which(ppOutR$marks=="luminal"& d1[ ,4]<Ndist)
  
  tcelldist=c(length(lx1a), length(lx1b),length(lx1c), length(lx1d))
  amx1=table(newdfR$type)
  Fracfound=tcelldist/amx1[-1]
  
  return(list(Fracfound=Fracfound, ppOut=ppOutR))
}