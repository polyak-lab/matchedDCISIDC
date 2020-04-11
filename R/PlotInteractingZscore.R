PlotInteractingZscore=function(Fracfound, DistMat, title="test", cols=palcols){
  Rx1=rbind(DistMat, Fracfound)
  Rx2=scale(Rx1)
  ## calculate a p value:
  pvals=2*pnorm(-abs(Rx2[1001, ]))
  
  px1=melt(Rx2[-1001, ])
  px2=melt(Rx2[1001,])
  px2$Var2=rownames(px2)
  
  g1=ggplot(px1, aes(x=Var2, y=value, fill=Var2, col=Var2)) + geom_violin()+
    stat_summary(geom="point", shape=23, color="black", size=2)+
    annotate("text", label = round(pvals*100)/100, size = 3, x = px2$Var2, y = 4)+
    stat_summary(data=px2, geom="point", shape=5, col="red", size=2)+labs(y="z-score", x="cell type")+ 
    ggtitle(title)+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_manual(values=cols)+scale_color_manual(values=cols)
  g1
}