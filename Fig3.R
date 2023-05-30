library(forestplot)

#figure3
types <- read.csv("plot3.csv", header = F)
instant<-types[1:3,]
de<-types[c(1,5,6),]
filter<-types[c(1,8,9),]
ground<-types[c(1,11:14),]
for (i in 1:4){
  names<-c("3a","3b","3c","3d")
  lists<-list(instant,de,filter,ground)
  png(paste0("forest_plot/Fig",names[i],".png"),res=600,width = 36,height = 8,units = "cm")
  data<-lists[[i]]
  forestplot(as.matrix(data[,c(1:4)]),
           mean = data$V5, 
           lower = data$V6,
           upper = data$V7,
           zero = 0,
           xlog =F,
           xlab="beta",
           fn.ci_norm = fpDrawCircleCI,
           boxsize = 0.15, is.summary=c(T,rep(F,5)),
           col=fpColors(line = "blue",box="#D55E00"),
           lty.ci = 7,lwd.ci = 3, 
           ci.vertices.height = 0.05,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.7),
                            xlab =gpar(cex = 0.8),cex = 1.2),
           lineheight = "auto",
           align="l",
           graph.pos = 4,
           hrzl_lines=T)
  dev.off()
}
