library(forestplot)

#figure1
data <- read.csv("forest_plot/plot1.csv", header = F)
png("forest_plot/Fig2.png",res=600,width = 38,height = 12,units = "cm")
forestplot(as.matrix(data[,c(1:5)]),
           mean = data$V6, 
           lower = data$V7,
           upper = data$V8,
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
        
           