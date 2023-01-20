library(ggplot2)
library(plotROC)

dnam <- read.csv("data_qPCR_ARVC_validation_clinical_roc.txt", sep = "\t", header = T)
row.names(dnam)<- make.names(dnam[,1], unique = T)
dnam<- dnam[,-1]
dnam=na.omit(dnam)
head(dnam)

test <- data.frame(D = dnam$D, M1=dnam$M1, M2=dnam$M2, M3=dnam$M3, M4=dnam$M4, M5=dnam$M5, M6=dnam$M6,
                   stringsAsFactors = FALSE)

longtest <- melt_roc(test, "D", c("M1", "M2","M3","M4","M5","M6"))
head(longtest)


p1 = ggplot(longtest, aes(d = D, m=M, color = name)) + geom_roc(n.cuts = 0, show.legend = F)+
  scale_color_manual(values = c("coral4", "green4", "blue4", "orange4","turquoise4","gray20"))

p1 + style_roc(major.breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
              minor.breaks = c(seq(0, 0.1, by = 0.01), seq(0.9, 1, by = 0.01)),
              guide = TRUE, xlab = "1 - Specificity (FPF)", ylab = "Sensitivity (TPF)") +
  
  annotate("text",x = .63, y = .55, label = paste("miRNA-16-5p = 0.970"),
           size =7,color="coral4")+
  annotate("text",x = .64, y = .47, label = paste("miRNA-92a-3p =", round(calc_auc(p1)$AUC, 3)[2]),
           size =7, color="green4")+
  annotate("text",x = .64, y = .39, label = paste("miRNA-19a-3p =", round(calc_auc(p1)$AUC, 3)[3]),
           size =7, color="blue4")+
  annotate("text",x = .64, y = .31, label = paste("miRNA-15a-5p =", round(calc_auc(p1)$AUC, 3)[4]),
           size =7, color="orange4")+
  annotate("text",x = .64, y = .23, label = paste("miRNA-145-5p = 0.700"),
           size =7, color="turquoise4")+
  annotate("text",x = .638, y = .15, label = paste("miRNA-29a-3p =", round(calc_auc(p1)$AUC, 3)[5]),
           size =7, color="gray20")+
  annotate("text",x = .79, y = .62, label = paste("AUC"),
           size =7, color="black")+
  theme(axis.title=element_text(size = 20),
        axis.text=element_text(size = 18))

