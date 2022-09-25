library(ggplot2)
getwd()
data<-read.csv("data_qPCR_ARVC_validation.txt", sep = "\t", header = T)
head(data)

### Subset by miRNA for validation
data_mir=subset(data, Target %in% c("cel-miR-39-3p", "hsa-miR-92a-3p", "hsa-miR-16-5p", 
                                    "hsa-miR-15a-5p", "hsa-miR-145-5p", "hsa-miR-19a-3p", 
                                    "hsa-miR-29a-3p", "hsa-miR-505-3p"))

### Boxplot CT for all miRNA expression distribution by samples
ggplot(data, aes(x=as.character(Sample), y=Cq.Mean, fill=Type)) +
  geom_boxplot(varwidth = TRUE, alpha=0.2) +
  theme(legend.position="none")

### Data transposition
my_splits_samples <- split(data_mir, data_mir$Sample)

l <- list()
for(i in 1:length(my_splits_samples)){
  df.now <- as.data.frame(t(my_splits_samples[[i]]))
  l[[i]] <- df.now
}

l_2=list()
l_3=list()
for(i in 1:length(l)){
  colnames(l[[i]])=l[[i]][2,]
  df.now <- l[[i]][-c(1,2),]
  l_2[[i]] <- df.now
  df.now <- l_2[[i]][-c(2),]
  l_3[[i]] <- df.now
  rownames(l_3[[i]])<- paste('Sample',l[[i]][1,1],sep = "_")
}

data_qPCR_ARVC=do.call(rbind,l_3)
head(data_qPCR_ARVC,15)
################################################################################

#NORM BY miRNA-39
data_qPCR_ARVC_num=as.data.frame(sapply(data_qPCR_ARVC, as.numeric))
delta39<- 2^-(data_qPCR_ARVC_num-data_qPCR_ARVC_num[,1])
rownames(delta39)=rownames(data_qPCR_ARVC)
head(delta39,15)
dim(delta39)

write.csv(delta39,"delta39.csv", row.names = T)

################################################################################

#Select sample only ARVC (TFC>3.5)

data_arvc<-read.csv("data_qPCR_ARVC_validation_clinical.txt", sep = "\t", header = T)
head(data_arvc)


#CV IN %
sapply(delta39, function(x) sd(x, na.rm=T) / mean(x, na.rm=T)*100)

ggplot(delta39, aes(y=delta39[,2])) +
  geom_boxplot(varwidth = TRUE, alpha=0.2) +
  theme(legend.position="none")


plist=list()
for (i in c(2:4)) {
  plist[[i]]=boxplot(data=delta39, delta39[,i])
}


##### Corr plot transfor in data matrix
library(corrplot)
data_arvc=data_arvc[complete.cases(data_arvc), ]
head(data_arvc)
mm<- data.matrix(data_arvc[,c(2:7)])
mm
dim(mm)
M<-cor(mm)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png(height=800, width=800, file="Correlation miRNAs ARVC samples.png", type = "cairo")
corrplot(M, method="circle", col=col(100),  
         diag=FALSE,
         type="lower", order="FPC", 
         addCoef.col = "black",
         mar=c(0,0,1,0),
         tl.srt = 45,
         number.cex= 2.1,
         tl.cex = 1.8,#tl.cex = legend texy size
         cl.cex = 2, number.digits = 2) #cl.cex = eixo x size text
dev.off()



#BOXPLOT AND M.W. TEST
data_arvc<-read.csv("data_qPCR_ARVC_validation_clinical.txt", sep = "\t", header = T)
dim(data_arvc)


#delta39_all_data_fil_type$type <- factor(delta39_all_data_fil_type$type,
                                         #levels=c("Control", "BrS")) #to put control first
names(data_arvc)
for (i in c(2:8)) {
boxplot(data=data_arvc, data_arvc[,i]~T_Nsus5)
}
       
for (i in c(2:18)) {
 print(wilcox.test(data_arvc[,i]~T_Nsus5, 
                   data=subset(data_arvc, 
                               T_Nsus5 %in% c("T1", "T3"))))
}





#HEATMAP
library(pheatmap)
data_arvc=data_arvc[complete.cases(data_arvc), ]
dim(data_arvc)
head(data_arvc)
mm<- data.matrix(data_arvc[,c(2:7)])
dim(mm)
mm
myannotation=as.data.frame(data_arvc$T_Nsus5)
myannotation

pheatmap(t(log(mm)), cutree_cols = 1, cluster_rows = T, cluster_cols = T,
        show_colnames =F,
        cellwidth = 5, cellheight = 18, fontsize = 8,
       # annotation_col=myannotation,
        main= "Heatmap 6 miRNAs diferently expression by ARVC Risk Score")



install.packages("circlize")
library(circlize)

circos.heatmap(log(mm), split = split, col = col_fun1, dend.side = "outside")
circos.clear()


#NORM BY 16

delta16_all_data_fil<- 2^-(all_data_fil-all_data_fil[,7])
delta16_all_data_fil<-delta16_all_data_fil[-2,] #samples off mir16

#BOxPLoT AND MW TEST
type<-c('Control', 'Control', rep('BrS',4), rep('Control',2), rep('BrS',2), 'Control')

delta16_all_data_fil_type<- cbind(delta16_all_data_fil, type)
delta16_all_data_fil_type$type <- factor(delta16_all_data_fil_type$type, levels=c("Control", "BrS"))
delta16_all_data_fil_type

for (i in c(2:18)){
  boxplot(data=delta16_all_data_fil_type, delta16_all_data_fil_type[,i]~type)
}

for (i in c(2:18)) {
  print(wilcox.test(delta16_all_data_fil_type[,i]~type, data=subset(delta16_all_data_fil_type, type %in% c("Control", "BrS"))))
}

#HEATMAPS

cat_df = data.frame("Classification"=delta16_all_data_fil_type$type)
head(cat_df)
rownames(cat_df) = row.names(delta16_all_data_fil)
cat_df
mm<- data.matrix(delta16_all_data_fil[,c(2:18)])
dim(mm)
mm
pheatmap(t(log(mm)), cutree_cols = 2, cluster_rows = T, cluster_cols = T,annotation_col = cat_df, show_colnames =F,
         cellwidth = 15, cellheight = 18, fontsize = 8, 
         main= "x")
