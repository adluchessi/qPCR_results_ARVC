library(ggplot2)
getwd()
data<-read.csv("data_qPCR_ARVC_validation.txt", sep = "\t", header = T)

### Subset by miRNA for validation
data_mir=subset(data, Target %in% c("cel-miR-39-3p", "hsa-miR-92a-3p", "hsa-miR-16-5p", 
                                    "hsa-miR-15a-5p", "hsa-miR-145-5p", "hsa-miR-19a-3p", 
                                    "hsa-miR-29a-3p", "hsa-miR-505-3p"))
tail(data_mir)
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
write.csv(delta39,"delta39.csv", row.names = T)
dim(delta39)
susVA_5=read.delim("susVA_5.txt", header = T, sep = "\t")
delta39=cbind(delta39, susVA_5)

#CV IN %
sapply(delta39, function(x) sd(x, na.rm=T) / mean(x, na.rm=T)*100)

ggplot(delta39, aes(y=delta39[,2])) +
  geom_boxplot(varwidth = TRUE, alpha=0.2) +
  theme(legend.position="none")


plist=list()
for (i in c(2:4)) {
  plist[[i]]=boxplot(data=delta39, delta39[,i])
}


#####transfor in data matrix
library(corrplot)
delta39=delta39[complete.cases(delta39), ]
mm<- data.matrix(delta39[,-1])
dim(mm)
M<-cor(mm)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png(height=800, width=800, file="Correlation O.D. All Dataset.png", type = "cairo")
corrplot(M, method="circle", col=col(100),  
         diag=FALSE,
         type="lower", order="FPC", 
         addCoef.col = "black",
         mar=c(0,0,1,0),
         tl.srt = 45,
         number.cex= 2.5,
         tl.cex = 2.3,
         cl.cex = 2, number.digits = 2)
dev.off()



#BOXPLOT AND M.W. TEST
type<-c(rep('BrS',4),rep('Control',3))

delta39_all_data_fil_type<- cbind(delta39_all_data_fil, type)
delta39_all_data_fil_type$type <- factor(delta39_all_data_fil_type$type,
                                         levels=c("Control", "BrS")) #to put control first

for (i in c(2:18)) {
boxplot(data=delta39_all_data_fil_type, delta39_all_data_fil_type[,i]~type)
}
       
for (i in c(2:18)) {
 print(wilcox.test(delta39_all_data_fil_type[,i]~type, 
                   data=subset(delta39_all_data_fil_type, 
                               type %in% c("Control", "BrS"))))
}
















#HEATMAP

cat_df = data.frame("Classification"=delta39_all_data_fil_type$type)
head(cat_df, 12)
rownames(cat_df) = row.names(delta39_all_data_fil_type)

mm<- data.matrix(delta39_all_data_fil[,c(2:18)])
dim(mm)

pheatmap(t(mm), cutree_cols = 2, cluster_rows = F, cluster_cols = T,
         annotation_col = cat_df, show_colnames =F,
         cellwidth = 15, cellheight = 18, fontsize = 8, 
         main= "Heatmap 17 miRNAs expression in IPSC-CM between BrS and Control")

       






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
