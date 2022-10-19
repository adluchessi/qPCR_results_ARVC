library(ggplot2)
getwd()
data<-read.csv("data_qPCR_ARVC_validation_control.txt", sep = "\t", header = T)
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