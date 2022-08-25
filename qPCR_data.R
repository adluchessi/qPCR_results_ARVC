getwd()
library(pheatmap)
library(ggplot2)
getwd()
data<-read.csv("miRNA_iPSC-CM_BrS_qPCRresults.txt", sep = "\t", header = T)
head(data)

data$Sample <- factor(data$Sample, levels=c("10_2991_CT","11_3002_CT","12_3014_CT",
"14_3036_CT","7_2917_CT","9_2985_BrS","1_2747_BrS","13_2995_BrS","2_2748_BrS",
"4_2775_BrS","5_2790_BrS","6_2914_BrS"))

ggplot(data, aes(x=Sample, y=Cq.Mean, fill=Type)) +
  geom_boxplot(varwidth = TRUE, alpha=0.2) +
  theme(legend.position="none")

samples1<- as.data.frame(t(data[c(1:24),3]))
colnames(samples1)<- data[c(1:24),2]
rownames(samples1)<- data[1,1]

samples2<- as.data.frame(t(data[c(25:48),3]))
colnames(samples2)<- data[c(1:24),2]
rownames(samples2)<- data[25,1]

samples3<- as.data.frame(t(data[c(49:72),3]))
colnames(samples3)<- data[c(1:24),2]
rownames(samples3)<- data[49,1]

samples4<- as.data.frame(t(data[c(73:96),3]))
colnames(samples4)<- data[c(1:24),2]
rownames(samples4)<- data[73,1]

samples5<- as.data.frame(t(data[c(97:120),3]))
colnames(samples5)<- data[c(1:24),2]
rownames(samples5)<- data[97,1]

samples6<- as.data.frame(t(data[c(121:144),3]))
colnames(samples6)<- data[c(1:24),2]
rownames(samples6)<- data[121,1]

samples7<- as.data.frame(t(data[c(145:168),3]))
colnames(samples7)<- data[c(1:24),2]
rownames(samples7)<- data[145,1]

samples8<- as.data.frame(t(data[c(169:192),3]))
colnames(samples8)<- data[c(1:24),2]
rownames(samples8)<- data[169,1]

samples9<- as.data.frame(t(data[c(193:216),3]))
colnames(samples9)<- data[c(1:24),2]
rownames(samples9)<- data[193,1]

samples10<- as.data.frame(t(data[c(217:240),3]))
colnames(samples10)<- data[c(1:24),2]
rownames(samples10)<- data[217,1]

samples11<- as.data.frame(t(data[c(241:264),3]))
colnames(samples11)<- data[c(1:24),2]
rownames(samples11)<- data[241,1]

samples12<- as.data.frame(t(data[c(265:288),3]))
colnames(samples12)<- data[c(1:24),2]
rownames(samples12)<- data[265,1]

all_data<- rbind(samples1, samples2, samples3, samples4, samples5, samples6,
                 samples7, samples8, samples9, samples10, samples11, samples12)
dim(all_data)

#FILT BY MIRNAs

all_data_fil<- all_data[,-c(10,12,16,17,19,21)]
dim(all_data_fil)
head(all_data_fil)

#NORM BY miRNA-39

delta39_all_data_fil<- 2^-(all_data_fil-all_data_fil[,1])
delta39_all_data_fil<- delta39_all_data_fil[-12,] #Filt by samples
dim(delta39_all_data_fil)
#CV IN %
sapply(all_data_fil, function(x) sd(x, na.rm=T) / mean(x, na.rm=T) * 100)

#BOXPLOT AND M.W. TEST
type<-c('BrS', 'Control', rep('BrS',4), 'Control','BrS',
        rep('Control',3))

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
