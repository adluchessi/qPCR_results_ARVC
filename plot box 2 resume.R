library(ggstatsplot)
library(hrbrthemes)
library(ggplot2)
library(ggpubr)


data_arvc<-read.csv("data_qPCR_ARVC_validation_clinical_all.txt", sep = "\t", header = T)
head(data_arvc)


plt_1 <- ggbetweenstats(
  data = data_arvc,
  x = T_Nsus5,
  y =hsa.miR.145.5p, type = "np",
  pairwise.display="s",
  p.adjust.method="none",
  results.subtitle = FALSE
)

plt_1 <- plt_1 + 
  # Add labels and title
  labs(
    x = "5-year event-free survival group",
    y = " hsa-mir-145-5p expression",
    #title = "Distribution of hsa-mir-154-5p expression across tercil Nsus5 risk"
  ) + 
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )
plt_1


#################################################################

library(tidyverse)
library(readxl)
library(ggplot2)
library(glue)
library(ggtext)
library(ggpubr)

data_arvc<-read.csv("data_qPCR_ARVC_validation_clinical_all.txt", sep = "\t", header = T)
summary(data_arvc)
str(data_arvc)

T1_n = data_arvc %>% 
  filter(T_Nsus5 == "T1") %>% 
  nrow()
T2_n = data_arvc %>% 
  filter(T_Nsus5 == "T2") %>% 
  nrow()
T3_n = data_arvc %>% 
  filter(T_Nsus5 == "T3") %>% 
  nrow()
T4_n = data_arvc %>% 
  filter(T_Nsus5 == "Control") %>% 
  nrow()

T1_color="#BEBEBE"
T2_color="#0000FF"
T3_color="#FF0000"
T4_color="gray"


p1 <- data_arvc %>%
  ggplot(aes(x=T_Nsus5, hsa.miR.145.5p, fill= T_Nsus5))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.5, width=0.6, coef=0)+
  geom_jitter(show.legend = F, width=0.25, shape=21, color="black")+
  stat_summary(fun = median, show.legend=F, 
               geom='crossbar', width=0.6, size=0.5, color='black')+
  labs(x=NULL, 
       y="hsa-miR-145-5p")+
  scale_x_discrete(breaks = c("T1", "T2", "T3", "Control"), 
                   labels = c(glue("T1 (N={T1_n})"),
                              glue("T2 (N={T2_n})"), 
                              glue("T3 (N={T3_n})"),
                              glue("Control (N={T4_n})"))
  )+
  scale_fill_manual(name = NULL,
                    values =c(T1_color, T2_color, T3_color,T4_color)
  )+
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black")
  )

compare_means(hsa.miR.145.5p ~ T_Nsus5, data = data_arvc, method = "wilcox.test")
compare_means(hsa.miR.145.5p ~ T_Nsus5, data = data_arvc, method = "kruskal.test")

my_comparisons <- list( c("T1", "T2"), c("T1", "T3"), c("T2", "T3"), c("Control", "T3"))
p1 = p1 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.61)     # Add global p-value
p1

p2 <- data_arvc %>%
  ggplot(aes(x=T_Nsus5, hsa.miR.16.5p, fill= T_Nsus5))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.5, width=0.6, coef=0)+
  geom_jitter(show.legend = F, width=0.25, shape=21, color="black")+
  stat_summary(fun = median, show.legend=F, 
               geom='crossbar', width=0.6, size=0.5, color='black')+
  labs(x=NULL, 
       y="hsa-miR-16-5p")+
  scale_x_discrete(breaks = c("T1", "T2", "T3", "Control"), 
                   labels = c(glue("T1 (N={T1_n})"),
                              glue("T2 (N={T2_n})"), 
                              glue("T3 (N={T3_n})"),
                              glue("Control (N={T4_n})"))
  )+
  scale_fill_manual(name = NULL,
                    values =c(T1_color, T2_color, T3_color, T4_color)
  )+
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black")
  )

compare_means(hsa.miR.16.5p ~ T_Nsus5, data = data_arvc, method = "wilcox.test")
compare_means(hsa.miR.16.5p ~ T_Nsus5, data = data_arvc, method = "kruskal.test")

my_comparisons <- list( c("T1", "T2"), c("T1", "T3"), c("T2", "T3"), c("Control", "T3"))
p2 = p2 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 14)     # Add global p-value
p2


p3 <- data_arvc %>%
  ggplot(aes(x=T_Nsus5, hsa.miR.15a.5p, fill= T_Nsus5))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.5, width=0.6, coef=0)+
  geom_jitter(show.legend = F, width=0.25, shape=21, color="black")+
  stat_summary(fun = median, show.legend=F, 
               geom='crossbar', width=0.6, size=0.5, color='black')+
  labs(x=NULL, 
       y="hsa-miR-15a-5p")+
  scale_x_discrete(breaks = c("T1", "T2", "T3","Control"), 
                   labels = c(glue("T1 (N={T1_n})"),
                              glue("T2 (N={T2_n})"), 
                              glue("T3 (N={T3_n})"),
                              glue("Control (N={T4_n})"))
  )+
  scale_fill_manual(name = NULL,
                    values =c(T1_color, T2_color, T3_color, T4_color)
  )+
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black")
  )

compare_means(hsa.miR.16.5p ~ T_Nsus5, data = data_arvc, method = "wilcox.test")
compare_means(hsa.miR.16.5p ~ T_Nsus5, data = data_arvc, method = "kruskal.test")

my_comparisons <- list( c("T1", "T2"), c("T1", "T3"), c("T2", "T3"),c("Control", "T3"))
p3 = p3 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6)     # Add global p-value
p3


p4 <- data_arvc %>%
  ggplot(aes(x=T_Nsus5, hsa.miR.92a.3p, fill= T_Nsus5))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.5, width=0.6, coef=0)+
  geom_jitter(show.legend = F, width=0.25, shape=21, color="black")+
  stat_summary(fun = median, show.legend=F, 
               geom='crossbar', width=0.6, size=0.5, color='black')+
  labs(x=NULL, 
       y="hsa-miR-92a-3p")+
  scale_x_discrete(breaks = c("T1", "T2", "T3","Control"), 
                   labels = c(glue("T1 (N={T1_n})"),
                              glue("T2 (N={T2_n})"), 
                              glue("T3 (N={T3_n})"),
                              glue("Control (N={T4_n})"))
  )+
  scale_fill_manual(name = NULL,
                    values =c(T1_color, T2_color, T3_color, T4_color)
  )+
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black")
  )

compare_means(hsa.miR.16.5p ~ T_Nsus5, data = data_arvc, method = "wilcox.test")
compare_means(hsa.miR.16.5p ~ T_Nsus5, data = data_arvc, method = "kruskal.test")

my_comparisons <- list( c("T1", "T2"), c("T1", "T3"), c("T2", "T3"),c("Control", "T3") )
p4 = p4 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 28)     # Add global p-value
p4

p5 <- data_arvc %>%
  ggplot(aes(x=T_Nsus5, hsa.miR.19a.3p, fill= T_Nsus5))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.5, width=0.6, coef=0)+
  geom_jitter(show.legend = F, width=0.25, shape=21, color="black")+
  stat_summary(fun = median, show.legend=F, 
               geom='crossbar', width=0.6, size=0.5, color='black')+
  labs(x=NULL, 
       y="hsa-miR-19a-3p")+
  scale_x_discrete(breaks = c("T1", "T2", "T3", "Control"), 
                   labels = c(glue("T1 (N={T1_n})"),
                              glue("T2 (N={T2_n})"), 
                              glue("T3 (N={T3_n})"),
                              glue("Control (N={T4_n})"))
  )+
  scale_fill_manual(name = NULL,
                    values =c(T1_color, T2_color, T3_color, T4_color)
  )+
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black")
  )

compare_means(hsa.miR.16.5p ~ T_Nsus5, data = data_arvc, method = "wilcox.test")
compare_means(hsa.miR.16.5p ~ T_Nsus5, data = data_arvc, method = "kruskal.test")

my_comparisons <- list( c("T1", "T2"), c("T1", "T3"), c("T2", "T3"), c("Control", "T3"))
p5 = p5 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1.3)     # Add global p-value
p5

p6 <- data_arvc %>%
  ggplot(aes(x=T_Nsus5, hsa.miR.29a.3p, fill= T_Nsus5))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha=0.5, width=0.6, coef=0)+
  geom_jitter(show.legend = F, width=0.25, shape=21, color="black")+
  stat_summary(fun = median, show.legend=F, 
               geom='crossbar', width=0.6, size=0.5, color='black')+
  labs(x=NULL, 
       y="hsa-miR-29a-3p")+
  scale_x_discrete(breaks = c("T1", "T2", "T3", "Control"), 
                   labels = c(glue("T1 (N={T1_n})"),
                              glue("T2 (N={T2_n})"), 
                              glue("T3 (N={T3_n})"),
                              glue("Control (N={T4_n})"))
  )+
  scale_fill_manual(name = NULL,
                    values =c(T1_color, T2_color, T3_color, T4_color)
  )+
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black")
  )

compare_means(hsa.miR.16.5p ~ T_Nsus5, data = data_arvc, method = "wilcox.test")
compare_means(hsa.miR.16.5p ~ T_Nsus5, data = data_arvc, method = "kruskal.test")

my_comparisons <- list( c("T1", "T2"), c("T1", "T3"), c("T2", "T3"), c("Control", "T3"))
p6 = p6 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.32)     # Add global p-value
p6

ggarrange(p1, p2, p3, p4, p5, p6, 
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)
