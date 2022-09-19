library(ggstatsplot)
library(hrbrthemes)
library(ggplot2)
library(ggpubr)


data_arvc<-read.csv("data_qPCR_ARVC_validation_clinical.txt", sep = "\t", header = T)
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
    x = "5-year event-free survival",
    y = "hsa-mir-145-5p",
    #title = "Distribution of hsa-mir-154-5p expression across tercil Nsus5 risk"
) + 
  theme(
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black", face = "bold"),
    
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
)





plt_2 <- ggbetweenstats(
  data = data_arvc,
  x = T_Nsus5,
  y =hsa.miR.15a.5p, type = "np",
  pairwise.display="s",
  p.adjust.method="none",
  results.subtitle = FALSE
)

plt_2 <- plt_2 + 
  # Add labels and title
  labs(
    x = "5-year event-free survival",
    y = "hsa-mir-15a-5p",
    #title = "Distribution of hsa-mir-15a-5p expression across tercil Nsus5 risk"
) + 
  theme(
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black", face = "bold"),
    
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
)



plt_3 <- ggbetweenstats(
  data = data_arvc,
  x = T_Nsus5,
  y =hsa.miR.19a.3p, type = "np",
  pairwise.display="s",
  p.adjust.method="none",
  results.subtitle = FALSE
)

plt_3 <- plt_3 + 
  # Add labels and title
  labs(
    x = "5-year event-free survival",
    y = "hsa-mir-19a-3p",
    #title = "Distribution of hsa-mir-19a-3p expression across tercil Nsus5 risk"
) + 
  theme(
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black", face = "bold"),
    
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
)






plt_4 <- ggbetweenstats(
  data = data_arvc,
  x = T_Nsus5,
  y =hsa.miR.92a.3p, type = "np",
  pairwise.display="s",
  p.adjust.method="none",
  results.subtitle = FALSE
)

plt_4 <- plt_4 + 
  # Add labels and title
  labs(
    x = "5-year event-free survival",
    y = "hsa-mir-92a-3p relative expression",
    #title = "Distribution of hsa-mir-92a-3p expression across tercil Nsus5 risk"
) + 
  theme(
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black", face = "bold"),
    
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
)





plt_5 <- ggbetweenstats(
  data = data_arvc,
  x = T_Nsus5,
  y =hsa.miR.16.5p, type = "np",
  pairwise.display="s",
  p.adjust.method="none",
  results.subtitle = FALSE
)

plt_5 <- plt_5 + 
  # Add labels and title
  labs(
    x = "5-year event-free survival",
    y = "hsa-mir-16-5p",
    #title = "Distribution of hsa-mir-16-5p expression across tercil Nsus5 risk"
) + 
  theme(
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black", face = "bold"),
    
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
)





plt_6 <- ggbetweenstats(
  data = data_arvc,
  x = T_Nsus5,
  y =hsa.miR.29a.3p, type = "np",
  pairwise.display="s",
  p.adjust.method="none",
  results.subtitle = FALSE
)

plt_6 <- plt_6 + 
  # Add labels and title
  labs(
    x = "5-year event-free survival",
    y = "hsa-mir-29a-3p",
    #title = "Distribution of hsa-mir-29a-3p expression across tercil Nsus5 risk"
) + 
  theme(
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black", face = "bold"),
    
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
)

ggarrange(plt_1, plt_2, plt_3, plt_4, plt_5, plt_6 + rremove("x.text"), 
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)
