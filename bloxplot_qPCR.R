library(ggstatsplot)
library(hrbrthemes)
library(ggplot2)


data_arvc<-read.csv("data_qPCR_ARVC_validation_clinical.txt", sep = "\t", header = T)
head(data_arvc)

plt <- ggbetweenstats(
  data = data_arvc,
  x = T_Nsus5,
  y =hsa.miR.505.3p, type = "np",
  pairwise.display="s",
  p.adjust.method="none"
)

plt <- plt + 
  # Add labels and title
  labs(
    x = "Nsus5 tercil",
    y = "hsa-mir-505-3p relative expression",
    title = "Distribution of hsa-mir-505-3p expression across tercil Nsus5 risk"
  ) + 
  # Customizations
  theme(
    # This is the new default font in the plot
    text = element_text(family = "TT Courier New", size = 14, color = "black"),
    plot.title = element_text(
      family = "Lobster Two", 
      size = 18,
      face = "bold",
      color = "black",
      hjust = 0.5
    ),
    # Statistical annotations below the main title
    plot.subtitle = element_text(
      family = "Roboto", 
      size = 14, 
      face = "bold",
      color="black"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14)
  )

plt <- plt  +
  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )

plt

names(data_arvc)
