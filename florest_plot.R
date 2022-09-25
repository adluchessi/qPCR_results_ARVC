df=read.delim("or_qpcr.txt", header = T, sep = "\t")

#dt=df[,1:6]
dt=df[,1:6]
head(dt,25)
#head(dt1,7)
#dt=dt1
dt$MODELS <- ifelse(is.na(dt$High), 
                    dt$MODELS,
                    paste0("   ", dt$MODELS))


dt$Low <- ifelse(is.na(dt$Low), "", dt$Low)
dt$High <- ifelse(is.na(dt$High), "", dt$High)
dt$se <- (log(dt$hi) - log(dt$est))/1.96

dt$` ` <- paste(rep(" ", 30), collapse = " ")

# Create confidence interval column to display
dt$`HR (95% CI)` <- ifelse(is.na(dt$se), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$est, dt$low, dt$hi))

head(dt, 10)


tm <- forest_theme(base_size = 10,
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 16,
                   ci_col = "#762a83",
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2, # Set an T end at the end of CI 
                   # Reference line width/type/color
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   # Vertical line width/type/color
                   vertline_lwd = 1,
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   # Change summary color for filling and borders
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   # Footnote font size/face/color
                   footnote_cex = 0.6,
                   footnote_fontface = "italic",
                   footnote_col = "blue")

p <- forest(dt[,c(1:3, 8:9)],
            est = dt$est,
            lower = dt$low, 
            upper = dt$hi,
            sizes = dt$se,
            ci_column = 4,
            ref_line = 1,
            arrow_lab = c("Low Risk", "High Risk"),
            xlim = c(0, 10),
            ticks_at = c(1, 2, 5, 7, 10),
            #footnote = "This is the demo data. Please feel free to change\nanything you want.",
            #is_summary = c(rep(FALSE, nrow(dt)-1), TRUE),
            theme=tm)
# Print plot
#plot(p)

g=add_underline(p, part = "header")
#g=insert_text(g, text = "Risk Score", col = 2:3, part = "header", gp= gpar(fontface = "bold"))
g=edit_plot(g, row = c(1, 5,9,13,15, 16, 8, 20, 21, 22, 25, 26, 27,29,30,31,32), gp=gpar(fontface = "bold"))


plot(g)

