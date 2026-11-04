library(DESeq2)
library(plot3D)
library(officer)
library(rvg)

TRAP_Sum <- read.delim("C:/Users/mw/Desktop/mw/Figure/Summarized figures/datasets/240903, TRAP_count.txt",
                       header = TRUE, skip = 1, check.names = FALSE)
Count_TRAP <- TRAP_Sum[, 7:ncol(TRAP_Sum)]
rownames(Count_TRAP) <- TRAP_Sum[, 1]
colnames(Count_TRAP) <- c("mRNA_GC", "mRNA_All", "Total_RNA", "S_GC_1", "S_GC_2", "S_GC_3",
                          "D_GC_1", "D_GC_2", "D_GC_3", "M_GC_1", "M_GC_2", "M_GC_3", "mRNA_PC")
Gene_ID <- rownames(Count_TRAP)

Count_all <- Count_TRAP[, c("Total_RNA", "mRNA_All", "mRNA_PC", "mRNA_GC",
                            "D_GC_1", "D_GC_2", "D_GC_3", "M_GC_1", "M_GC_2", "M_GC_3",
                            "S_GC_1", "S_GC_2", "S_GC_3"), drop = FALSE]

Count_GC <- Count_TRAP[, c("D_GC_1", "D_GC_2", "D_GC_3", "M_GC_1", "M_GC_2", "M_GC_3",
                           "S_GC_1", "S_GC_2", "S_GC_3"), drop = FALSE]

entire_sample <- data.frame(condition = factor(c("Total_RNA", "mRNA_All", "mRNA_PC", "mRNA_GC", "D_GC", "D_GC", "D_GC", "M_GC", "M_GC", "M_GC", "S_GC", "S_GC", "S_GC"), levels = c("Total_RNA", "mRNA_All", "mRNA_PC", "mRNA_GC", "D_GC", "M_GC", "S_GC")), row.names = colnames(Count_all))
Condition <- entire_sample$condition

# Figure 2b, 3b, MDS plot for all samples or GCs groups

dds_all <- DESeqDataSetFromMatrix(countData = Count_all, colData = entire_sample, design = ~ 1)
vsd_all <- vst(dds_all, blind = TRUE)
mds_all <- cmdscale(dist(t(assay(vsd_all))), k = 3, eig = TRUE)
dimension_percent_all <- round(100 * mds_all$eig / sum(mds_all$eig), 1)

GC_sample <- entire_sample[colnames(Count_GC), , drop = FALSE]
GC_sample$condition <- droplevels(GC_sample$condition)
Condition_GCs <- GC_sample$condition

dds_GCs <- DESeqDataSetFromMatrix(countData = Count_GC, colData = GC_sample, design = ~ 1)
vsd_GCs <- vst(dds_GCs, blind = TRUE)
mds_GCs <- cmdscale(dist(t(assay(vsd_GCs))), k = 3, eig = TRUE)
dimension_percent_GCs <- round(100 * mds_GCs$eig / sum(mds_GCs$eig), 1)

par(lwd = 2, mar = c(2, 2, 2, 4), xpd = NA)
scatter3D(x = mds_all$points[, 1], y = mds_all$points[, 2], z = mds_all$points[, 3],
          colvar = NULL, col = c(Total_RNA = "darkgray", mRNA_All = "black", mRNA_PC = "orange", mRNA_GC = "purple",
                                 D_GC = "red", M_GC = "blue", S_GC = "#008000")[as.character(entire_sample$condition)],
          pch = 16, cex = 0.8, cex.lab = 0.7, cex.axis = 0.7, theta = 320, phi = 25,
          bty = "u", col.grid = "lightgray", col.axis = "black", ticktype = "detailed",
          xlim = c(-20, 110), ylim = c(-80, 70), zlim = c(-30, 30),
          xlab = paste0("\nDimension 1 (", dimension_percent_all[1], "%)"),
          ylab = paste0("\nDimension 2 (", dimension_percent_all[2], "%)"),
          zlab = paste0("\nDimension 3 (", dimension_percent_all[3], "%)"))
legend("right", inset = c(-0.3, 0),
       legend = expression("Total RNA", "All", "PCs", "GCs", "D-GCs", "M-GCs", "S-GCs"),
       col = c(Total_RNA = "darkgray", mRNA_All = "black", mRNA_PC = "orange",
               mRNA_GC = "purple", D_GC = "red", M_GC = "blue", S_GC = "#008000"),
       pch = 16, pt.cex = 0.7, cex = 0.8, bty = "n")

par(lwd = 2, mar = c(2, 2, 2, 4), xpd = NA)
scatter3D(x = mds_GCs$points[, 1], y = mds_GCs$points[, 2], z = mds_GCs$points[, 3],
          colvar = NULL, col = c(D_GC = "red", M_GC = "blue", S_GC = "#008000")[as.character(Condition_GCs)],
          pch = 16, cex = 2, cex.lab = 0.7, cex.axis = 0.7, theta = 240, phi = 0,
          bty = "u", col.grid = "lightgray", col.axis = "black", ticktype = "detailed",
          xlim = c(-20, 30), ylim = c(-30, 30), zlim = c(-30, 30),
          xlab = paste0("\nDimension 1 (", dimension_percent_GCs[1], "%)"),
          ylab = paste0("\nDimension 2 (", dimension_percent_GCs[2], "%)"),
          zlab = paste0("\nDimension 3 (", dimension_percent_GCs[3], "%)"))
legend("top", inset = c(0, -0.1),
       legend = c("D-GCs", "M-GCs", "S-GCs"),
       col = c(D_GC = "red", M_GC = "blue", S_GC = "#008000"), pch = 16, pt.cex = 0.7, cex = 0.8, bty = "n")

ppt <- read_pptx()
ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")

ppt <- ph_with(ppt, dml(code = {
  par(lwd = 2, mar = c(2, 2, 2, 4), xpd = NA)
  scatter3D(x = mds_all$points[, 1], y = mds_all$points[, 2], z = mds_all$points[, 3],
            colvar = NULL, col = c(Total_RNA = "darkgray", mRNA_All = "black", mRNA_PC = "orange", mRNA_GC = "purple",
                                   D_GC = "red", M_GC = "blue", S_GC = "#008000")[as.character(entire_sample$condition)],
            pch = 16, cex = 0.8, cex.lab = 0.7, cex.axis = 0.7, theta = 320, phi = 25,
            bty = "u", col.grid = "lightgray", col.axis = "black", ticktype = "detailed",
            xlim = c(-20, 110), ylim = c(-80, 70), zlim = c(-30, 30),
            xlab = paste0("\nDimension 1 (", dimension_percent_all[1], "%)"),
            ylab = paste0("\nDimension 2 (", dimension_percent_all[2], "%)"),
            zlab = paste0("\nDimension 3 (", dimension_percent_all[3], "%)"))
  legend("right", inset = c(-0.3, 0),
         legend = expression("Total RNA", "All", "PCs", "GCs", "D-GCs", "M-GCs", "S-GCs"),
         col = c(Total_RNA = "darkgray", mRNA_All = "black", mRNA_PC = "orange",
                 mRNA_GC = "purple", D_GC = "red", M_GC = "blue", S_GC = "#008000"),
         pch = 16, pt.cex = 0.7, cex = 0.8, bty = "n")
}), location = ph_location(left = 0, top = 1.25, width = 4, height = 4))

ppt <- ph_with(ppt, dml(code = {
  par(lwd = 2, mar = c(2, 2, 2, 4), xpd = NA)
  scatter3D(x = mds_GCs$points[, 1], y = mds_GCs$points[, 2], z = mds_GCs$points[, 3],
            colvar = NULL, col = c(D_GC = "red", M_GC = "blue", S_GC = "#008000")[as.character(Condition_GCs)],
            pch = 16, cex = 2, cex.lab = 0.7, cex.axis = 0.7, theta = 240, phi = 0,
            bty = "u", col.grid = "lightgray", col.axis = "black", ticktype = "detailed",
            xlim = c(-20, 30), ylim = c(-30, 30), zlim = c(-30, 30),
            xlab = paste0("\nDimension 1 (", dimension_percent_GCs[1], "%)"),
            ylab = paste0("\nDimension 2 (", dimension_percent_GCs[2], "%)"),
            zlab = paste0("\nDimension 3 (", dimension_percent_GCs[3], "%)"))
  legend("top", inset = c(0, -0.1),
         legend = c("D-GCs", "M-GCs", "S-GCs"),
         col = c(D_GC = "red", M_GC = "blue", S_GC = "#008000"), pch = 16, pt.cex = 0.7, cex = 0.8, bty = "n")
}), location = ph_location(left = 5, top = 2.25, width = 3.5, height = 3.5))

print(ppt, target = "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Figure 2b, 3d 3D MDS plot_All and GCs.pptx")
