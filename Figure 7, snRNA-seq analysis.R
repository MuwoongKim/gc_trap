library(Seurat)
library(dplyr)
library(tibble)
library(readxl)
library(writexl)
library(ggplot2)
library(rvg)
library(officer)

granule <- readRDS("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Datasets/granule.RDS")
granule <- UpdateSeuratObject(granule)
DefaultAssay(granule) <- "RNA"
granule$subcluster <- factor(granule$subcluster, levels = c("Granule_1", "Granule_2", "Granule_3", "Granule_4", "Granule_5"))
color_code <- c("Granule_1" = "#00B0F6" , "Granule_2" = "#E76BF3", "Granule_3" = "#A3A500" , "Granule_4" = "#00BF7D" , "Granule_5" = "#F8766D")
UMAP_cerebellum <- DimPlot(granule, reduction = "umap", group.by = "subcluster", cols = color_code,
                           label = FALSE, label.size = 3, label.box = FALSE, raster = TRUE, raster.dpi = c(500, 500), pt.size = 1) +
  theme(plot.title = element_blank(), text = element_text(size = 5), legend.position = "none") + NoAxes()
UMAP_cerebellum <- LabelClusters(UMAP_cerebellum, id = "subcluster", size = 20 / .pt)
UMAP_cerebellum
ggsave("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/261103, GC clusters (Kozareva).tiff", UMAP_cerebellum, width = 5, height = 5, dpi = 600, compression = "lzw")
ggsave("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/261103, GC clusters (Kozareva).jpeg", UMAP_cerebellum, width = 5, height = 5, dpi = 600, quality = 100)

# Dotplot of marker genes (Confirmation)
genes_gran = c("Rasgrf1", "Gprin3", "Ebf1", "Galntl6", "Nr4a2", "Chrm3")
Idents(granule) <- "subcluster"
Marker_genes <- DotPlot(granule, features = genes_gran, group.by = "subcluster", dot.scale = 5, scale.min = 0.1) +
  guides(size = guide_legend(title = "Exp.(%)"), color = guide_colorbar(title = "Scaled exp")) +
  theme(plot.title = element_blank(), text = element_text(size = 7.5), axis.text = element_text(size = 7.5),
        axis.title = element_blank(), legend.position = "right", legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        axis.line = element_line(linewidth = 0.1, colour = "black"),
        axis.ticks = element_line(linewidth = 0.1)) + RotatedAxis()
Marker_genes

# Selection process
Candidate_genes <- read_excel("C:/Users/mw/Desktop/mw/1.TRAP/260116, GC_DEGs_Total_Comparison.xlsx")

D_GC_Genes <- Candidate_genes %>% filter(`Log2FC_D.MS` > log(1.5,2) & `p_D.MS` < 0.05 & `Count_D_GC` > 100) %>% pull(Gene_Symbol)
S_GC_Genes <- Candidate_genes %>% filter(`Log2FC_S.DM` > log(1.5,2) & `p_S.DM` < 0.05 & `Count_S_GC` > 100) %>% pull(Gene_Symbol)

# Cluster comparisons
granule_score <- AddModuleScore(object = granule, features = list(D_GC_Genes), name = "D_GC_score")
colnames(granule_score@meta.data)[ncol(granule_score@meta.data)] <- "D-GC Module"
granule_score <- AddModuleScore(object = granule_score, features = list(S_GC_Genes), name = "S_GC_score")
colnames(granule_score@meta.data)[ncol(granule_score@meta.data)] <- "S-GC Module"
granule_score@meta.data

D_GC_missing_genes <- setdiff(D_GC_Genes, rownames(granule))
S_GC_missing_genes <- setdiff(S_GC_Genes, rownames(granule))

missing_genes <- list(D_GC_missing_genes, S_GC_missing_genes)
max_length <- max(lengths(missing_genes))
n <- max(length(D_GC_missing_genes), length(S_GC_missing_genes))
missing_genes <- data.frame(D_GC = D_GC_missing_genes[seq_len(n)],
                            S_GC = S_GC_missing_genes[seq_len(n)])

D_GC_Feature <- FeaturePlot(granule_score, features = "D-GC Module", cols = c("#E8E9EB", "#FF0000"), raster = TRUE, label = FALSE,
                            raster.dpi = c(500, 500), pt.size = 1, label.size = 2.5, min.cutoff = "q05", max.cutoff = "q95") +
  guides(color = guide_colorbar(title = "Module score")) +
  theme(plot.title = element_blank(), legend.position = c(0.03, 0.85),
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.line = element_line(linewidth = 0.1, colour = "black"),
        axis.ticks = element_line(linewidth = 0.1)) & NoAxes()

D_GC_Feature <- FeaturePlot(granule_score, features = "D-GC Module", cols = c("#E8E9EB", "#FF0000"), raster = TRUE, label = FALSE,
                            raster.dpi = c(500, 500), pt.size = 1, label.size = 2.5, min.cutoff = "q05", max.cutoff = "q95") +
  theme(plot.title = element_blank(), legend.position = "none") & NoAxes()

D_GC_Feature

ggsave("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/261103, D_GC_features.tiff", D_GC_Feature, width = 5, height = 5, dpi = 600, compression = "lzw")
ggsave("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/261103, D_GC_features.jpeg", D_GC_Feature, width = 5, height = 5, dpi = 600, quality = 100)

S_GC_Feature <- FeaturePlot(granule_score, features = "S-GC Module", cols = c("#E8E9EB", "#008000"), raster = TRUE, label = FALSE,
                            raster.dpi = c(500, 500), pt.size = 1, label.size = 2.5, min.cutoff = "q05", max.cutoff = "q95") +
  guides(color = guide_colorbar(title = "Module score")) +
  theme(plot.title = element_blank(), legend.position = c(0.03, 0.85),
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        axis.line = element_line(linewidth = 0.1, colour = "black"),
        axis.ticks = element_line(linewidth = 0.1)) & NoAxes()

S_GC_Feature <- FeaturePlot(granule_score, features = "S-GC Module", cols = c("#E8E9EB", "#008000"), raster = TRUE, label = FALSE,
                            raster.dpi = c(500, 500), pt.size = 1, label.size = 2.5, min.cutoff = "q05", max.cutoff = "q95") +
  theme(plot.title = element_blank(), legend.position = "none") & NoAxes()

S_GC_Feature

ggsave("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/261103, S_GC_features.tiff", S_GC_Feature, width = 5, height = 5, dpi = 600, compression = "lzw")
ggsave("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/261103, S_GC_features.jpeg", S_GC_Feature, width = 5, height = 5, dpi = 600, quality = 100)

Idents(granule) <- "subcluster"
D_GC_high_genes = c("Fgf1", "Ryr3", "Igfbp5", "Kcnip1", "Vipr2", "Ppm1e", "Kcnq3", "Scn1a")
D_GC_DotPlot <- DotPlot(granule, features = D_GC_high_genes, group.by = "subcluster",
                        dot.scale = 5, scale.min = 0.1) +
  guides(size = guide_legend(title = "Expression\n(%)", order = 1), color = guide_colorbar(title = "Scaled\nexpression", order = 2)) +
  theme(plot.title = element_blank(), text = element_text(size = 8), axis.text = element_text(size = 8, face = "italic"),
        axis.title = element_blank(), legend.position = "right", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3, colour = "black"),
        axis.ticks = element_line(linewidth = 0.3)) + RotatedAxis()
D_GC_DotPlot

D_GC_FeaturePlot <- FeaturePlot(granule, features = D-GC_high_genes, cols = c("#E8E9EB", "#FF0000"), raster = TRUE, 
                                raster.dpi = c(500, 500), pt.size = 1, order = TRUE, min.cutoff = "q05", max.cutoff = "q95", ncol = 4) &
  theme(plot.title = element_text(size = 12, face = "italic"),
        legend.position = c(0.03, 0.85), legend.key.size = unit(0.15, "cm"),
        legend.text = element_text(size = 4, margin = margin(l = 0.5, unit = "pt")), legend.title = element_blank()) & NoAxes()
D_GC_FeaturePlot


for (g in D_GC_high_genes) {
  p <- FeaturePlot(granule, features = g, cols = c("#E8E9EB", "#FF0000"),
                   raster = TRUE, raster.dpi = c(500, 500), pt.size = 1, order = TRUE,
                   min.cutoff = "q05", max.cutoff = "q95") &
    theme(plot.title = element_blank(), aspect.ratio = 1, legend.position = "none") & NoAxes()
  ggsave(paste0("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/", "261103, ", g, ".tiff"), p, width = 2.5, height = 2.5, dpi = 600, compression = "lzw")
  }

for (g in D_GC_high_genes) {
  p <- FeaturePlot(granule, features = g, cols = c("#E8E9EB", "#FF0000"),
                   raster = TRUE, raster.dpi = c(500, 500), pt.size = 1, order = TRUE,
                   min.cutoff = "q05", max.cutoff = "q95") &
    theme(plot.title = element_blank(), aspect.ratio = 1, legend.position = "none") & NoAxes()
  ggsave(paste0("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/", "261103, ", g, ".jpeg"), p, width = 2.5, height = 2.5, dpi = 600, quality = 100)
}


# S_GC specific genes
S_GC_high_genes = c("Atp2b4", "Marcks", "Ptpre", "Synpr", "Igsf3", "Necab3", "Fxyd6", "Erbb4")
S_GC_DotPlot <- DotPlot(granule, features = S_GC_high_genes, group.by = "subcluster",
                        dot.scale = 5, scale.min = 0.1) +
  guides(size = guide_legend(title = "Expression\n(%)", order = 1), color = guide_colorbar(title = "Scaled\nexpression", order = 2)) +
  theme(plot.title = element_blank(), text = element_text(size = 8), axis.text = element_text(size = 8, face = "italic"),
        axis.title = element_blank(), legend.position = "right", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3, colour = "black"),
        axis.ticks = element_line(linewidth = 0.3)) + RotatedAxis()
S_GC_DotPlot

S_GC_FeaturePlot <- FeaturePlot(granule, features = S_GC_high_genes, cols = c("#E8E9EB", "#008000"), raster = TRUE, 
                                raster.dpi = c(500, 500), pt.size = 1, order = TRUE, min.cutoff = "q05", max.cutoff = "q95", ncol = 4) &
  theme(plot.title = element_text(size = 12, face = "italic"),
        legend.position = c(0.03, 0.85), legend.key.size = unit(0.15, "cm"),
        legend.text = element_text(size = 4, margin = margin(l = 0.5, unit = "pt")), legend.title = element_blank()) & NoAxes()
S_GC_FeaturePlot


for (g in S_GC_high_genes) {
  p <- FeaturePlot(granule, features = g, cols = c("#E8E9EB", "#008000"),
                   raster = TRUE, raster.dpi = c(500, 500), pt.size = 1, order = TRUE,
                   min.cutoff = "q05", max.cutoff = "q95") &
    theme(plot.title = element_blank(), aspect.ratio = 1, legend.position = "none") & NoAxes()
  ggsave(paste0("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/", "261103, ", g, ".tiff"), p, width = 2.5, height = 2.5, dpi = 600, compression = "lzw")
}

for (g in S_GC_high_genes) {
  p <- FeaturePlot(granule, features = g, cols = c("#E8E9EB", "#008000"),
                   raster = TRUE, raster.dpi = c(500, 500), pt.size = 1, order = TRUE,
                   min.cutoff = "q05", max.cutoff = "q95") &
    theme(plot.title = element_blank(), aspect.ratio = 1, legend.position = "none") & NoAxes()
  ggsave(paste0("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/", "261103, ", g, ".jpeg"), p, width = 2.5, height = 2.5, dpi = 600, quality = 100)
}

# Save plot
fig2j_k <- read_pptx()
fig2j_k <- add_slide(fig2j_k, layout = "Blank", master = "Office Theme")
fig2j_k <- ph_with(fig2j_k, value = dml(ggobj = UMAP_cerebellum), location = ph_location(width = 2, height = 2.5, left = 0, top = 0))
fig2j_k <- ph_with(fig2j_k, value = dml(ggobj = D_GC_Feature), location = ph_location(width = 2, height = 2.5, left = 2.5, top = 0))
fig2j_k <- ph_with(fig2j_k, value = dml(ggobj = S_GC_Feature), location = ph_location(width = 2, height = 2.5, left = 5, top = 0))
fig2j_k <- add_slide(fig2j_k, layout = "Blank", master = "Office Theme")
fig2j_k <- ph_with(fig2j_k, value = dml(ggobj = D_GC_FeaturePlot), location = ph_location(width = 3.8, height = 2.8, left = 0, top = 0))
fig2j_k <- ph_with(fig2j_k, value = dml(ggobj = D_GC_DotPlot), location = ph_location(width = 3.5, height = 2, left = 0, top = 3))
fig2j_k <- ph_with(fig2j_k, value = dml(ggobj = S_GC_FeaturePlot), location = ph_location(width = 3.8, height = 2.8, left = 3.8, top = 0))
fig2j_k <- ph_with(fig2j_k, value = dml(ggobj = S_GC_DotPlot), location = ph_location(width = 3.5, height = 2, left = 3.8, top = 3))
print(fig2j_k, target = "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Figure 7,s7, 261102, GC scRNA-seq subclustering.pptx")

write_xlsx(missing_genes, "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Figure 7,s7, 261102, GC scRNA-seq subclustering missing genes.xlsx")
