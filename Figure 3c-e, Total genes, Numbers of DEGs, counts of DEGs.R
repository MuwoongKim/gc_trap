library(readxl)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(writexl)
library(officer)
library(rvg)

# Opening files, selection of genes
TRAP <- read_excel("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Figure 2, D-GC vs M-GC vs S-GC.xlsx")
TRAP <- TRAP[!is.na(TRAP$Gene_Symbol) & TRAP$Gene_Symbol != "", , drop = FALSE]
FC_criteria <- abs(TRAP$Log2FC_S.D) > log2(1.5) | abs(TRAP$Log2FC_S.M) > log2(1.5) | abs(TRAP$Log2FC_M.D) > log2(1.5)
count_criteria <- TRAP$Count_D_GC > 100 | TRAP$Count_M_GC > 100 | TRAP$Count_S_GC > 100
pvalue_criteria <- TRAP$p_S.D < 0.05 | TRAP$p_S.M < 0.05 | TRAP$p_M.D < 0.05
selected_genes <- list(unique(TRAP$Gene_Symbol[which(FC_criteria)]), unique(TRAP$Gene_Symbol[which(count_criteria)]),
                       unique(TRAP$Gene_Symbol[which(pvalue_criteria)]))
names(selected_genes) <- c("FC > |1.5|", "Count > 100", "p value < 0.05")
TRAP_DEGs <- TRAP %>% filter(FC_criteria & count_criteria & pvalue_criteria) %>% arrange(p_S.D)

# Figure 3C, Selection of DEGs (Venn-Diagram)
Total_gene <- paste0("Detected genes (", prettyNum(length(unique(TRAP$Gene_Symbol)), big.mark = ","), ")")
TRAP[duplicated(TRAP$Gene_Symbol), ]                              # Confirmation
TRAP %>% filter(is.na(Gene_Symbol) | trimws(Gene_Symbol) == "")   # Confirmation
Venn_data <- process_data(Venn(selected_genes))
Circle <- venn_setedge(Venn_data)
Area   <- venn_regionedge(Venn_data)
Region <- venn_regionlabel(Venn_data)

Region$label <- paste0(prettyNum(Region$count, big.mark = ","), "\n", "(", sprintf("%.1f", 100 * Region$count / length(unique(TRAP$Gene_Symbol))), "%)")
fill_colors <- c("1" = "#BF00FF", "2" = "#FFBF00", "3" = "#00FFBF", "1/2" = "#FFBFFF", "1/3" = "#BFFFFF",
                 "2/3" = "#FFFFBF","1/2/3" = "#FFFFFF")

Gene_selection <- ggplot() +
  geom_polygon(data = Area, aes(X, Y, group = id, fill = id)) +
  geom_polygon(data = Circle, aes(X, Y, group = id), color = "black", fill = NA, linewidth = 0.3) +
  geom_text(data = Region, aes(X, Y, label = label, size = ifelse(id == "1/2/3", 10 / .pt, 8 / .pt),
                               fontface = ifelse(id == "1/2/3", "bold", "plain")), lineheight = 1) +
  scale_size_identity() +
  annotate("text", x = c(-2, 5.5), y = 4, label = names(selected_genes)[1:2],
           color = c("#BF00FF", "#B08200"), vjust = -0.4, size = 10 / .pt, fontface = "bold") +
  annotate("text", x = 2, y = -8, label = names(selected_genes)[3],
           color = "#00B888", vjust = 1.3, size = 10 / .pt, fontface = "bold") +
  annotate("text", x = -5.5, y = -1.5, label = "DEGs",
           color = "black", vjust = -0.4, size = 12 / .pt, fontface = "bold") +
  annotate("segment", x = 1, y = -1.3, xend = -4, yend = -1, color = "black", linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = c("1" = "#BF00FF", "2" = "#FFBF00", "3" = "#00FFBF", "1/2" = "#FFBFFF", "1/3" = "#BFFFFF",
                 "2/3" = "#FFFFBF","1/2/3" = "#FFFFFF")) +
  scale_color_manual(values = "black") +
  coord_equal(xlim = c(-4, 8), ylim = c(-8, 4), expand = FALSE, clip = "off") +
  labs(title = Total_gene) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 5, unit = "mm")),
        plot.margin = margin(t = 5, r = 0, b = 5, l = 25, unit = "mm"))
Gene_selection

# Figure 3D, Selection of DEGs (Lollipop plot)
TRAP_counts_Combined <- tibble(x = 1:6,
    Type = factor(c("D-GCs", "S-GCs", "D-GCs", "M-GCs", "M-GCs", "S-GCs"), levels = c("D-GCs", "M-GCs", "S-GCs")),
    Count = c(sum(TRAP_DEGs$Log2FC_S.D < -log2(1.5) & TRAP_DEGs$p_S.D < 0.05, na.rm = TRUE),
              sum(TRAP_DEGs$Log2FC_S.D >  log2(1.5) & TRAP_DEGs$p_S.D < 0.05, na.rm = TRUE),
              sum(TRAP_DEGs$Log2FC_M.D < -log2(1.5) & TRAP_DEGs$p_M.D < 0.05, na.rm = TRUE),
              sum(TRAP_DEGs$Log2FC_M.D >  log2(1.5) & TRAP_DEGs$p_M.D < 0.05, na.rm = TRUE),
              sum(TRAP_DEGs$Log2FC_S.M < -log2(1.5) & TRAP_DEGs$p_S.M < 0.05, na.rm = TRUE),
              sum(TRAP_DEGs$Log2FC_S.M >  log2(1.5) & TRAP_DEGs$p_S.M < 0.05, na.rm = TRUE)))

TRAP_counts_plot <- ggplot(TRAP_counts_Combined, aes(x, Count, color = Type)) +
  geom_segment(aes(xend = x, y = 0, yend = Count), linewidth = 0.7) +
  geom_point(size = 2.1) +
  scale_color_manual(values = c("D-GCs" = "red", "M-GCs" = "blue", "S-GCs" = "#008000"),
                     labels = c("D-GCs" = "High in D-GCs", "M-GCs" = "High in M-GCs", "S-GCs" = "High in S-GCs")) +
  geom_vline(xintercept = c(2.5, 4.5), linetype = "dashed", color = "black", linewidth = 0.3) +
  scale_x_continuous(breaks = c(1.5, 3.5, 5.5), labels = c("D-GCs\nvs\nS-GCs", "D-GCs\nvs\nM-GCs", "M-GCs\nvs\nS-GCs"),
                     limits = c(0.5, 6.5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 450)) +
  labs(x = NULL, y = paste0("Number of", "\n", "up-regulated DEGs"), color = NULL) +
  theme_classic() +
  theme(aspect.ratio = 0.5, axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, color = "black", hjust = 0.5, lineheight = 0.8),
        axis.text.y = element_text(size = 8),
        legend.position = "inside", legend.position.inside = c(0.98, 1.00), legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = NA), legend.direction = "vertical",
        legend.text = element_text(size = 8), legend.key.height = grid::unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(0, 1, 3, 3))
TRAP_counts_plot

# Figure 3E, Selection of DEGs (Venn-Diagram)
DEGs_count <- list(`D-GCs vs S-GCs` = TRAP_DEGs$Gene_Symbol[which(abs(TRAP_DEGs$Log2FC_S.D) > log2(1.5) & TRAP_DEGs$p_S.D < 0.05)],
                   `D-GCs vs M-GCs` = TRAP_DEGs$Gene_Symbol[which(abs(TRAP_DEGs$Log2FC_M.D) > log2(1.5) & TRAP_DEGs$p_M.D < 0.05)],
                   `M-GCs vs S-GCs` = TRAP_DEGs$Gene_Symbol[which(abs(TRAP_DEGs$Log2FC_S.M) > log2(1.5) & TRAP_DEGs$p_S.M < 0.05)])

TRAP_VennDiagram <- ggVennDiagram(DEGs_count, category.names = c("", "", ""),
                                  label = "both", label_geom = "text", label_alpha = 1, label_size = 7 / .pt,
                                  label_color = ifelse(venn_regionlabel(process_data(Venn(DEGs_count)))$count >= 300, "white", "black"),
                                  edge_size = 0.3, label_percent_digit = 1, edge_lty = "solid") +
  scale_fill_gradient(low = "white", high = "#363737") +
  annotate("text", x = c(-5.7, 9.7), y = -3, label = c("D-GCs\nvs\nS-GCs", "D-GCs\nvs\nM-GCs"), vjust = -0.4, size = 8 / .pt) +
  annotate("text", x = 2, y = -8, label = "M-GCs vs S-GCs", vjust = 1.3, size = 8 / .pt) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  coord_equal(clip = "off") +
  theme(panel.grid = element_blank(), panel.background = element_blank(), plot.background = element_blank(),
        legend.text = element_text(size = 6, margin = margin(l = 0.5, unit = "mm")),
        legend.title = element_text(size = 8, face = "bold", margin = margin(b = 0.5, unit = "mm")),
        legend.key.size = grid::unit(0.2, "cm"), legend.position = "inside", legend.position.inside = c(0.85, 0.3))

TRAP_VennDiagram

# Export data to excel
D_vs_M <- DEGs_count[["D-GCs vs M-GCs"]]
S_vs_D <- DEGs_count[["D-GCs vs S-GCs"]]
S_vs_M <- DEGs_count[["M-GCs vs S-GCs"]]

max_len <- max(length(D_vs_M), length(S_vs_D), length(S_vs_M))

length(D_vs_M) <- max_len
length(S_vs_D) <- max_len
length(S_vs_M) <- max_len

TRAP_Comparisons_Fig3D <- data.frame(D_vs_M, S_vs_D, S_vs_M)

print(Gene_selection)
print(TRAP_counts_plot)
print(TRAP_VennDiagram)

write_xlsx(TRAP_DEGs, "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Supplementary Data 4, List of DEGs.xlsx")
write_xlsx(TRAP_Comparisons_Fig3D, "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Figure 3c-e, DEGs used in Venndiagram, Lollipop plot.xlsx")

# Export plots to powerpoint
Figure_3c_e <- read_pptx()
Figure_3c_e <- add_slide(Figure_3c_e, layout = "Blank", master = "Office Theme")
Figure_3c_e <- ph_with(Figure_3c_e, dml(ggobj = Gene_selection),
                       location = ph_location(left = 2, top = 1, width = 3, height = 3))
Figure_3c_e <- ph_with(Figure_3c_e, dml(ggobj = TRAP_counts_plot),
                       location = ph_location(left = 5.2, top = 0.8, width = 2, height = 2))
Figure_3c_e <- ph_with(Figure_3c_e, dml(ggobj = TRAP_VennDiagram),
                       location = ph_location(left = 5.2, top = 2, width = 2.2, height = 2.2))
print(Figure_3c_e, target = "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Figure 3c-e, Venndiagram, Lollipop plot.pptx")
