library(readxl)
library(openxlsx)
library(officer)
library(rvg)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(Orthology.eg.db)

# Figure 4d, GO/KEGG analysis

# Reading the data from trend analysis
TRAP_selected <- read_excel("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Supplement Data 4, List of DEGs, pattern analysis.xlsx")
id_map <- bitr(TRAP_selected$Gene_Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
TRAP_selected <- TRAP_selected %>% left_join(id_map, by = c("Gene_Symbol" = "SYMBOL"))
group_levels_combined <- c("D-GC High", "S-GC High")
TRAP_selected_combined <- TRAP_selected %>%
  mutate(Group = case_when(Class %in% c("Continuous Decrease", "D-GC High", "S-GC Low") ~ "D-GC High",
                           Class %in% c("Continuous Increase", "S-GC High", "D-GC Low") ~ "S-GC High")) %>%
  filter(!is.na(Group), !is.na(ENTREZID))

D_genes <- TRAP_selected_combined %>% filter(Group == "D-GC High") %>% pull(ENTREZID)
S_genes <- TRAP_selected_combined %>% filter(Group == "S-GC High") %>% pull(ENTREZID)
D_genes[duplicated(D_genes)] # confirmation
S_genes[duplicated(S_genes)] # Confirmation
intersect(D_genes, S_genes) # Confirmation

# D-GC High: GO BP/CC/MF & KEGG
GO_D_BP <- enrichGO(gene = D_genes, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                    pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE)
GO_D_BP <- clusterProfiler::simplify(GO_D_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
D_BP_df <- GO_D_BP@result %>% mutate(Source = "BP", Group = "D-GC High") %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, Source, Group)

GO_D_CC <- enrichGO(gene = D_genes, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH",
                    pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE)
GO_D_CC <- clusterProfiler::simplify(GO_D_CC, cutoff = 0.7, by = "p.adjust", select_fun = min)
D_CC_df <- GO_D_CC@result %>% mutate(Source = "CC", Group = "D-GC High") %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, Source, Group)

GO_D_MF <- enrichGO(gene = D_genes, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH",
                    pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE)
GO_D_MF <- clusterProfiler::simplify(GO_D_MF, cutoff = 0.7, by = "p.adjust", select_fun = min)
D_MF_df <- GO_D_MF@result %>% mutate(Source = "MF", Group = "D-GC High") %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, Source, Group)

KEGG_D <- enrichKEGG(gene = D_genes, organism = "mmu", pAdjustMethod = "BH", pvalueCutoff = 0.1, qvalueCutoff = 0.2)
KEGG_D <- setReadable(KEGG_D, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
D_KEGG_df <- KEGG_D@result %>% mutate(Source = "KEGG", Group = "D-GC High") %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, Source, Group)

GO_D_results <- bind_rows(D_BP_df, D_CC_df, D_MF_df, D_KEGG_df)

# S-GC High: GO BP/CC/MF & KEGG
GO_S_BP <- enrichGO(gene = S_genes, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                    pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE)
GO_S_BP <- clusterProfiler::simplify(GO_S_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
S_BP_df <- GO_S_BP@result %>% mutate(Source = "BP", Group = "S-GC High") %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, Source, Group)

GO_S_CC <- enrichGO(gene = S_genes, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH",
                    pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE)
GO_S_CC <- clusterProfiler::simplify(GO_S_CC, cutoff = 0.7, by = "p.adjust", select_fun = min)
S_CC_df <- GO_S_CC@result %>% mutate(Source = "CC", Group = "S-GC High") %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, Source, Group)

GO_S_MF <- enrichGO(gene = S_genes, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH",
                    pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE)
GO_S_MF <- clusterProfiler::simplify(GO_S_MF, cutoff = 0.7, by = "p.adjust", select_fun = min)
S_MF_df <- GO_S_MF@result %>% mutate(Source = "MF", Group = "S-GC High") %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, Source, Group)

KEGG_S <- enrichKEGG(gene = S_genes, organism = "mmu", pAdjustMethod = "BH", pvalueCutoff = 0.1, qvalueCutoff = 0.2)
KEGG_S <- setReadable(KEGG_S, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
S_KEGG_df <- KEGG_S@result %>% mutate(Source = "KEGG", Group = "S-GC High") %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, Source, Group)

GO_S_results <- bind_rows(S_BP_df, S_CC_df, S_MF_df, S_KEGG_df)
GO_combined_results <- bind_rows(GO_D_results, GO_S_results)

GO_combined_results_selected <- GO_combined_results %>% filter(p.adjust < 0.05)
GO_combined_results_keyword_selected <- GO_combined_results_selected %>%
  filter(str_detect(Description, regex("synap|neur|axon|potential|plasticity|vesicle|dendrit|channel|glutam|gaba|sodium|calcium|inhibitory|excitatory|LTP|LTD|Long-term|pruning", ignore_case = TRUE)))



# Dotplot for GO/KEGG analysis
GO_combined_results_keyword_selected <- GO_combined_results_keyword_selected %>% filter(p.adjust < 0.01) %>%
    separate(GeneRatio, into = c("A", "B"), sep = "/", remove = FALSE) %>%
    mutate(GeneRatioNum = as.numeric(A) / as.numeric(B),
           log_p_adj = -log10(p.adjust),
           Group = factor(Group, levels = c("D-GC High", "S-GC High")),
           Source = factor(Source, levels = c("BP", "CC", "MF", "KEGG")),
           Description = reorder(Description, p.adjust)) %>% dplyr::select(-A, -B)
  
GO_KEGG_Combined_DotPlot <- ggplot(GO_combined_results_keyword_selected, aes(x = Group, y = Description)) +
    geom_point(aes(size = GeneRatioNum, color = log_p_adj)) +
    facet_grid(Source ~ ., scales = "free", space = "free") +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(FDR)") +
    scale_size_continuous(range = c(1, 5), name = "Gene Ratio") +
    scale_x_discrete(limits = c("D-GC High", "S-GC High")) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(color = "black", size = 8, angle = 30, hjust = 1),
          axis.text.y = element_text(size = 8, color = "black"),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          strip.text = element_text(face = "bold", size = 8, color = "white"),
          strip.background = element_rect(fill = "black"),
          panel.grid.minor = element_blank(), panel.spacing = unit(0.45, "lines"),
          legend.position = "right", legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(GO_KEGG_Combined_DotPlot)
  

# Supplement Figure 3E, SynGO analysis
syngo_raw <- read.xlsx("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Datasets/ontologies.xlsx")
# https://www.syngoportal.org/download, release 1.3, ontologies.xlsx
TERM2GENE <- syngo_raw %>% dplyr::select(id, entrez_id) %>% mutate(entrez_id = as.character(entrez_id)) %>%
  separate_rows(entrez_id, sep = ", ")
TERM2NAME <- syngo_raw %>% dplyr::select(id, name)
ortholog_map_all <- AnnotationDbi::select(Orthology.eg.db,
                                          keys = as.character(AnnotationDbi::keys(Orthology.eg.db, keytype = "Mus.musculus")),
                                          columns = "Homo.sapiens", keytype = "Mus.musculus") %>%
  dplyr::rename(ENTREZID = Mus.musculus, Human_ENTREZID = Homo.sapiens) %>%
  mutate(ENTREZID = as.character(ENTREZID), Human_ENTREZID = as.character(Human_ENTREZID))

ortholog_pairs <- ortholog_map_all %>% filter(!is.na(Human_ENTREZID)) %>% distinct(ENTREZID, Human_ENTREZID) %>%
  add_count(ENTREZID, name = "Human_matches") %>% add_count(Human_ENTREZID, name = "Mouse_matches")
ortholog_pairs %>% filter(Human_matches != 1 | Mouse_matches != 1) # confirmation

SynGO_mapping_all <- TRAP_selected_combined %>% dplyr::select(Gene_Symbol, ENTREZID, Group) %>%
  left_join(ortholog_pairs, by = "ENTREZID")

SynGO_input <- SynGO_mapping_all %>% filter(!is.na(Human_ENTREZID))
SynGO_excluded <- SynGO_mapping_all %>% filter(is.na(Human_ENTREZID))

D_genes_SynGO <- SynGO_input %>% filter(Group == "D-GC High") %>% pull(Human_ENTREZID)
S_genes_SynGO <- SynGO_input %>% filter(Group == "S-GC High") %>% pull(Human_ENTREZID)

SynGO_D <- clusterProfiler::enricher(gene = D_genes_SynGO, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                    minGSSize = 10, maxGSSize = 500, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME)
SynGO_D_results <- SynGO_D@result %>% mutate(Source = "SynGO", Group = "D-GC High") %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, Source, Group)
SynGO_S <- clusterProfiler::enricher(gene = S_genes_SynGO, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", 
                    minGSSize = 10, maxGSSize = 500, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME)
SynGO_S_results <- SynGO_S@result %>% mutate(Source = "SynGO", Group = "S-GC High") %>% 
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, Source, Group)

SynGO_all_combined <- bind_rows(SynGO_D_results, SynGO_S_results)
SynGO_combined_selected <- SynGO_all_combined %>% filter(pvalue < 0.05)


SynGO_plot_data <- SynGO_combined_selected %>%
  separate(GeneRatio, into = c("A", "B"), sep = "/", remove = FALSE) %>%
  mutate(GeneRatioNum = as.numeric(A) / as.numeric(B), log_p_value = -log10(pvalue),
         Group = factor(Group, levels = group_levels_combined), Source = factor(Source, levels = "SynGO"),
         Description = reorder(Description, -pvalue)) %>% dplyr::select(-A, -B)

SynGO_combined_DotPlot <- ggplot(SynGO_plot_data, aes(x = Group, y = Description)) +
  geom_point(aes(size = GeneRatioNum, color = log_p_value)) +
  facet_grid(Source ~ ., scales = "free", space = "free") +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p-value)") +
  scale_size_continuous(range = c(1, 5), name = "Gene Ratio") +
  scale_x_discrete(limits = group_levels_combined) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(color = "black", size = 8, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(face = "bold", size = 8, color = "white"),
    strip.background = element_rect(fill = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(0.45, "lines"),
    legend.position = "right",
    legend.key.size = grid::unit(0.3, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8))

print(SynGO_combined_DotPlot)

Figure_4d_s3e <- read_pptx()
Figure_4d_s3e <- add_slide(Figure_4d_s3e, layout = "Blank", master = "Office Theme")
Figure_4d_s3e <- ph_with(Figure_4d_s3e, value = dml(ggobj = GO_KEGG_Combined_DotPlot), location = ph_location(width = 4.5, height = 5.5, left = 0, top = 0))
Figure_4d_s3e <- add_slide(Figure_4d_s3e, layout = "Blank", master = "Office Theme")
Figure_4d_s3e <- ph_with(Figure_4d_s3e, value = dml(ggobj = SynGO_combined_DotPlot),
                     location = ph_location(width = 5.8, height = 3.2, left = 0, top = 0))
print(Figure_4d_s3e, target = "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Figure 4d,s3e, GO_KEGG_DS_keywords_dotplot.pptx")

write.xlsx(list(GO_KEGG_results = GO_combined_results_selected, SynGO_results = SynGO_combined_selected,
                SynGO_Mapping_excluded = SynGO_excluded, Mapping_all = SynGO_mapping_all),
  "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Supplement Data 5, GO_KEGG_SynGO_results.xlsx", overwrite = TRUE)