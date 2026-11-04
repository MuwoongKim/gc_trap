library(DESeq2)
library(writexl)
library(biomaRt)

TRAP_Sum <- read.delim("C:/Users/mw/Desktop/mw/Figure/Summarized figures/datasets/240903, TRAP_count.txt",
                       header = TRUE, skip = 1, check.names = FALSE)
Count_TRAP <- TRAP_Sum[, 7:ncol(TRAP_Sum)]
rownames(Count_TRAP) <- TRAP_Sum[, 1]
colnames(Count_TRAP) <- c("mRNA_GC", "mRNA_All", "Total_RNA", "S_GC_1", "S_GC_2", "S_GC_3",
                          "D_GC_1", "D_GC_2", "D_GC_3", "M_GC_1", "M_GC_2", "M_GC_3", "mRNA_PC")
Gene_ID <- rownames(Count_TRAP)
Count_all <- Count_TRAP[, c("Total_RNA", "mRNA_All", "mRNA_PC", "mRNA_GC",
                            "D_GC_1", "D_GC_2", "D_GC_3", "M_GC_1", "M_GC_2", "M_GC_3", "S_GC_1", "S_GC_2", "S_GC_3"), drop = FALSE]
Count_GC <- Count_TRAP[, c("D_GC_1", "D_GC_2", "D_GC_3", "M_GC_1", "M_GC_2", "M_GC_3", "S_GC_1", "S_GC_2", "S_GC_3"), drop = FALSE]


# Because of no internet connection to ensembl server (v112), I used downloaded gtf file instead.
#mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", version = 112)
#Gene_mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = unique(Gene_ID), mart = mart)
#colnames(Gene_mapping) <- c("Gene_ID", "Gene_Symbol")

BioMart_v112 <- rtracklayer::import("C:/Users/mw/Desktop/mw/Figure/Summarized figures/datasets/Mus_musculus.GRCm39.112.gtf.gz",
                                    format = "gtf", feature.type = "gene")
BioMart_gene_id_name <- data.frame(ensembl_gene_id = as.character(BioMart_v112$gene_id),
                                   gene_name = as.character(BioMart_v112$gene_name), stringsAsFactors = FALSE)

Gene_position <- match(Gene_ID, BioMart_gene_id_name$ensembl_gene_id)
Gene_Symbol <- BioMart_gene_id_name$gene_name[Gene_position]
Gene_mapping <- data.frame(Gene_ID = Gene_ID, Gene_Symbol = Gene_Symbol, stringsAsFactors = FALSE)

entire_sample <- data.frame(condition = factor(c("Total_RNA", "mRNA_All", "mRNA_PC", "mRNA_GC", "D_GC", "D_GC", "D_GC", "M_GC", "M_GC", "M_GC", "S_GC", "S_GC", "S_GC"), levels = c("Total_RNA", "mRNA_All", "mRNA_PC", "mRNA_GC", "D_GC", "M_GC", "S_GC")), row.names = colnames(Count_all))


# CAG/PC vs D-GCs/M-GCs/S-GCs
dds_entire <- DESeqDataSetFromMatrix(countData = Count_all, colData = entire_sample, design = ~ condition)
dds_entire <- DESeq(dds_entire)

CAG_vs_D <- results(dds_entire, contrast = c("condition", "mRNA_All", "D_GC"))
CAG_vs_M <- results(dds_entire, contrast = c("condition", "mRNA_All", "M_GC"))
CAG_vs_S <- results(dds_entire, contrast = c("condition", "mRNA_All", "S_GC"))
PC_vs_D  <- results(dds_entire, contrast = c("condition", "mRNA_PC", "D_GC"))
PC_vs_M  <- results(dds_entire, contrast = c("condition", "mRNA_PC", "M_GC"))
PC_vs_S  <- results(dds_entire, contrast = c("condition", "mRNA_PC", "S_GC"))
Norm_entire <- counts(dds_entire, normalized = TRUE)

PC_CAG_GC_entire <- data.frame(Gene_mapping,
  Log2FC_A.D = CAG_vs_D$log2FoldChange, p_A.D = CAG_vs_D$pvalue, adj_p_A.D = CAG_vs_D$padj,
  Log2FC_A.M = CAG_vs_M$log2FoldChange, p_A.M = CAG_vs_M$pvalue, adj_p_A.M = CAG_vs_M$padj,
  Log2FC_A.S = CAG_vs_S$log2FoldChange, p_A.S = CAG_vs_S$pvalue, adj_p_A.S = CAG_vs_S$padj,
  Log2FC_P.D = PC_vs_D$log2FoldChange, p_P.D = PC_vs_D$pvalue, adj_p_P.D = PC_vs_D$padj,
  Log2FC_P.M = PC_vs_M$log2FoldChange, p_P.M = PC_vs_M$pvalue, adj_p_P.M = PC_vs_M$padj,
  Log2FC_P.S = PC_vs_S$log2FoldChange, p_P.S = PC_vs_S$pvalue, adj_p_P.S = PC_vs_S$padj,
  Count_CAG = Norm_entire[, "mRNA_All"],
  Count_PC = Norm_entire[, "mRNA_PC"],
  Count_D_GC = rowMeans(Norm_entire[, entire_sample$condition == "D_GC", drop = FALSE]),
  Count_M_GC = rowMeans(Norm_entire[, entire_sample$condition == "M_GC", drop = FALSE]),
  Count_S_GC = rowMeans(Norm_entire[, entire_sample$condition == "S_GC", drop = FALSE]),
  row.names = NULL)

Count_for_summary <- Count_all[, c("mRNA_All", "mRNA_PC", "D_GC_1", "D_GC_2", "D_GC_3", "M_GC_1", "M_GC_2", "M_GC_3",
                                   "S_GC_1", "S_GC_2", "S_GC_3"), drop = FALSE]
PC_CAG_GC_expression <- PC_CAG_GC_entire[!is.na(Gene_mapping$Gene_Symbol) & Gene_mapping$Gene_Symbol != "" & rowSums(Count_for_summary) > 0, ,drop = FALSE]

# D-GC vs M-GC vs S-GC
GC_comparison <- data.frame(condition = factor(c("D_GC", "D_GC", "D_GC", "M_GC", "M_GC", "M_GC", "S_GC", "S_GC", "S_GC"), levels = c("D_GC", "M_GC", "S_GC")), row.names = colnames(Count_GC))
dds_GC <- DESeqDataSetFromMatrix(countData = Count_GC, colData = GC_comparison, design = ~ condition)
dds_GC <- DESeq(dds_GC)

S_vs_D <- results(dds_GC, contrast = c("condition", "S_GC", "D_GC"))
S_vs_M <- results(dds_GC, contrast = c("condition", "S_GC", "M_GC"))
M_vs_D <- results(dds_GC, contrast = c("condition", "M_GC", "D_GC"))
Norm_GC <- counts(dds_GC, normalized = TRUE)

GC_mean_counts <- data.frame(Count_D_GC = rowMeans(Norm_GC[, GC_comparison$condition == "D_GC", drop = FALSE]),
                             Count_M_GC = rowMeans(Norm_GC[, GC_comparison$condition == "M_GC", drop = FALSE]),
  Count_S_GC = rowMeans(Norm_GC[, GC_comparison$condition == "S_GC", drop = FALSE]),
  row.names = NULL)

GC_pairwise_all <- data.frame(Gene_mapping, Log2FC_S.D = S_vs_D$log2FoldChange, p_S.D = S_vs_D$pvalue, adj_p_S.D = S_vs_D$padj,
                              Log2FC_S.M = S_vs_M$log2FoldChange, p_S.M = S_vs_M$pvalue, adj_p_S.M = S_vs_M$padj,
                              Log2FC_M.D = M_vs_D$log2FoldChange, p_M.D = M_vs_D$pvalue, adj_p_M.D = M_vs_D$padj,
                              GC_mean_counts, row.names = NULL)

GC_pairwise_summary <- GC_pairwise_all[!is.na(Gene_mapping$Gene_Symbol) & Gene_mapping$Gene_Symbol != "" & (GC_pairwise_all$Count_D_GC > 0 | GC_pairwise_all$Count_M_GC > 0 | GC_pairwise_all$Count_S_GC > 0), , drop = FALSE]

# D-GC/M-GC/S-GC vs Other GCs
D_vs_MS <- data.frame(condition = factor(c("D_GC","D_GC","D_GC","MS_GC","MS_GC","MS_GC","MS_GC","MS_GC","MS_GC"), levels = c("MS_GC", "D_GC")), row.names = colnames(Count_GC))
dds_D_MS <- DESeqDataSetFromMatrix(countData = Count_GC, colData = D_vs_MS, design = ~ condition)
dds_D_MS <- DESeq(dds_D_MS)
D_vs_MS <- results(dds_D_MS, contrast = c("condition", "D_GC", "MS_GC"))

M_vs_DS <- data.frame(condition = factor(c("DS_GC", "DS_GC", "DS_GC", "M_GC", "M_GC", "M_GC", "DS_GC", "DS_GC", "DS_GC"), levels = c("DS_GC", "M_GC")), row.names = colnames(Count_GC))
dds_M_DS <- DESeqDataSetFromMatrix(countData = Count_GC, colData = M_vs_DS, design = ~ condition)
dds_M_DS <- DESeq(dds_M_DS)
M_vs_DS <- results(dds_M_DS, contrast = c("condition", "M_GC", "DS_GC"))

S_vs_DM <- data.frame(condition = factor(c("DM_GC", "DM_GC", "DM_GC", "DM_GC", "DM_GC", "DM_GC", "S_GC", "S_GC", "S_GC"), levels = c("DM_GC", "S_GC")), row.names = colnames(Count_GC))
dds_S_DM <- DESeqDataSetFromMatrix(countData = Count_GC, colData = S_vs_DM, design = ~ condition)
dds_S_DM <- DESeq(dds_S_DM)
S_vs_DM <- results(dds_S_DM, contrast = c("condition", "S_GC", "DM_GC"))

GC_vs_others_all <- data.frame(Gene_mapping,
  Log2FC_D.MS = D_vs_MS$log2FoldChange, p_D.MS = D_vs_MS$pvalue, adj_p_D.MS = D_vs_MS$padj,
  Log2FC_M.DS = M_vs_DS$log2FoldChange, p_M.DS = M_vs_DS$pvalue, adj_p_M.DS = M_vs_DS$padj,
  Log2FC_S.DM = S_vs_DM$log2FoldChange, p_S.DM = S_vs_DM$pvalue, adj_p_S.DM = S_vs_DM$padj,
  GC_mean_counts, row.names = NULL)

GC_vs_others_summary <- GC_vs_others_all[!is.na(Gene_mapping$Gene_Symbol) & Gene_mapping$Gene_Symbol != "" & (GC_vs_others_all$Count_D_GC > 0 | GC_vs_others_all$Count_M_GC > 0 | GC_vs_others_all$Count_S_GC > 0), , drop = FALSE]

write_xlsx(list(Summary = PC_CAG_GC_expression, All_genes = PC_CAG_GC_entire),
           "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Figure 2, All vs GC or PC vs GC.xlsx")
write_xlsx(list(Summary = GC_pairwise_summary, All_genes = GC_pairwise_all),
           "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Figure 2, D-GC vs M-GC vs S-GC.xlsx")
write_xlsx(list(Summary = GC_vs_others_summary, All_genes = GC_vs_others_all),
           "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Figure 2, GC vs other GCs.xlsx")