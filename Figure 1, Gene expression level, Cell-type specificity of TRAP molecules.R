library(tidyverse)
library(ggpubr)
library(broom)
library(ggh4x)
library(DescTools)
library(readxl)
library(writexl)
library(officer)
library(rvg)

## Figure 1D - TRAP yield (GAPDH expression)
# 1. Sample information
TRAP_GAPDH_data <- read_excel("C:/Users/mw/Desktop/mw/Figure/Summarized Figures/Datasets/Figure 1D, 260402, TRAP molecule expression_modified.xlsx") %>% pivot_longer(everything(), names_to = "Sample", values_to = "Data") %>%
  filter(Sample %in% GAPDH_sample_levels, !is.na(Data), Data > 0) %>%
  mutate(Data = log2(Data), Sample = factor(Sample, levels = GAPDH_sample_levels), x.position = GAPDH_x[as.integer(Sample)])

GAPDH_sample_levels <- c("GABRa6_L10a", "sSyn_L10a", "GABRa6_L22", "CAG_L22", "GABRa6_DIO_rev_L22(Atoh1_Cre)", "sSyn_DIO_rev_L22(Atoh1_Cre)", "sSyn_DIO_rev_L22(Pcp2_Cre)", "GfaABC1D_L22")
GAPDH_x <- c(1, 2, 3, 4, 6.5, 7.5, 8.5, 9.5)
GAPDH_colors <- setNames(c("#948A54", "#706827", "#947A6B", "#595959", "#633503", "#7030A0", "#E46C0A", "#DAB200"), c("GABRa6_L10a", "sSyn_L10a", "GABRa6_L22", "CAG_L22", "GABRa6_DIO_rev_L22(Atoh1_Cre)", "sSyn_DIO_rev_L22(Atoh1_Cre)", "sSyn_DIO_rev_L22(Pcp2_Cre)", "GfaABC1D_L22"))

# 2. Divide into two groups
GAPDH_i_iv <- TRAP_GAPDH_data %>% filter(Sample %in% c("GABRa6_L10a", "sSyn_L10a", "GABRa6_L22", "CAG_L22")) %>% droplevels()
GAPDH_v_viii <- TRAP_GAPDH_data %>% filter(Sample %in% c("GABRa6_DIO_rev_L22(Atoh1_Cre)", "sSyn_DIO_rev_L22(Atoh1_Cre)",
                                                         "sSyn_DIO_rev_L22(Pcp2_Cre)", "GfaABC1D_L22")) %>% droplevels()
# 3. Statistics
model_i_iv <- aov(Data ~ Sample, data = GAPDH_i_iv)
anova_model_i_iv <- anova(model_i_iv)
model_v_viii <- aov(Data ~ Sample, data = GAPDH_v_viii)
anova_model_v_viii <- anova(model_v_viii)

GAPDH_anova <- tibble(Block = c("i-iv", "v-viii"), Effect = "Sample",
                      DFn = c(anova_model_i_iv$Df[1], anova_model_v_viii$Df[1]),
                      DFd = c(anova_model_i_iv$Df[2], anova_model_v_viii$Df[2]),
                      SSn = c(anova_model_i_iv$`Sum Sq`[1], anova_model_v_viii$`Sum Sq`[1]),
                      SSd = c(anova_model_i_iv$`Sum Sq`[2], anova_model_v_viii$`Sum Sq`[2]),
                      F = c(anova_model_i_iv$`F value`[1], anova_model_v_viii$`F value`[1]),
                      p = c(anova_model_i_iv$`Pr(>F)`[1], anova_model_v_viii$`Pr(>F)`[1])) %>%
  mutate(eta2 = SSn / (SSn + SSd), omega2 = DFn * (F - 1) / (DFn * F + DFd + 1))

Scheffe_i_iv <- ScheffeTest(model_i_iv, conf.level = 0.95)$Sample
Scheffe_v_viii <- ScheffeTest(model_v_viii, conf.level = 0.95)$Sample

GAPDH_mean_sd <- TRAP_GAPDH_data %>% group_by(Sample) %>% summarise(mean = mean(Data), sd = sd(Data), .groups = "drop")
GAPDH_n <- TRAP_GAPDH_data %>% dplyr::count(Sample, name = "n")
GAPDH_post_hoc <- bind_rows("i-iv" = as.data.frame(Scheffe_i_iv) %>% rownames_to_column("comparison"),
                            "v-viii" = as.data.frame(Scheffe_v_viii) %>% rownames_to_column("comparison"), .id = "Block") %>%
  separate(comparison, into = c("group2", "group1"), sep = "-") %>%
  rename(estimate = diff, conf.low = lwr.ci, conf.high = upr.ci, p.adj = pval) %>%
  mutate(method = "Scheffe", p.adj.signif = as.character(cut(p.adj, c(-Inf, 0.001, 0.01, 0.05, Inf),
                                                             labels = c("***", "**", "*", "ns")))) %>%
  left_join(GAPDH_n %>% rename(group1 = Sample, n1 = n), by = "group1") %>%
  left_join(GAPDH_n %>% rename(group2 = Sample, n2 = n), by = "group2")

GAPDH_post_hoc <- GAPDH_post_hoc %>%
  left_join(GAPDH_mean_sd %>% rename(group1 = Sample, mean1 = mean, sd1 = sd), by = "group1") %>%
  left_join(GAPDH_mean_sd %>% rename(group2 = Sample, mean2 = mean, sd2 = sd), by = "group2") %>%
  mutate(Mean_SD = sprintf("%.2f ± %.2f vs %.2f ± %.2f", mean1, sd1, mean2, sd2))

GAPDH_y_lower <- floor(min(TRAP_GAPDH_data$Data, 0))
GAPDH_y_upper <- ceiling(max(TRAP_GAPDH_data$Data, 0))

GAPDH_post_hoc_sig <- GAPDH_post_hoc %>% filter(p.adj < 0.05) %>%
  mutate(xmin = GAPDH_x[match(group1, GAPDH_sample_levels)], xmax = GAPDH_x[match(group2, GAPDH_sample_levels)]) %>%
  group_by(Block) %>% arrange(abs(xmax - xmin), xmin, xmax, .by_group = TRUE) %>%
  mutate(y.position = max(12, GAPDH_y_upper + 0.5) + (row_number() - 1) * 1) %>%
  ungroup()

GAPDH_dividers <- mean(GAPDH_x[4:5])
GAPDH_footer <- tibble(xmin = GAPDH_x[c(1, 3, 5)] - 0.3, xmax = GAPDH_x[c(2, 4, 8)] + 0.3, label = c("L10a-GFP", "L22-HA", "L22-HA")) %>% mutate(x = (xmin + xmax) / 2)

# 4. Export boxplot by ggplot2
TRAP_Expression <- ggplot(TRAP_GAPDH_data, aes(x = x.position, y = Data)) +
  geom_boxplot(aes(group = Sample, color = Sample), fill = "white", linewidth = 0.3, outlier.shape = NA, width = 0.7) +
  geom_jitter(aes(color = Sample), width = 0.15, height = 0, size = 0.5, alpha = 0.6) +
  scale_color_manual(values = GAPDH_colors) +
  scale_x_continuous(breaks = GAPDH_x, labels = c("i", "ii", "iii", "iv", "v", "vi", "vii", "viii"), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = seq(GAPDH_y_lower, GAPDH_y_upper, by = 2), expand = expansion(mult = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  annotate("segment", x = GAPDH_dividers, xend = GAPDH_dividers, y = GAPDH_y_lower, yend = GAPDH_y_upper,
           linetype = "solid", linewidth = 0.3) +
  geom_segment(data = GAPDH_footer, aes(x = xmin, xend = xmax), y = GAPDH_y_lower - 1.5, yend = GAPDH_y_lower - 1.5,
               linewidth = 0.3, inherit.aes = FALSE) +
  geom_text(data = GAPDH_footer, aes(x = x, label = label), y = GAPDH_y_lower - 2.5, size = 8 / .pt, inherit.aes = FALSE) +
  coord_cartesian(xlim = c(GAPDH_x[1] - 0.6, GAPDH_x[8] + 0.6), ylim = c(GAPDH_y_lower, GAPDH_y_upper), clip = "off") +
  labs(x = NULL, y = "log2(FC normalized by control)") +
  theme_classic(base_size = 10) +
  theme(text = element_text(color = "black"), axis.title = element_text(size = 10), axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3), legend.position = "none", plot.margin = margin(t = 60, r = 5, b = 40, l = 5))

if (nrow(GAPDH_post_hoc_sig) > 0) {TRAP_Expression <- TRAP_Expression + stat_pvalue_manual(
      GAPDH_post_hoc_sig, label = "p.adj.signif", xmin = "xmin", xmax = "xmax", y.position = "y.position", step.increase = 0,
      tip.length = 0, size = 12 / .pt, bracket.size = 0.6, vjust = 0.52)}

TRAP_Expression

## Figure 1E
# 1. Sample information
TRAP_Specificity <- read.csv("C:/Users/mw/Desktop/mw/Figure/Summarized Figures/Datasets/Figure 1E, 251125, Cell_type_specificity_modified.csv", check.names = FALSE)
TRAP_Specificity_data <- tibble(Group = rep(c("Gabra6", "Pcp2", "Pvalb", "Prkcd", "Gdf10"), ncol(TRAP_Specificity) - 1),
                                Subgroup = rep(colnames(TRAP_Specificity)[-1],
                                               each = length(c("Gabra6", "Pcp2", "Pvalb", "Prkcd", "Gdf10"))),
                                Data = as.numeric(unlist(TRAP_Specificity[1:5, -1], use.names = FALSE))) %>%
  filter(!is.na(Data), Data > 0) %>%
  mutate(Data = log2(Data), Group = factor(Group, levels = c("Gabra6", "Pcp2", "Pvalb", "Prkcd", "Gdf10")),
  Subgroup = factor(Subgroup, levels = c("GAB_DIO_L22 (Atoh1)", "sSyn_DIO_L22 (Atoh1)", "sSyn_DIO_L22 (Pcp2)", "Gfab1d_L22"))) %>%
  filter(!is.na(Subgroup))

Specificity_summary <- TRAP_Specificity_data %>%
  group_by(Group, Subgroup) %>%
  summarise(n = n(), mean = mean(Data), sd = sd(Data), y.position = max(Data) -0.1, .groups = "drop")

# 2. Statistics
Specificity_subgroup_levels <- c("GAB_DIO_L22 (Atoh1)", "sSyn_DIO_L22 (Atoh1)", "sSyn_DIO_L22 (Pcp2)", "Gfab1d_L22")
Specificity_t_test <- TRAP_Specificity_data %>%
  group_by(Group, Subgroup) %>%
  group_modify(~ tidy(t.test(.x$Data, mu = 0, alternative = "two.sided"))) %>%
  ungroup() %>%
  rename(t = statistic, df = parameter, p = p.value) %>%
  left_join(Specificity_summary, by = c("Group", "Subgroup")) %>%
  group_by(Group) %>%
  mutate(p.adj = p.adjust(p, method = "holm", n = length(Specificity_subgroup_levels)),
         method = "One-sample two-sided t-test", null.value = 0, p.adjust.method = "holm",
         p.adj.signif = case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**",
                                  p.adj < 0.05 ~ "*", TRUE ~ "ns"),
         Mean_SD = sprintf("%.2f ± %.2f", mean, sd)) %>% ungroup()
Specificity_y_axis <- TRAP_Specificity_data %>%
  group_by(Group) %>%
  summarise(lower = floor(min(Data, 0)), upper = ceiling(max(Data, 0)), .groups = "drop")
Specificity_y_scales <- map2(Specificity_y_axis$lower, Specificity_y_axis$upper,
                             ~ scale_y_continuous(limits = c(.x, .y + 0.8), breaks = seq(.x, .y, 1), expand = expansion(mult = 0)))

# 3. Export boxplot by ggplot2
TRAP_specificity <- ggplot(TRAP_Specificity_data, aes(x = Subgroup, y = Data)) +
  geom_boxplot(aes(color = Subgroup), fill = "white", linewidth = 0.3, outlier.shape = NA, width = 0.7) +
  geom_jitter(aes(color = Subgroup), width = 0.15, height = 0, size = 0.5, alpha = 0.6) +
  scale_color_manual(values = setNames(c("#633503", "#7030A0", "#E46C0A", "#DAB200"), Specificity_subgroup_levels)) +
  scale_x_discrete(limits = Specificity_subgroup_levels, labels = c("v", "vi", "vii", "viii")) +
  facet_wrap(~Group, scales = "free_y", strip.position = "bottom", nrow = 1) +
  facetted_pos_scales(y = Specificity_y_scales) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_text(data = filter(Specificity_t_test, p.adj < 0.05),
            aes(x = Subgroup, y = y.position, label = p.adj.signif),
            size = 12 / .pt, vjust = 0, inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "log2(FC normalized by All)") +
  theme_classic(base_size = 10) +
  theme(text = element_text(color = "black"), axis.title = element_text(size = 10),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3), legend.position = "none",
        strip.text = element_text(size = 10, face = "italic"), strip.placement = "outside",
        strip.background = element_blank(), panel.spacing = unit(0.5, "lines"),
        plot.margin = margin(t = 50, r = 5, b = 5, l = 5))
TRAP_specificity

# Export to powerpoint or excel
Figure_1d_e <- read_pptx() %>%
  add_slide(layout = "Blank", master = "Office Theme") %>%
  ph_with(dml(ggobj = TRAP_Expression), ph_location(width = 3.3, height = 3.5, left = 3.5, top = 0)) %>%
  ph_with(dml(ggobj = TRAP_specificity), ph_location(width = 7, height = 3.5, left = 0, top = 3))
print(Figure_1d_e, target = "C:/Users/mw/Desktop/mw/Figure/Summarized Figures/Figure, Statistics/Figure 1d_e, TRAP_Yield_Specificity.pptx")

Figure_1d_e_statistics <- bind_rows(ANOVA_Fig1d = GAPDH_anova, Scheffe_Fig1d = GAPDH_post_hoc,
                                    T_test_Fig1e = Specificity_t_test, .id = "Test")
write_xlsx(Figure_1d_e_statistics, "C:/Users/mw/Desktop/mw/Figure/Summarized Figures/Figure, Statistics/Figure 1d_e, TRAP_Yield_specificity_statistics.xlsx")
