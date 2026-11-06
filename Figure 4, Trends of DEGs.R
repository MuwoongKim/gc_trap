library(readxl)
library(writexl)
library(tidyverse)
library(ggrepel)
library(ggforce)
library(shadowtext)
library(officer)
library(rvg)

TRAP <- read_excel("C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Figure 2, D-GC vs M-GC vs S-GC.xlsx")

# Classification of DEGs based on gene expression patterns
FC_criteria <- abs(TRAP$Log2FC_S.D) > log2(1.5) | abs(TRAP$Log2FC_S.M) > log2(1.5) | abs(TRAP$Log2FC_M.D) > log2(1.5)
count_criteria <- TRAP$Count_D_GC > 100 | TRAP$Count_M_GC > 100 | TRAP$Count_S_GC > 100
pvalue_criteria <- TRAP$p_S.D < 0.05 | TRAP$p_S.M < 0.05 | TRAP$p_M.D < 0.05
TRAP_DEGs <- TRAP %>% filter(FC_criteria & count_criteria & pvalue_criteria) %>% arrange(p_S.D)

# Group 1: Continuous increase
Type_1_1 <- TRAP_DEGs %>% filter(Log2FC_M.D > 0, p_M.D < 0.05, Log2FC_S.M > 0, p_S.M < 0.05) %>% mutate(Type = "1-1; D < M < S")
Type_1_2 <- TRAP_DEGs %>% filter(Log2FC_M.D > 0, p_M.D >= 0.05, Log2FC_S.M > 0, p_S.M >= 0.05, Log2FC_S.D > 0, p_S.D < 0.05) %>% mutate(Type = "1-2; D ≤ M ≤ S, D < S")

# Group 2: Continuous decrease
Type_2_1 <- TRAP_DEGs %>% filter(Log2FC_M.D < 0, p_M.D < 0.05, Log2FC_S.M < 0, p_S.M < 0.05) %>% mutate(Type = "2-1; D > M > S")
Type_2_2 <- TRAP_DEGs %>% filter(Log2FC_M.D < 0, p_M.D >= 0.05, Log2FC_S.M < 0, p_S.M >= 0.05, Log2FC_S.D < 0, p_S.D < 0.05) %>% mutate(Type = "2-2; D ≥ M ≥ S, D > S")

# ── Group 3: High in S-GCs
Type_3_1 <- TRAP_DEGs %>% filter(Log2FC_M.D > 0, p_M.D >= 0.05, Log2FC_S.M > 0, p_S.M < 0.05) %>% mutate(Type = "3-1; D ≤ M < S")
Type_3_2 <- TRAP_DEGs %>% filter(Log2FC_M.D < 0, p_M.D >= 0.05, Log2FC_S.M > 0, p_S.M < 0.05) %>% mutate(Type = "3-2; D ≥ M < S")

# ── Group 4: Low in S-GCs
Type_4_1 <- TRAP_DEGs %>% filter(Log2FC_M.D < 0, p_M.D >= 0.05, Log2FC_S.M < 0, p_S.M < 0.05) %>% mutate(Type = "4-1; D ≥ M > S")
Type_4_2 <- TRAP_DEGs %>% filter(Log2FC_M.D > 0, p_M.D >= 0.05, Log2FC_S.M < 0, p_S.M < 0.05) %>% mutate(Type = "4-2; D ≤ M > S")

# ── Group 5: Low in D-GCs
Type_5_1 <- TRAP_DEGs %>% filter(Log2FC_M.D > 0, p_M.D < 0.05, Log2FC_S.M > 0, p_S.M >= 0.05) %>% mutate(Type = "5-1; D < M ≤ S")
Type_5_2 <- TRAP_DEGs %>% filter(Log2FC_M.D > 0, p_M.D < 0.05, Log2FC_S.M < 0, p_S.M >= 0.05) %>% mutate(Type = "5-2; D < M ≥ S")

# ── Group 6: High in D-GCs
Type_6_1 <- TRAP_DEGs %>% filter(Log2FC_M.D < 0, p_M.D < 0.05, Log2FC_S.M < 0, p_S.M >= 0.05) %>% mutate(Type = "6-1; D > M ≥ S")
Type_6_2 <- TRAP_DEGs %>% filter(Log2FC_M.D < 0, p_M.D < 0.05, Log2FC_S.M > 0, p_S.M >= 0.05) %>% mutate(Type = "6-2; D > M ≤ S")

# ── Group 7, 8: High in M-GCs
Type_7 <- TRAP_DEGs %>% filter(Log2FC_M.D > 0, p_M.D < 0.05, Log2FC_S.M < 0, p_S.M < 0.05) %>% mutate(Type = "7; D < M > S")
Type_8 <- TRAP_DEGs %>% filter(Log2FC_M.D < 0, p_M.D < 0.05, Log2FC_S.M > 0, p_S.M < 0.05) %>% mutate(Type = "8; D > M < S")

Cont_inc      <- bind_rows(Type_1_1, Type_1_2) %>% mutate('Class' = "Continuous increase")
Cont_dec      <- bind_rows(Type_2_1, Type_2_2) %>% mutate('Class' = "Continuous decrease")
S_GC_high     <- bind_rows(Type_3_1, Type_3_2) %>% mutate('Class' = "High in S-GCs")
S_GC_low      <- bind_rows(Type_4_1, Type_4_2) %>% mutate('Class' = "Low in S-GCs")
D_GC_low      <- bind_rows(Type_5_1, Type_5_2) %>% mutate('Class' = "Low in D-GCs")
D_GC_high     <- bind_rows(Type_6_1, Type_6_2) %>% mutate('Class' = "High in D-GCs")
M_GC_high     <- Type_7 %>% mutate('Class' = "High in M-GCs")
M_GC_low      <- Type_8 %>% mutate('Class' = "Low in M-GCs")

Cont_inc_dec  <- bind_rows(Cont_inc,  Cont_dec)
S_GC_high_low <- bind_rows(S_GC_high, S_GC_low)
D_GC_high_low <- bind_rows(D_GC_low, D_GC_high)
M_GC_high_low <- bind_rows(M_GC_high, M_GC_low)

Total <- bind_rows(Cont_inc_dec, S_GC_high_low, D_GC_high_low, M_GC_high_low)
nrow(Total) == nrow(TRAP_DEGs) # Confirmation
sum(duplicated(Total$Gene_Symbol)) # Confirmation

# Figure 4c, Line plot for trend analysis
colors_Cont_inc_dec <- c("1-1; D < M < S" = "#2E7D32", "1-2; D ≤ M ≤ S, D < S" = "#00BFA5",
                         "2-1; D > M > S" = "#E53935", "2-2; D ≥ M ≥ S, D > S"  = "#9C27B0")
colors_S_GC_high_low <- c("3-1; D ≤ M < S" = "#66E0C2", "3-2; D ≥ M < S" = "#00695C",
                          "4-1; D ≥ M > S" = "#7A1E3A", "4-2; D ≤ M > S" = "#FF8C9E")
colors_D_GC_high_low <- c("5-1; D < M ≤ S" = "#1B5E20", "5-2; D < M ≥ S" = "#9CCC65",
                          "6-1; D > M ≥ S" = "#FF6A3D", "6-2; D > M ≤ S" = "#B11226")
colors_M_GC_high_low <- c("7; D < M > S" = "#31859C", "8; D > M < S" = "#1A237E")

colors_combined_Trends <- c(colors_Cont_inc_dec, colors_S_GC_high_low, colors_D_GC_high_low, colors_M_GC_high_low)

Trend_plot_source <- Total %>% select(Gene_Symbol, Type, Log2FC_M.D, Log2FC_S.D) %>% mutate(D_GC = 0) %>%
  pivot_longer(cols = c(D_GC, Log2FC_M.D, Log2FC_S.D), names_to = "Condition", values_to = "Log2FC") %>%
  mutate(Time = case_when(Condition == "D_GC" ~ 0, Condition == "Log2FC_M.D" ~ 1, Condition == "Log2FC_S.D" ~ 2))

Average_line_all <- Trend_plot_source %>% group_by(Time, Type) %>% summarise(Average_Log2FC = mean(Log2FC, na.rm = TRUE), .groups = "drop")

Continuous_inc_dec_trends <- Trend_plot_source %>% filter(Type %in% names(colors_Cont_inc_dec))
Continuous_inc_dec_avg <- Average_line_all %>% filter(Type %in% names(colors_Cont_inc_dec))

Continuous_S_GC_high_low_trends <- Trend_plot_source %>% filter(Type %in% names(colors_S_GC_high_low))
Continuous_S_GC_high_low_avg <- Average_line_all %>% filter(Type %in% names(colors_S_GC_high_low))

Continuous_D_GC_high_low_trends <- Trend_plot_source %>% filter(Type %in% names(colors_D_GC_high_low))
Continuous_D_GC_high_low_avg <- Average_line_all %>% filter(Type %in% names(colors_D_GC_high_low))

Continuous_M_GC_high_low_trends <- Trend_plot_source %>% filter(Type %in% names(colors_M_GC_high_low))
Continuous_M_GC_high_low_avg <- Average_line_all %>% filter(Type %in% names(colors_M_GC_high_low))

trends_theme <- theme_bw() + theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), panel.grid = element_blank(),
    legend.position = c(0.15, 0.15), legend.background = element_blank(), legend.key = element_blank(),
    axis.title = element_text(size = 10), axis.text = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 10, color = "black"), legend.title = element_blank(), legend.key.spacing.y = unit(0, "cm"),
    legend.key.height = unit(0.4, "cm"))

Continuous_inc_dec_trends_graph <- ggplot(Continuous_inc_dec_trends,
                                          aes(x = Time, y = Log2FC, group = interaction(Type, Gene_Symbol), color = Type)) +
  geom_line(alpha = 0.6) +
  geom_point(aes(shape = Type), stroke = 1, size = 2, fill = "white") +
  geom_line(data = Continuous_inc_dec_avg, aes(x = Time, y = Average_Log2FC, group = Type),
            color = "black", linewidth = 1, linetype = "solid", lineend = "round", inherit.aes = FALSE) +
  geom_line(data = Continuous_inc_dec_avg, aes(x = Time, y = Average_Log2FC, group = Type, color = Type),
            linewidth = 0.5, linetype = "dashed", lineend = "round", inherit.aes = FALSE) +
  geom_hline(yintercept = c(0, log2(1.5), -log2(1.5)),
             linewidth = 0.5, colour = "darkgray", linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 1, 2), labels = c("D-GCs", "M-GCs", "S-GCs")) +
  scale_y_continuous(breaks = seq(-2, 3, 1), limits = c(-2.3, 2.9)) +
  scale_color_manual(values = colors_Cont_inc_dec, drop = FALSE, labels = c("1-1", "1-2", "2-1", "2-2")) +
  scale_shape_manual(values = c(21, 22, 23, 24), drop = FALSE, labels = c("1-1", "1-2", "2-1", "2-2")) +
  labs(x = "GC Subtype", y = "log2(FC normalized by D-GCs)", color = "Type", shape = "Type") +
  trends_theme

Continuous_S_GC_high_low_trends_graph <- ggplot(Continuous_S_GC_high_low_trends,
                                                aes(x = Time, y = Log2FC, group = interaction(Type, Gene_Symbol), color = Type)) +
  geom_line(alpha = 0.6) +
  geom_point(aes(shape = Type), stroke = 1, size = 2, fill = "white") +
  geom_line(data = Continuous_S_GC_high_low_avg, aes(x = Time, y = Average_Log2FC, group = Type),
            color = "black", linewidth = 1, linetype = "solid", lineend = "round", inherit.aes = FALSE) +
  geom_line(data = Continuous_S_GC_high_low_avg, aes(x = Time, y = Average_Log2FC, group = Type, color = Type),
            linewidth = 0.5, linetype = "dashed", lineend = "round", inherit.aes = FALSE) +
  geom_hline(yintercept = c(0, log2(1.5), -log2(1.5)), linewidth = 0.5, colour = "darkgray", linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 1, 2), labels = c("D-GCs", "M-GCs", "S-GCs")) +
  scale_y_continuous(breaks = seq(-4, 2, 1), limits = c(-4.8, 2)) +
  scale_color_manual(values = colors_S_GC_high_low, drop = FALSE, labels = c("3-1", "3-2", "4-1", "4-2")) +
  scale_shape_manual(values = c(21, 22, 23, 24), drop = FALSE, labels = c("3-1", "3-2", "4-1", "4-2")) +
  labs(x = "GC Subtype", y = "log2(FC normalized by D-GCs)", color = "Type", shape = "Type") +
  trends_theme

Continuous_D_GC_high_low_trends_graph <- ggplot(Continuous_D_GC_high_low_trends,
                                                aes(x = Time, y = Log2FC, group = interaction(Type, Gene_Symbol), color = Type)) +
  geom_line(alpha = 0.6) +
  geom_point(aes(shape = Type), stroke = 1, size = 2, fill = "white") +
  geom_line(data = Continuous_D_GC_high_low_avg, aes(x = Time, y = Average_Log2FC, group = Type),
            color = "black", linewidth = 1, linetype = "solid", lineend = "round", inherit.aes = FALSE) +
  geom_line(data = Continuous_D_GC_high_low_avg,
            aes(x = Time, y = Average_Log2FC, group = Type, color = Type),
            linewidth = 0.5, linetype = "dashed", lineend = "round", inherit.aes = FALSE) +
  geom_hline(yintercept = c(0, log2(1.5), -log2(1.5)), linewidth = 0.5, colour = "darkgray", linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 1, 2), labels = c("D-GCs", "M-GCs", "S-GCs")) +
  scale_y_continuous(breaks = seq(-4, 2, 1), limits = c(-3.8, 2.4)) +
  scale_color_manual(values = colors_D_GC_high_low, drop = FALSE, labels = c("5-1", "5-2", "6-1", "6-2")) +
  scale_shape_manual(values = c(21, 22, 23, 24), drop = FALSE, labels = c("5-1", "5-2", "6-1", "6-2")) +
  labs(x = "GC Subtype", y = "log2(FC normalized by D-GCs)", color = "Type", shape = "Type") +
  trends_theme

Continuous_M_GC_high_low_trends_graph <- ggplot(Continuous_M_GC_high_low_trends,
                                                aes(x = Time, y = Log2FC, group = interaction(Type, Gene_Symbol), color = Type)) +
  geom_line(alpha = 0.6) +
  geom_point(aes(shape = Type), stroke = 1, size = 2, fill = "white") +
  geom_line(data = Continuous_M_GC_high_low_avg, aes(x = Time, y = Average_Log2FC, group = Type),
            color = "black", linewidth = 1, linetype = "solid", lineend = "round", inherit.aes = FALSE) +
  geom_line(data = Continuous_M_GC_high_low_avg, aes(x = Time, y = Average_Log2FC, group = Type, color = Type),
            linewidth = 0.5, linetype = "dashed", lineend = "round", inherit.aes = FALSE) +
  geom_hline(yintercept = c(0, log2(1.5), -log2(1.5)), linewidth = 0.5, colour = "darkgray", linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 1, 2), labels = c("D-GCs", "M-GCs", "S-GCs")) +
  scale_y_continuous(breaks = seq(-1, 1, 1), limits = c(-1.1, 1)) +
  scale_color_manual(values = colors_M_GC_high_low, drop = FALSE, labels = "8") +
  scale_shape_manual(values = c(21, 22, 23, 24), drop = FALSE, labels = "8") +
  labs(x = "GC Subtype", y = "log2(FC normalized by D-GCs)", color = "Type", shape = "Type") +
  trends_theme

# Figure 4b, Pie chart and Sunburst chart for trend analysis
counts_all <- Total %>% group_by(Class, Type) %>% summarise(Count = n_distinct(Gene_Symbol), .groups = "drop")

colors_group <- c("Continuous decrease" = "#C63A62", "High in D-GCs" = "#D03E32", "Low in S-GCs" = "#B4556B",
                  "Continuous increase" = "#169E6C", "High in S-GCs" = "#33B8A0", "Low in D-GCs" = "#5A9543",
                  "High in M-GCs" = "#31859C", "Low in M-GCs" = "#1A237E")
group_order <- c("Continuous decrease", "High in D-GCs", "Low in S-GCs", "Continuous increase",
                 "High in S-GCs", "Low in D-GCs", "High in M-GCs", "Low in M-GCs")
type_order <- c("2-1; D > M > S", "2-2; D ≥ M ≥ S, D > S", "6-2; D > M ≤ S", "6-1; D > M ≥ S", "4-1; D ≥ M > S", "4-2; D ≤ M > S",
                "1-1; D < M < S", "1-2; D ≤ M ≤ S, D < S", "3-2; D ≥ M < S", "3-1; D ≤ M < S", "5-1; D < M ≤ S", "5-2; D < M ≥ S",
                "7; D < M > S", "8; D > M < S")
sunburst_source <- counts_all %>% rename(Group = Class) %>% select(Type, Count, Group) %>%
                          mutate(Group = factor(Group, levels = group_order), Type = factor(Type, levels = type_order)) %>%
                          arrange(Group, Type)
                                                           
sunburst_outer_circle <- sunburst_source %>% mutate(total = sum(Count), fraction = Count / total, ymax     = cumsum(fraction), ymin = lag(ymax, default = 0), mid_y = (ymin + ymax) / 2)

group_num_map <- c("Continuous decrease" = "2", "High in D-GCs" = "6", "Low in S-GCs" = "4", "Continuous increase" = "1",
                   "High in S-GCs" = "3", "Low in D-GCs" = "5", "High in M-GCs" = "7", "Low in M-GCs" = "8")

sunburst_inner_circle <- sunburst_outer_circle %>% group_by(Group) %>% summarise(ymin  = min(ymin), ymax  = max(ymax),
                                                                    Count = sum(Count),  .groups = "drop") %>%
  mutate(mid_y = (ymin + ymax) / 2, percentage   = Count / sum(Count) * 100,
         inner_label  = paste0(group_num_map[as.character(Group)], "\n(", round(percentage, 1), "%)"))
sunburst_outer_circle <- sunburst_outer_circle %>% mutate(outer_label = sub(";.*", "", as.character(Type))) %>% ungroup()

sunburst_inner_circle_plot <- sunburst_inner_circle %>% mutate(fill_var = as.character(Group))
sunburst_outer_circle_plot <- sunburst_outer_circle %>% mutate(fill_var = as.character(Type))
all_colors <- c(colors_group, colors_combined_Trends)

sunburst_plot <- ggplot() + geom_rect(data = sunburst_inner_circle_plot, aes(xmin = 1.0, xmax = 2.3, ymin = ymin, ymax = ymax,
                                                                             fill = fill_var), color = "white", linewidth = 0.6) +
  geom_rect(data = sunburst_outer_circle_plot, aes(xmin = 2.3, xmax = 3.5, ymin = ymin, ymax = ymax, fill = fill_var),
            color = "white", linewidth = 0.3) +
  geom_shadowtext(data = sunburst_inner_circle_plot, aes(x = 1.65, y = mid_y, label = inner_label), size = 3, 
                  lineheight = 0.9, fontface = "bold", color = "black", bg.color = "white", bg.r = 0.055) +
  geom_shadowtext(data = sunburst_outer_circle_plot %>% filter(!outer_label %in% c("7", "8")),
                  aes(x = 2.9, y = mid_y, label = outer_label), size = 2.8, fontface = "bold", color = "black",
                  bg.color = "white", bg.r = 0.055) +
  geom_text_repel(data = sunburst_outer_circle_plot %>% filter(outer_label %in% c("7", "8")),
                  aes(x = 3.5, y = mid_y, label = outer_label), nudge_x = 0.7, direction = "y", size = 2.8,
                  fontface = "bold", color = "black", bg.color = "white", bg.r = 0.055, segment.color = "gray50", segment.size  = 0.3) +
  scale_fill_manual(values = all_colors) +
  coord_polar(theta = "y") +
  xlim(c(0, 5.5)) +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

print(sunburst_plot)

pie_map <- c("2" = "D-GC High", "4" = "D-GC High", "6" = "D-GC High", "1" = "S-GC High",
             "3" = "S-GC High", "5" = "S-GC High", "8" = "M-GC Low")

colors_pie <- c("D-GC High" = "#D03E32", "S-GC High" = "#169E6C",  "M-GC Low" = "#1A237E")

pie_data <- sunburst_inner_circle %>% mutate(num = group_num_map[as.character(Group)], Category = pie_map[num]) %>%
  filter(!is.na(Category)) %>% group_by(Category) %>% summarise(Count = sum(Count), .groups = "drop") %>%
  mutate(Category = factor(Category, levels = c("D-GC High", "S-GC High", "M-GC Low"))) %>%
  arrange(Category) %>% mutate(pct = Count / sum(Count) * 100, end = 2 * pi * cumsum(Count) / sum(Count),
                               start = lag(end, default = 0), mid = (start + end) / 2,
                               type_label = case_when(Category == "D-GC High" ~ "(types 2, 4, 6)",
                                                      Category == "S-GC High" ~ "(types 1, 3, 5)",
                                                      Category == "M-GC Low"  ~ "(type 8)"),
                               label = paste0(Category, "\n(", sprintf("%.1f", pct), "%)\n", type_label))
print(pie_data)

r_pie   <- 1.0
r_label <- 0.6
pie_inner_lab <- pie_data %>% filter(Category != "M-GC Low") %>% mutate(lx = r_label * sin(mid),ly = r_label * cos(mid))
m_row   <- pie_data %>% filter(Category == "M-GC Low")
arrow_x <- 0.92 * sin(m_row$mid)
arrow_y <- 0.92 * cos(m_row$mid)
top_y   <- arrow_y + (1.35 - arrow_y) / 2
m_label <- paste0("M-GC Low\n(", sprintf("%.1f", m_row$pct), "%)\n(type 8)")

pie_plot <- ggplot() + geom_arc_bar(data = pie_data, aes(x0 = 0, y0 = 0, r0 = 0, r = r_pie, start = start,
                                                         end = end, fill = Category), color = "white", linewidth = 0.1) +
  geom_text(data = pie_inner_lab, aes(x = lx, y = ly, label = label), size = 3.15,
            lineheight = 0.9, fontface = "plain", color = "black") +
  annotate("segment", x = arrow_x, xend = arrow_x, y = top_y, yend = arrow_y, color = "black", linewidth = 0.5) +
  annotate("shadowtext", x = arrow_x, y = top_y + 0.03, label = m_label, vjust = 0, size = 3.3, fontface = "plain",
           color = "black", bg.color = "white", bg.r = 0.055) +
  scale_fill_manual(values = colors_pie) +
  coord_fixed(xlim = c(-1.3, 1.3), ylim = c(-1.1, 1.6), clip = "off") +
  theme_void() +
  theme(legend.position = "none")


print(pie_plot)

# Export plots to powerpoint and excel
Figure_4 <- read_pptx() %>%
  add_slide(layout = "Blank", master = "Office Theme") %>%
  ph_with(value = dml(ggobj = Continuous_inc_dec_trends_graph),
          location = ph_location(width = 2.6, height = 2.7, left = 1, top = 0)) %>%
  ph_with(value = dml(ggobj = Continuous_S_GC_high_low_trends_graph),
          location = ph_location(width = 2.6, height = 2.7, left = 3.6, top = 0)) %>%
  ph_with(value = dml(ggobj = Continuous_D_GC_high_low_trends_graph),
          location = ph_location(width = 2.6, height = 2.7, left = 1, top = 3)) %>%
  ph_with(value = dml(ggobj = Continuous_M_GC_high_low_trends_graph),
          location = ph_location(width = 2.6, height = 2.7, left = 3.6, top = 3)) %>%
  add_slide(layout = "Blank", master = "Office Theme") %>%
  ph_with(value = dml(ggobj = pie_plot),      location = ph_location(width = 2.5, height = 2.5, left = 2, top = 0)) %>% 
  ph_with(value = dml(ggobj = sunburst_plot), location = ph_location(width = 4, height = 4, left = 5, top = 0))

print(Figure_4, target = "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Figure 4, GC Trends Figures.pptx")
write_xlsx(Total, "C:/Users/mw/Desktop/mw/Figure/Summarized figures/Figure, Statistics/Supplement Data 4, List of DEGs, pattern analysis.xlsx")