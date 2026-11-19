# Figure 3. Cross-tissue colocalization of GTEx muscle H3K27ac haQTLs
# with mQTLs in nine GTEx tissues (panels A-E).
# Run from the directory containing the input files. TIFFs are written there.
# Input tables are precomputed; their upstream generation is not included.
# Strong colocalization: PP4 >= 0.8.

library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)

# -----------------------------------------------------------------
# Panel A
# Input: panelA_adjust_number_for_plot.txt (Tissue, sample_size, count_coloced_haqtl).
# -----------------------------------------------------------------
data=read.table("panelA_adjust_number_for_plot.txt",head=T,as.is=T)

custom_colors <- c("#99ff00", "#cc9955", "#ffaaff", "#dddddd", "#33cccc", 
                   "#22ffdd", "#aaaaff", "#ff00bb", "#aaaaaa")

data$Tissue <- factor(data$Tissue,
                      levels = c("9_Lung","8_Colon","7_Ovary","6_Prostate",
                                 "5_Breast","2_Kidney","1_Muscle","4_Blood","3_Testis"))

tiff("Panel_A_scatter_plot_dev3.tiff", width = 7, height = 6, units = "in", res = 300)

p <- ggplot(data, aes(x = sample_size, y = count_coloced_haqtl, color = Tissue)) +
  geom_point(size = 6, alpha = 0.6) +
  labs(
    title = "",
    x = "mQTL sample size in each tissue",
    y = "muscle haQTLs colocalized with mQTLs",
    color = "Tissues (# of gCpGs)"
  ) +
  scale_color_manual(
    values = custom_colors,
    breaks =  c("9_Lung","8_Colon","7_Ovary","6_Prostate",
                "5_Breast","2_Kidney","1_Muscle","4_Blood","3_Testis"),
    labels = c("lung (206,802)","colon (202,258)","ovary (185,408)","prostate (161,320)",
               "breast (141,631)","kidney (138,748)","muscle (123,705)","whole blood (116,509)","testis (108,844)")
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.75,0.3),
    legend.title = element_text(size = 20, color = "black", family = "Arial"),
    legend.text  = element_text(size = 20, color = "black", family = "Arial"),
    axis.title.x = element_text(size = 22, color = "black", family = "Arial", margin = margin(t = 15)),
    axis.title.y = element_text(size = 22, color = "black", family = "Arial", margin = margin(r = 15)),
    axis.text.y  = element_text(size = 18, color = "black", family = "Arial"),
    axis.text.x  = element_text(size = 18, color = "black", family = "Arial")
  )

print(p)
dev.off()

# -----------------------------------------------------------------
# Panel B
# Published counts for 485 muscle gHMs.
# -----------------------------------------------------------------

data <- data.frame(
  group = LETTERS[1:5],
  value = c(80, 37, 173, 145, 50)
)

purple_pal <- colorRampPalette(c("#ffffff", "#aaaaff", "#000066"))(5)

tiff("Panel_B_pie_chart2.tiff",
     width = 14, height = 8, units = "in", res = 300)

p <- ggplot(data, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  geom_text(aes(label = sprintf("%d\n(%d%%)", value, round(100 * value / sum(value))),
                color = group), position = position_stack(vjust = 0.5),
            family = "Arial", size = 5, show.legend = FALSE) +
  scale_color_manual(values = c("black", "black", "black", "white", "white")) +
  coord_polar("y", start = 0) +
  scale_fill_manual(
    values = purple_pal,
    labels = c(
      "no colocalization test",
      "no strong colocalization",
      "colocalized (PP4\u22650.8) with mQTLs in 1-3 tissues",
      "colocalized (PP4\u22650.8) with mQTLs in 4-6 tissues",
      "colocalized (PP4\u22650.8) with mQTLs in 7-9 tissues"
    )
  ) +
  labs(fill = "Colocalization status for muscle haQTL") +
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 16, color = "black", family = "Arial"),
    legend.text  = element_text(size = 14, color = "black", family = "Arial")
  ) +
  guides(
    fill = guide_legend(
      nrow = 5,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5
    )
  )

print(p)
dev.off()

# -----------------------------------------------------------------
# Panel C
# Published percentages; upstream numerator/denominator calculations not included.
# -----------------------------------------------------------------
data <- data.frame(
  name = factor(c("muscle", "kidney", "testis", "whole blood", "ovary", "lung", "colon", "prostate", "breast"),
                levels = c("muscle", "kidney", "testis", "whole blood", "ovary", "lung", "colon", "prostate", "breast")),
  value = c(2.27, 0.34, 0.48, 0.53, 0.69, 0.72, 0.73, 0.81, 0.85)
)

tiff("Panel_C_10neg6.tiff", width = 10, height = 6, units = "in", res = 300)

p <- ggplot(data, aes(x = name, y = value, fill = name)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = value), vjust = -0.5, size = 7) +

  scale_y_continuous(
    limits = c(0, 2.5),
    expand = c(0, 0)
  ) +
  
  scale_fill_manual(values = c("#aaaaff", "#22ffdd", "#aaaaaa", "#ff00bb", "#ffaaff", 
                               "#99ff00", "#cc9955", "#dddddd", "#33cccc")) +
  
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.x = element_text(size = 18, color = "black", family = "Arial",
                                margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, color = "black", family = "Arial",
                                margin = margin(r = 10)),
    axis.text.y = element_text(size = 18, color = "black", family = "Arial"),
    axis.text.x = element_text(size = 18, color = "black", family = "Arial",
                               angle = 45, hjust = 1),
    legend.position = "none",

    plot.margin = margin(t = 40, r = 10, b = 10, l = 10)
  ) +
  
  labs(
    x = "",
    y = "% mQTLs colocalized with muscle haQTLs"
  )

print(p)
dev.off()

# -----------------------------------------------------------------
# Panel D
# Input: n_coloc_tissues_by_type_gAREs.txt (type, n_coloced_tissues). Percentages within type.
# -----------------------------------------------------------------
a <- read.table("n_coloc_tissues_by_type_gAREs.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

summary_data <- a %>%
  count(type, n_coloced_tissues, name = "count") %>%
  group_by(type) %>%
  mutate(percentage = 100 * count / sum(count)) %>%
  ungroup()

fill_colors <- c("haQTL-shared_gAREs" = "white",
                 "muscle-specific_gAREs" = "#aaaaff")

summary_data$type <- factor(summary_data$type,
                            levels = c("muscle-specific_gAREs", "haQTL-shared_gAREs"))

tiff("Panel_newD_percentage.tiff", width = 10, height = 7, units = "in", res = 300)

p <- ggplot(summary_data, aes(x = factor(n_coloced_tissues, levels = 0:9),
                         y = percentage, fill = type)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.9),
           color = "black", width = 0.7) +
  scale_fill_manual(
    values = fill_colors,
    breaks = c("muscle-specific_gAREs", "haQTL-shared_gAREs"),
    labels = c("muscle-specific gHMs", "tissue-shared gHMs")
  ) +
  labs(x = "# of tissues in which mQTLs colocalized with muscle haQTLs",
       y = "percentage of gHMs",
       fill = " ") +
  scale_y_continuous(limits = c(0, 30), expand = expansion(mult = c(0, 0.05))) +
  theme_classic(base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text  = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 24, color = "black"),
    legend.title = element_text(size = 24, color = "black"),
    legend.text  = element_text(size = 20, color = "black"),

    legend.position = c(0.65, 0.85),
    legend.background = element_rect(fill = "white", color = NA)
  )

print(p)
dev.off()

# -----------------------------------------------------------------
# Panel E
# Input: panel_D_data_new.txt (status, cate). Pearson chi-square test.
# -----------------------------------------------------------------
data=read.table("panel_D_data_new.txt", head=T,as.is=T)

tab <- table(data$status, data$cate)
print(tab)
print(chisq.test(tab))

base_family <- "Arial"

data$cate   <- factor(data$cate, levels = c("no_test", "no_coloc", "1-3", "4-6", "7-9"))
data$status <- factor(data$status, levels = c("enhancer", "promoter"))

plot_data <- data %>%
  count(status, cate, name = "n") %>%
  tidyr::complete(status, cate, fill = list(n = 0))

spacer <- tibble(
  status = factor("spacer", levels = c(levels(data$status), "spacer")),
  cate   = factor("no_test", levels = levels(data$cate)),
  n = NA_real_
)
plot_data_spaced <- bind_rows(plot_data, spacer)

plot_data_spaced$group <- interaction(plot_data_spaced$status, plot_data_spaced$cate, sep = "_", drop = TRUE)

x_levels <- c(
  paste("enhancer", levels(data$cate), sep = "_"),
  "spacer_no_test",
  paste("promoter", levels(data$cate), sep = "_")
)

plot_data_spaced$group <- factor(plot_data_spaced$group, levels = x_levels)

purple_pal  <- colorRampPalette(c("#ffffff", "#aaaaff", "#000066"))(5)
legend_order <- c("no_test", "no_coloc", "1-3", "4-6", "7-9")
legend_labels <- c(
  "no colocalization test",
  "no strong colocalization",
  "colocalized with mQTLs in 1-3 tissues",
  "colocalized with mQTLs in 4-6 tissues",
  "colocalized with mQTLs in 7-9 tissues"
)

tiff("Panel_E_dev2.tiff", width = 8, height = 3.5, units = "in", res = 300)

p <- ggplot(plot_data_spaced, aes(x = group, y = n, fill = cate)) +
  geom_bar(stat = "identity", color = "black", na.rm = TRUE) +
  scale_x_discrete(limits = x_levels, labels = c("", "", "enhancer", "", "", "", "", "", "promoter", "", ""), drop = FALSE) +
  scale_fill_manual(
    values = purple_pal,
    breaks = legend_order,
    labels = legend_labels,
    na.translate = FALSE
  ) +
  labs(
    x = NULL,
    y = "# of muscle gHMs",
    fill = "colocalization status for muscle haQTL"
  ) +
  theme_classic(base_family = base_family) +
  coord_cartesian(clip = "off") +
  theme(

    legend.position = c(0.6, 0.95),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.title = element_text(size = 16, color = "black"),
    legend.text  = element_text(size = 14, color = "black"),
    legend.key.size = unit(12, "pt"),
    legend.spacing.y = unit(2, "pt"),
    
    plot.margin = margin(t = 10, r = 40, b = 10, l = 10),
    
    axis.title.y = element_text(size = 20, color = "black"),
    axis.text.y  = element_text(size = 20, color = "black"),
    axis.text.x  = element_text(size = 16, color = "black"),
    plot.title   = element_text(size = 18, color = "black", hjust = 0.5)
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      byrow = TRUE,
      title.position = "top"
    )
  )

print(p)
dev.off()
