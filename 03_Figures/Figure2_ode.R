# Figure 2. Colocalization of GTEx skeletal muscle H3K27ac haQTLs
# with skeletal muscle mQTLs (primary haQTL threshold P < 1e-6).
# Run from the directory containing the prepared input tables below.
# Requires ggplot2, dplyr and Arial installed on the plotting system.
# Panels are exported separately for assembly into the manuscript figure.
# Prepared tables must correspond to the primary coloc.abf analysis.
library(ggplot2)
library(dplyr)
library(grid)

out_dir <- 'Figure2_output'
dir.create(out_dir, showWarnings = FALSE)
read_input <- function(path) read.table(path, header = TRUE, stringsAsFactors = FALSE)
save_panel <- function(p, filename, width, height, units = 'in', res = 300) {
  tiff(file.path(out_dir, filename), width = width, height = height,
       units = units, res = res)
  on.exit(dev.off())
  print(p)
}
status_levels <- c('1_strong_coloc','2_suggestive_coloc','3_no_coloc','4_no_coloc_test')
status_colors <- setNames(c('#aaaaff','#33dd33','#ffaa00','#aaaaaa'), status_levels)
status_labels <- c('strong colocalization (PP4\u22650.8)',
                   'suggestive colocalization (0.5\u2264PP4<0.8)',
                   'no colocalization (PP4<0.5)', 'no colocalization test')
axis_theme <- function(title_size, text_size, line_width = 0.5, ratio = 1) {
  theme(panel.background = element_rect(fill = 'white', color = 'white'),
        axis.line = element_line(color = 'black', linewidth = line_width),
        axis.title = element_text(size = title_size, color = 'black', family = 'Arial'),
        axis.text = element_text(size = text_size, color = 'black', family = 'Arial'),
        aspect.ratio = ratio, legend.background = element_rect(fill = 'white', color = 'white'))
}

# ============================================================
# Figure 2A: Colocalization status of 485 haQTLs
# Counts transcribed from the primary manuscript figure, not SuSiE.
# ============================================================
a <- data.frame(group = factor(status_levels, levels = status_levels),
                value = c(169,54,37,225))
a$label <- sprintf('%d\n(%d%%)', a$value, round(100*a$value/sum(a$value)))
p_a <- ggplot(a, aes(x = '', y = value, fill = group)) +
  geom_col(width = 1, color = 'white') +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), family = 'Arial') +
  coord_polar('y', start = 0) +
  scale_fill_manual(values = status_colors, labels = status_labels) +
  labs(fill = 'haQTL colocalization status') + theme_void(base_family = 'Arial') +
  theme(legend.position = 'top', legend.box = 'vertical', legend.direction = 'vertical',
        legend.title.position = 'top', legend.title = element_text(size = 18, color = 'black'),
        legend.text = element_text(size = 16, color = 'black'), legend.text.align = 0,
        legend.box.margin = margin(l = 20), legend.box.spacing = unit(-4,'pt'),
        plot.margin = margin(10,40,10,40)) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE, title.position = 'top', title.hjust = 0.05))
save_panel(p_a, 'Figure2A.tiff', 12, 5)

# ============================================================
# Figure 2B: Number of CpGs impacted per strongly colocalized haQTL
# Counts supplied in the original plotting script (sum = 169).
# ============================================================
b <- data.frame(name = factor(c('1','2','3','4','5','6','7-26'),
                             levels = c('1','2','3','4','5','6','7-26')),
                value = c(76,31,15,11,11,10,15))
p_b <- ggplot(b, aes(name,value)) + geom_col(fill = '#aaaaff') +
  geom_text(aes(label = value), vjust = -0.5, size = 6.5, family = 'Arial') +
  axis_theme(22,20) + labs(x = '# of CpGs impacted per haQTL', y = 'count of haQTL')
save_panel(p_b, 'Figure2B.tiff', 10, 10)

# ============================================================
# Figure 2C: Colocalization status by enhancer/promoter annotation
# Input: Panel_newC_data.txt, one row per haQTL; status, coloc_status.
# Original plotting code; 2x2 Pearson test evaluates test availability.
# ============================================================
c <- read_input('Panel_newC_data.txt')
c$status <- factor(c$status, levels = c('enhancer','promoter'))
c$coloc_status <- factor(c$coloc_status, levels = status_levels)
stopifnot(!anyNA(c$status), !anyNA(c$coloc_status))
c_summary <- c %>% count(status, coloc_status, name = 'count')
print(table(c$status,c$coloc_status))
print(chisq.test(table(c$status,c$coloc_status)))
c$tested <- ifelse(c$coloc_status == '4_no_coloc_test','No','Yes')
print(chisq.test(table(c$status,c$tested)))
print(fisher.test(table(c$status,c$tested)))
p_c <- ggplot(c_summary, aes(status, count, fill = coloc_status)) +
  geom_col(position = position_dodge2(width = 0.7, preserve = "single"), color = "black") +
  theme_classic(base_family = "Arial") +
  labs(x = NULL, y = "count of haQTL", fill = "haQTL colocalization status") +
  scale_fill_manual(values = status_colors, labels = status_labels) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  theme(plot.margin = margin(1, 5.5, 1, 1, "cm"),
        legend.position = "top", legend.justification = c(0.2, 1),
        legend.box = "vertical", legend.direction = "horizontal",
        legend.title.position = "top", legend.title = element_text(size = 16),
        legend.text = element_text(size = 14), axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20, color = "black"),
        axis.line = element_line(color = "black"), panel.grid = element_blank())
save_panel(p_c, "Figure2C.tiff", 10, 10)

# ============================================================
# Figure 2D: Distance and effect direction for colocalized pairs
# Input: Panel_newD_data.txt. Assumes direction coding 0=opposite,
# 1=same (verify against upstream preparation).
# ============================================================
d <- read_input('Panel_newD_data.txt')
d <- subset(d, abs(gare_cpg_dist) < 50000 & palindromic == 0 & !is.na(same_direct_haqtl_mqtl))
stopifnot(all(as.character(d$same_direct_haqtl_mqtl) %in% c('0','1')))
d$direction <- factor(as.character(d$same_direct_haqtl_mqtl), levels = c('0','1'))
d$gare_cpg_dist_kb <- d$gare_cpg_dist/1000
n_direction <- table(d$direction)
print(n_direction)
# Manuscript legend counts; warn if the prepared data do not reproduce them.
if (!identical(as.integer(n_direction), c(310L, 71L)))
  warning("Panel D counts differ from manuscript (310 opposite, 71 same).")
p_d <- ggplot(d,aes(gare_cpg_dist_kb,fill = direction)) +
  geom_histogram(position = 'identity',alpha = 0.5,color = 'grey',bins = 30) +
  scale_fill_manual(values = c('0'='#aaaaff','1'='#33dd33'), breaks = c('0','1'),
                    labels = c('opposite direction (310 pairs)', 'same direction (71 pairs)')) +
  guides(fill = guide_legend(title = 'Effect directions of haQTL and mQTL',nrow = 2,byrow = TRUE)) +
  labs(x = 'distance between gHM and gCpG (kb)',y = 'count of colocalizations') +
  axis_theme(29,28,1) + theme(plot.margin = margin(30,30,15,15),
    axis.title.y = element_text(size = 30), legend.position = 'top',legend.direction = 'vertical',
    legend.box = 'vertical',legend.box.margin = margin(3,10,3,10),
    legend.title = element_text(size = 30,family = 'Arial',hjust = 0.5),
    legend.text = element_text(size = 30,family = 'Arial'))
save_panel(p_d,'Figure2D.tiff',10,9)

# ============================================================
# Figure 2E: Pearson correlation of effect estimates
# Input: Panel_newF_data.txt, prepared opposite-direction pairs
# with alleles already aligned upstream. Exclude palindromic SNPs.
# ============================================================
e <- read_input('Panel_newF_data.txt')
e <- subset(e,palindromic == 0 & is.finite(haqtl_beta) & is.finite(mqtl_slope))
# Verify input has the intended allele orientation before interpreting.
if (any(e$haqtl_beta * e$mqtl_slope >= 0)) warning('Panel E input contains non-opposite effects; verify upstream selection.')
print(nrow(e))
print(cor.test(e$haqtl_beta,e$mqtl_slope,method = 'pearson'))
p_e <- ggplot(e,aes(haqtl_beta,mqtl_slope)) + geom_point() + geom_smooth(method = 'lm') +
  labs(x = 'beta coefficient of haQTL',y = 'beta coefficient of mQTL') +
  axis_theme(24,22,1.2,1.2) + theme(axis.line = element_line(color = 'grey30',linewidth = 1.2))
save_panel(p_e,'Figure2E.tiff',10,10)

# ============================================================
# Figure 2F: Methylation levels and allelic effects near gHMs
# Exact final plotting version, using prepared new_data_for_plot.txt.
# The original preprocessing-to-long-format step was not supplied.
# Required columns: cpg_gare_dist, new_mval_mean, gARE_match.
# ============================================================
new <- read_input("new_data_for_plot.txt")
new$new_mval_mean <- as.numeric(new$new_mval_mean)
new <- new[!is.na(new$new_mval_mean), ]
new$new_beta <- 2^new$new_mval_mean / (1 + 2^new$new_mval_mean)
new <- new[!is.na(new$new_beta), ]

# Assign dot sizes
new$dot_size <- ifelse(new$gARE_match == "other_cpg", 0.8, 1.8)

# Calculate mean_beta for non-grey groups
mean_lines <- new %>%
  filter(gARE_match != "other_cpg") %>%
  group_by(gARE_match) %>%
  summarise(mean_beta = mean(new_beta, na.rm = TRUE)) %>%
  mutate(x = -1000, xend = 1000)  # Add x positions for segment lines

# Plot

new$gARE_match <- factor(new$gARE_match,
                         levels = c("other_cpg", "coloced_cpg", "upper", "lower"))

p_f <- ggplot(new, aes(x = cpg_gare_dist, y = new_beta, color = gARE_match)) +
  geom_point(aes(size = dot_size), alpha = 0.6) +
  scale_size_identity() +
  geom_smooth(data = subset(new, gARE_match == "other_cpg"),
              method = "loess", se = FALSE) +
  geom_segment(data = mean_lines, inherit.aes = FALSE,
               aes(x = x, xend = xend, y = mean_beta, yend = mean_beta, color = gARE_match),
               linewidth = 1, linetype = "solid") +
  theme_classic(base_family = "Arial") +
  labs(
    title = " ",
    x = "CpG-gHM distance",
    y = "CpG beta mean",
    color = " "
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.x = element_text(size = 26, color = "black",family = "Arial"),
    axis.title.y = element_text(size = 26, color = "black",family = "Arial"),
    axis.text.x  = element_text(size = 24, color = "black",family = "Arial"),
    axis.text.y  = element_text(size = 24, color = "black",family = "Arial"),
    legend.position = "top",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 16,face="bold"),
legend.margin = margin(t = -15, b = 0),          # pulls legend closer vertically
  legend.box.margin = margin(0, 0, -15, 0)         # reduces space between plot and legend
  ) +
  scale_color_manual(
    values = c(
      "coloced_cpg" = "red",
      "upper"       = "darkorchid1",
      "lower"       = "deepskyblue",
      "other_cpg"   = "grey60"
    ),
    breaks = c("coloced_cpg", "upper", "lower", "other_cpg"),
    labels = c("Colocalized CpGs", "Upper", "Lower", "Other CpGs")
  )


save_panel(p_f, "Figure2F.tiff", 800, 600, units = "px", res = 72)

# ============================================================
# Figure 2G: Measured CpGs within 10 kb of each gHM
# Input: Panel_newE_data.txt. The supplied code was for P < 1e-7;
# use the full primary P < 1e-6 table here, ONLY after confirming it
# contains all 485 primary haQTLs and that coloc coding is as expected.
# ============================================================
g <- read_input('Panel_newE_data.txt')
stopifnot(nrow(g) == 485L, !anyDuplicated(g$gare))
print(table(g$coloc, useNA = 'ifany'))
# Preserve original grouping: two or more, one, uncertain, no test.
# Verify that the factor level order matches these four groups.
g$coloc <- factor(g$coloc)
stopifnot(nlevels(g$coloc) == 4L, !anyNA(g$coloc))
print(levels(g$coloc))
# Compare each of the first three groups with the no-test group.
for (category in levels(g$coloc)[1:3]) {
  print(t.test(g$n_cpg_10kb_eachside[g$coloc == category],
               g$n_cpg_10kb_eachside[g$coloc == levels(g$coloc)[4]], var.equal = TRUE))
}
p_g <- ggplot(g,aes(coloc,n_cpg_10kb_eachside,fill = coloc)) +
  geom_violin(trim = FALSE,color = 'black',alpha = 0.7) +
  geom_boxplot(width = 0.1,outlier.shape = NA,color = 'black',fatten = 0) +
  theme_classic(base_family = 'Arial') +
  theme(plot.margin = margin(2,2,2,2,'cm'),legend.position = 'none',
        axis.title.x = element_text(size = 18,margin = margin(t = 20)),
        axis.title.y = element_text(size = 18,margin = margin(r = 10)),
        axis.text.y = element_text(size = 16,color = 'black'),
        axis.text.x = element_blank(),
        axis.line = element_line(color = 'black')) +
  labs(x = NULL,y = 'nearby measured CpGs per gHM (within 10 kb)') +
  scale_fill_manual(values = unname(status_colors),
    labels = c("haQTLs colocalized with 2+ mQTLs", "haQTLs colocalized with 1 mQTL",
               "haQTLs with uncertain colocalizations", "haQTLs with no colocalization test")) +
  stat_summary(fun = mean,geom = 'point',size = 5,color = 'purple',shape = 18)
save_panel(p_g,'Figure2G.tiff',10,10)
