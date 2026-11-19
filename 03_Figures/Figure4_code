# Figure 4. ROSMAP DLPFC H3K9ac haQTL and 450K mQTL colocalization
# Run from the directory containing the panel input files.
# Requires ggplot2, dplyr, and Arial installed on the system.
library(ggplot2)
library(dplyr)

status_levels <- c('1_strong_coloc', '2_suggestive_coloc', '3_no_coloc', '4_no_coloc_test')
status_labels <- c('strong colocalization (PP4\u22650.8)',
                   'suggestive colocalization (0.5\u2264PP4<0.8)',
                   'no colocalization (PP4<0.5)', 'no colocalization test')
palette <- c('#aaaaff', '#33dd33', '#ffaa00', '#aaaaaa')
read_panel <- function(path) read.table(path, header = TRUE, stringsAsFactors = FALSE)
save_panel <- function(p, path, width = NULL, height = NULL) {
  if (is.null(width)) tiff(path) else
    tiff(path, width = width, height = height, units = 'in', res = 300)
  on.exit(dev.off())
  print(p)
}

# ---- Figure 4A: Colocalization status ----
# Select the highest PP4 per tested peak; report category counts.
a <- read_panel('PanelA_data.txt') %>%
  group_by(Peak) %>% slice_max(pp4, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(category = case_when(pp4 >= 0.8 ~ status_levels[1],
                              pp4 >= 0.5 ~ status_levels[2],
                              pp4 < 0.5 ~ status_levels[3]))
print(count(a, category))
# Published figure counts, including untested peaks.
a_plot <- data.frame(category = factor(status_levels, levels = status_levels),
                     count = c(1396, 56, 285, 632))
a_plot$label <- sprintf('%d\n(%d%%)', a_plot$count, round(100*a_plot$count/sum(a_plot$count)))
p <- ggplot(a_plot, aes(x = '', y = count, fill = category)) +
  geom_col(width = 1, color = 'white') +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), family = 'Arial', size = 5) +
  coord_polar(theta = 'y') +
  scale_fill_manual(values = palette, labels = status_labels) +
  labs(fill = 'haQTL colocalization status') + theme_void(base_family = 'Arial') +
  theme(legend.position = 'top', legend.direction = 'vertical',
        legend.title = element_text(size = 18), legend.text = element_text(size = 16)) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))
save_panel(p, 'Figure4A.tiff', 8, 7)

# ---- Figure 4B: Number of CpGs impacted per haQTL ----
# Input is the prepared set of strongly colocalized peak-CpG pairs.
b <- read_panel('PanelB_data.txt') %>% count(Peak, name = 'n_cpg')
print(count(b, n_cpg, name = 'n_peaks'))
# Published panel B counts sum to 1391, not panel A's 1396.
# Preserve the actual published bars; the five-peak difference needs source-data reconciliation.
b_plot <- data.frame(name = factor(c('1','2','3','4','5','6','7-36'),
                                   levels = c('1','2','3','4','5','6','7-36')),
                     value = c(437,305,171,143,81,82,172))
if (sum(b_plot$value) != a_plot$count[1]) warning('Figure 4B published counts differ from Figure 4A by five peaks.')
p <- ggplot(b_plot, aes(name, value)) + geom_col(fill = '#EEEE00') +
  geom_text(aes(label = value), vjust = -0.5, size = 6, family = 'Arial') +
  theme_classic(base_family = 'Arial') +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18, color = 'black'), aspect.ratio = 1) +
  labs(x = '# of CpGs impacted per haQTL', y = 'count of haQTL')
save_panel(p, 'Figure4B.tiff')

# ---- Figure 4C: Chromatin category ----
c_data <- read_panel('Panel_newC_data_dev2.txt') %>% mutate(
  chromatin_18cate2 = factor(chromatin_18cate2, levels = c('enhancer','quies','reprPC','tss','tx'),
    labels = c('Enhancer','Quiescent','Repressed PC','TSS','Transcription')),
  coloc_stat = factor(coloc_stat, levels = status_levels, labels = status_labels))
tab <- table(c_data$chromatin_18cate2, c_data$coloc_stat)
print(tab)
print(chisq.test(tab))
p <- ggplot(count(c_data, chromatin_18cate2, coloc_stat, name = 'count'),
            aes(chromatin_18cate2, count, fill = coloc_stat)) +
  geom_col(position = position_dodge(width = 0.8), color = 'black', width = 0.7) +
  geom_text(aes(label = count), position = position_dodge(width = 0.8), vjust = -0.3, size = 6, family = 'Arial') +
  scale_fill_manual(values = palette) + scale_y_continuous(limits = c(0,500)) +
  labs(x = 'chromatin category', y = 'count of gHMs', fill = 'haQTL colocalization status') +
  theme_minimal(base_family = 'Arial') + theme(axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text = element_text(size = 18, color = 'black'), axis.title = element_text(size = 20, color = 'black'),
    axis.line = element_line(color = 'black'), panel.grid = element_blank(),
    legend.title = element_text(size = 18), legend.text = element_text(size = 18),
    legend.position = c(0.80, 1.01), legend.justification = c('right','top'),
    legend.background = element_rect(fill = 'white', color = NA))
save_panel(p, 'Figure4C.tiff', 9, 7)

# ---- Figure 4D: Signed distance and effect directions ----
# Exclude palindromic SNPs and restrict absolute distance to <50 kb.
d <- read_panel('panelD_dist_by_effect.txt') %>%
  filter(abs(peak_cpg_dist_kb) < 50, palindromic == 0) %>%
  mutate(direction = factor(same_effect_direct))
print(table(d$direction))
if (!identical(as.integer(table(d$direction)), c(2886L,711L))) warning('Check effect direction coding and manuscript counts.')
p <- ggplot(d, aes(peak_cpg_dist_kb, fill = direction)) +
  geom_histogram(position = 'identity', alpha = 0.5, color = 'grey', bins = 30) +
  xlim(-50,50) + ylim(0,2200) +
  scale_fill_manual(values = palette[1:2], labels = c('opposite direction (2886 pairs)','same direction (711 pairs)')) +
  labs(x = 'distance between gHM and gCpG (kb)', y = 'count of colocalizations', fill = 'Effect directions of haQTL and mQTL') +
  theme_classic(base_family = 'Arial') + theme(plot.margin = margin(1,2,1,1,'cm'),
    legend.title = element_text(size = 18), legend.text = element_text(size = 18),
    legend.position = c(0.5,1), axis.title = element_text(size = 19),
    axis.text = element_text(size = 17, color = 'black'), aspect.ratio = 1,
    legend.background = element_rect(fill = 'white', color = 'white'))
save_panel(p, 'Figure4D.tiff')

# ---- Figure 4E: Effect size correlation ----
# Prepared input contains opposite-direction colocalized pairs.
e <- read_panel('haqtl_mqtl_beta_corr.txt') %>% filter(palindromic == 0)
print(cor.test(e$betahaqtl, e$beta, method = 'pearson'))
print(sum(complete.cases(e[,c('betahaqtl','beta')])))
p <- ggplot(e, aes(betahaqtl, beta)) + geom_point() + geom_smooth(method = 'lm') +
  labs(x = 'beta coefficient of haQTL', y = 'beta coefficient of mQTL') +
  theme_classic(base_family = 'Arial') + theme(axis.title = element_text(size = 22),
    axis.text = element_text(size = 22, color = 'black'), aspect.ratio = 1)
save_panel(p, 'Figure4E.tiff')

# ---- Figure 4F: Nearby measured CpGs within 10 kb ----
f <- read_panel('Panel_newE_data.txt')
# The corrected input must contain counts for the 10-kb window.
# Do not silently relabel counts calculated over a different window.
stopifnot('count_cpg_within10kb_per_peak' %in% names(f))
f$coloc_cate <- factor(as.character(f$coloc_cate), levels = c(status_levels[1:3], '4_no_test', '4_no_coloc_test'),
                       labels = c('strong coloc','suggestive coloc','no coloc','no coloc test','no coloc test'))
stopifnot(!anyNA(f$coloc_cate))
p <- ggplot(f, aes(coloc_cate, count_cpg_within10kb_per_peak, fill = coloc_cate)) +
  geom_violin(trim = FALSE, color = 'black', alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = 'black', fatten = 0) +
  stat_summary(fun = median, geom = 'point', size = 5, color = 'blue', shape = 18) +
  labs(x = NULL, y = 'nearby measured CpGs per gHM (within 10 kb)') +
  scale_fill_manual(values = palette) + theme_classic(base_family = 'Arial') +
  theme(plot.margin = margin(5,2,2,2,'cm'), legend.position = 'none',
    axis.title.y = element_text(size = 18, color = 'black', margin = margin(r = 10)),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.text.x = element_text(size = 18, color = 'black', angle = 45, hjust = 1, vjust = 1))
save_panel(p, 'Figure4F.tiff', 7, 7)
# Two-sided unpaired Student's t-tests against the no-test category.
reference <- f$count_cpg_within10kb_per_peak[f$coloc_cate == 'no coloc test']
for (category in c('strong coloc','suggestive coloc','no coloc')) {
  test <- t.test(f$count_cpg_within10kb_per_peak[f$coloc_cate == category], reference,
                 alternative = 'two.sided', paired = FALSE, var.equal = TRUE)
  cat(category, 'vs no coloc test: P =', format(test$p.value, digits = 8), '\n')
}
