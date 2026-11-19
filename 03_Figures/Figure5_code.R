# Figure 5: ROSMAP DLPFC H3K9ac haQTLs and NIMH WGBS mQTLs
# Panels A-F; panel G is not included.
# Run from a directory containing the named input files.
library(data.table)
library(ggplot2)
colors <- c("#aaaaff", "#33dd33", "#ffaa00", "#aaaaaa")
status_labels <- c("strong colocalization (PP4≥0.8)", "suggestive colocalization (0.5≤PP4<0.8)", "no colocalization (PP4<0.5)", "no colocalization test")
status_levels <- c("1_strong_coloc", "2_suggestive_coloc", "3_no_coloc", "4_no_test")
all_coloc <- fread("coloc_abf_chr1_22_all_results.tsv")
coloc <- all_coloc[n_observed > 0 & n_common >= 50]
peak_pp4 <- coloc[!is.na(PP.H4), .(max_PP4=max(PP.H4)), by=peak]
peak_pp4[, status := fcase(max_PP4 >= .8, status_levels[1], max_PP4 >= .5, status_levels[2], default=status_levels[3])]
# Panel A: derive counts from latest results, compare to manuscript.
counts <- c(as.integer(table(factor(peak_pp4$status, levels=status_levels[1:3]))), 2369L-nrow(peak_pp4))
if (any(counts < 0)) stop("More tested peaks than manuscript total")
if (!identical(as.integer(counts), c(1127L,61L,74L,1107L))) warning("Panel A counts differ from manuscript figure; verify input version")
a <- data.frame(group=factor(status_levels, levels=status_levels), value=counts)
pa <- ggplot(a,aes(x="",y=value,fill=group))+geom_col(width=1,color="white")+coord_polar(theta="y")+scale_fill_manual(values=colors,labels=status_labels)+labs(fill="haQTL colocalization status")+theme_void(base_family="Arial")+theme(legend.position="top",legend.direction="vertical",legend.title=element_text(size=16),legend.text=element_text(size=14))+guides(fill=guide_legend(ncol=1))
ggsave("PanelA_pie_chart.tiff",pa,width=12,height=5,units="in",dpi=300,compression="lzw")
# Panel B: count unique strong CpGs per tested peak, including zero.
dat <- coloc[, .(
  n_PP4_ge_p8 = uniqueN(cpg[!is.na(PP.H4) & PP.H4 >= 0.8])
), by = peak]

# Plotting position:
# n_PP4 = 0      -> bar centered at -10  (covers -20 to 0)
# n_PP4 = 1-20   -> bar centered at 10   (covers 0 to 20)
# n_PP4 = 21-40  -> bar centered at 30   (covers 20 to 40)
# etc.
dat[, x_plot := fifelse(
  n_PP4_ge_p8 == 0,
  -10,
  ((n_PP4_ge_p8 - 1) %/% 20) * 20 + 10
)]

# Count number of peaks in each bin
plot_dat <- dat[, .(
  n_peaks = .N
), by = x_plot]

setorder(plot_dat, x_plot)

# Check counts
print(plot_dat)

cat("Total unique haQTL peaks:", nrow(dat), "\n")
cat("Peaks with 0 colocalized CpGs:",
    sum(dat$n_PP4_ge_p8 == 0), "\n")
cat("Peaks with >=1 colocalized CpG:",
    sum(dat$n_PP4_ge_p8 > 0), "\n")

# Plot
p <- ggplot(plot_dat, aes(x = x_plot, y = n_peaks)) +
  geom_col(
    fill = "#EEEE00",
    color = "black",
    width = 20
  ) +
  scale_x_continuous(
    breaks = seq(0, max(dat$n_PP4_ge_p8) + 100, by = 100),
    limits = c(-20, max(dat$n_PP4_ge_p8) + 20),
    expand = c(0, 0)
  ) +
  labs(
    x = "# colocalized CpGs per haQTL",
    y = "count of haQTL peaks"
  ) +
  theme_classic(base_size = 16, base_family = "Arial") +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

print(p)

ggsave(
  "ssimp_panelB_hist.tiff",
  p,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

# Panel C: chromatin categories (final annotated peak-level input).
cdata <- fread("Panel_newC_data_dev2.txt")
if (!"coloc_status" %in% names(cdata)) cdata[, coloc_status := coloc_stat]
cdata[, coloc_status := factor(sub("4_no_coloc_test", "4_no_test", coloc_status), levels=status_levels)]
cdata[, chromatin_18cate2 := factor(chromatin_18cate2, levels=c("enhancer","quies","reprPC","tss","tx"), labels=c("Enhancer","Quiescent","Repressed PC","TSS","Transcription"))]
print(chisq.test(table(cdata$chromatin_18cate2,cdata$coloc_status)))
cs <- cdata[, .(count=.N), by=.(chromatin_18cate2,coloc_status)]
pc <- ggplot(cs,aes(chromatin_18cate2,count,fill=coloc_status))+geom_col(position=position_dodge(.9),color="black",width=.9)+geom_text(aes(label=count),position=position_dodge(.9),vjust=-.3,size=5.5)+scale_y_continuous(limits=c(0,700))+scale_fill_manual(values=colors,labels=status_labels,drop=FALSE)+labs(x="chromatin category",y="count of gHMs",fill="haQTL colocalization status")+guides(fill=guide_legend(ncol=1))+theme_classic(base_family="Arial")+theme(plot.margin=margin(1,1,1,1,"cm"),legend.position="top",legend.text=element_text(size=18),legend.title=element_text(size=18),axis.title=element_text(size=18),axis.text.y=element_text(size=18,color="black"),axis.text.x=element_text(size=14,color="black",angle=20,hjust=1))
ggsave("PanelC_chromatin.tiff",pc,width=9,height=8,units="in",dpi=300,compression="lzw")
# Panel D: exact peak/lead SNP and CpG/lead SNP matching.
# The distance formula retains the user-confirmed latest code.
# IMPORTANT: haqtl_peakPos must be in hg38 to match mqtl_cpgPos_hg38.
coloc2 <- coloc[
  !is.na(PP.H4) & PP.H4 >= 0.8,
  .(peak, cpg, lead_haSNP, PP.H4)
]

# Remove spaces before creating combined variables
coloc2[, peak_leadhaSNP := paste0(
  gsub("\\s+", "", peak),
  "_",
  gsub("\\s+", "", lead_haSNP)
)]

coloc2[, cpg_haSNP := paste0(
  gsub("\\s+", "", cpg),
  "_",
  gsub("\\s+", "", lead_haSNP)
)]

# Preselected observed WGBS mQTL data, including CpG-SNP key.
mqtl <- readRDS("WGBS_mQTL_selected_panelD.rds")
setDT(mqtl)
mqtl <- mqtl[gsub("\\s+", "", cpg) %chin% unique(gsub("\\s+", "", coloc2$cpg))]
ha_all <- readRDS("haQTL_allchr_1Mb.rds")


###############merging###############


## Make sure all are data.table
setDT(coloc2)
setDT(ha_all)
setDT(mqtl)


############################################################
## 1. Prepare haQTL data
############################################################

## Remove t because we don't need it
if ("t" %in% names(ha_all)) ha_all[, t := NULL]

## Rename haQTL-specific variables
setnames(
  ha_all,
  old = c(
    "rsID", "chr", "pos",
    "A1", "A2",
    "peakPos",
    "beta", "se", "p"
  ),
  new = c(
    "haqtl_snp", "haqtl_chr", "haqtl_snpPos",
    "haqtl_A1", "haqtl_A2",
    "haqtl_peakPos",
    "haqtl_beta", "haqtl_se", "haqtl_p"
  )
)


############################################################
## 2. Link coloc2 with haQTL using peak + lead haQTL SNP
############################################################

## Clean matching keys just in case there are spaces
coloc2[, peak_leadhaSNP := gsub("\\s+", "", peak_leadhaSNP)]
ha_all[, peak_snp := gsub("\\s+", "", peak_snp)]

cat("Before haQTL merge:\n")
cat("  coloc rows:", nrow(coloc2), "\n")
cat("  unique peaks:", uniqueN(coloc2$peak), "\n\n")

## Inner merge = retain matched rows only
coloc_ha <- merge(
  coloc2,
  ha_all,
  by.x = "peak_leadhaSNP",
  by.y = "peak_snp",
  all = FALSE
)

cat("After haQTL merge:\n")
cat("  rows:", nrow(coloc_ha), "\n")
cat("  unique peaks:", uniqueN(coloc_ha$peak.x), "\n\n")


############################################################
## 3. Clean up duplicated peak variable from merge
############################################################

## peak.x is the coloc peak
## peak.y is the haQTL peak
## They should be identical

cat("Peak mismatches:",
    sum(coloc_ha$peak.x != coloc_ha$peak.y, na.rm = TRUE),
    "\n")

## Keep one peak variable
coloc_ha[, peak := peak.x]
coloc_ha[, c("peak.x", "peak.y") := NULL]


############################################################
## 4. Prepare mQTL data
############################################################

setnames(
  mqtl,
  old = c(
    "methChr", "methPos_hg38",
    "snpRsNum", "snpChr_hg38", "snpPos_hg38",
    "snpCounted", "snpAlt",
    "beta", "pvalue"
  ),
  new = c(
    "mqtl_cpgChr", "mqtl_cpgPos_hg38",
    "mqtl_snp", "mqtl_snpChr_hg38", "mqtl_snpPos_hg38",
    "mqtl_snpCounted", "mqtl_snpAlt",
    "mqtl_beta", "mqtl_p"
  )
)

## Clean mQTL matching key
mqtl[, cpg_snpRsNum := gsub("\\s+", "", cpg_snpRsNum)]
coloc_ha[, cpg_haSNP := gsub("\\s+", "", cpg_haSNP)]


############################################################
## 5. Link coloc-haQTL data with mQTL
##    using CpG + lead haQTL SNP
############################################################

panelD <- merge(
  coloc_ha,
  mqtl,
  by.x = "cpg_haSNP",
  by.y = "cpg_snpRsNum",
  all = FALSE
)

cat("\nAfter mQTL merge:\n")
cat("  rows:", nrow(panelD), "\n")
cat("  unique peaks:", uniqueN(panelD$peak), "\n")
cat("  unique CpGs:", uniqueN(panelD$cpg.x), "\n")


############################################################
## 6. Check CpG agreement
############################################################

cat("CpG mismatches:",
    sum(panelD$cpg.x != panelD$cpg.y, na.rm = TRUE),
    "\n")

## Keep one cpg variable
panelD[, cpg := cpg.x]
panelD[, c("cpg.x", "cpg.y") := NULL]


############################################################
## Final check
############################################################

cat("\nFINAL PANEL D DATA:\n")
cat("Rows:", nrow(panelD), "\n")
cat("Unique peaks:", uniqueN(panelD$peak), "\n")
cat("Unique CpGs:", uniqueN(panelD$cpg), "\n")


##########################################################################################################
# Remove palindromic SNPs (A/T, T/A, C/G, G/C).
# Harmonize mQTL beta to the haQTL A1 effect allele:- haqtl_A1 == mqtl_snpCounted and haqtl_A2 == mqtl_snpAlt → keep mQTL beta.
# - reversed → multiply mQTL beta by −1.
# - anything else → remove.
# 
# Calculate signed peak-to-CpG distance: (CpG - peak)/1000.
# Classify haQTL/mQTL effects as same vs opposite direction and print frequencies.
##########################################################################################################


############################################################
## 1. Standardize alleles
############################################################

panelD[, `:=`(
  haqtl_A1        = toupper(trimws(haqtl_A1)),
  haqtl_A2        = toupper(trimws(haqtl_A2)),
  mqtl_snpCounted = toupper(trimws(mqtl_snpCounted)),
  mqtl_snpAlt     = toupper(trimws(mqtl_snpAlt))
)]


############################################################
## 2. Remove palindromic SNPs
##    A/T, T/A, C/G, G/C
############################################################

panelD[, haqtl_alleles := paste0(haqtl_A1, haqtl_A2)]

cat("Before removing palindromic SNPs:\n")
cat("Rows:", nrow(panelD), "\n")
cat("Unique peaks:", uniqueN(panelD$peak), "\n\n")

panelD2 <- panelD[
  !haqtl_alleles %in% c("AT", "TA", "CG", "GC")
]

cat("After removing palindromic SNPs:\n")
cat("Rows:", nrow(panelD2), "\n")
cat("Unique peaks:", uniqueN(panelD2$peak), "\n\n")


############################################################
## 3. Determine allele orientation
##
## haQTL A1 = effect allele
## mQTL snpCounted = effect allele
############################################################

panelD2[, allele_match := fcase(
  
  ## Same effect allele orientation
  haqtl_A1 == mqtl_snpCounted &
    haqtl_A2 == mqtl_snpAlt,
  "same",
  
  ## Effect alleles reversed
  haqtl_A1 == mqtl_snpAlt &
    haqtl_A2 == mqtl_snpCounted,
  "reverse",
  
  default = "mismatch"
)]

cat("Allele matching:\n")
print(table(panelD2$allele_match, useNA = "ifany"))


############################################################
## 4. Remove unresolved allele mismatches
############################################################

panelD2 <- panelD2[allele_match != "mismatch"]

cat("\nAfter removing allele mismatches:\n")
cat("Rows:", nrow(panelD2), "\n")
cat("Unique peaks:", uniqueN(panelD2$peak), "\n")
cat("Unique CpGs:", uniqueN(panelD2$cpg), "\n\n")


############################################################
## 5. Harmonize mQTL beta to haQTL A1 effect allele
############################################################

panelD2[, mqtl_beta_harmonized :=
          fifelse(
            allele_match == "same",
            mqtl_beta,
            -mqtl_beta
          )
]


############################################################
## 6. Calculate signed peak-to-CpG distance
##
## positive = CpG is downstream/right of peak
## negative = CpG is upstream/left of peak
############################################################

panelD2[, peak_cpg_distance_kb :=
          (mqtl_cpgPos_hg38 - haqtl_peakPos) / 1000
]


############################################################
## 7. Classify effect direction
############################################################

panelD2[, effect_direction := fcase(
  
  haqtl_beta > 0 & mqtl_beta_harmonized > 0,
  "same",
  
  haqtl_beta < 0 & mqtl_beta_harmonized < 0,
  "same",
  
  haqtl_beta > 0 & mqtl_beta_harmonized < 0,
  "opposite",
  
  haqtl_beta < 0 & mqtl_beta_harmonized > 0,
  "opposite",
  
  default = NA_character_
)]


############################################################
## 8. Frequency of effect directions
############################################################

cat("\nEffect direction counts:\n")

print(
  table(
    panelD2$effect_direction,
    useNA = "ifany"
  )
)

cat("\nEffect direction percentages:\n")

print(
  round(
    prop.table(
      table(panelD2$effect_direction)
    ) * 100,
    2
  )
)


## Use current Panel D data
down <- copy(panelD2[!is.na(effect_direction) & !is.na(peak_cpg_distance_kb)])

## Factor order: opposite first, same second
down$same_direct_haqtl_mqtl <- factor(
  down$effect_direction,
  levels = c("opposite", "same")
)

tiff(
  "PanelD_gare_cpg_dist.tiff",
  width = 6,
  height = 6,
  units = "in",
  res = 300
)

print(ggplot(
  down,
  aes(
    x = peak_cpg_distance_kb,
    fill = same_direct_haqtl_mqtl
  )
) +
  geom_histogram(
    position = "identity",
    alpha = 0.5,
    color = "grey"
  ) +
  xlim(-50, 50) +
  theme(
    plot.margin = margin(0.3, 0.5, 0.3, 0.3, "cm"),
    panel.background = element_rect(fill = "white", color = "white"),
    
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.direction = "vertical",
    
    legend.title = element_text(
      size = 15,
      color = "black",
      family = "Arial",
      hjust = 0
    ),
    legend.text = element_text(
      size = 14,
      color = "black",
      family = "Arial"
    ),
    
    ## Thicker black x- and y-axis lines
    axis.line = element_line(
      color = "black",
      linewidth = 0.8
    ),
    
    ## Thicker black tick marks
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.8
    ),
    
    ## Longer tick marks
    axis.ticks.length = unit(0.15, "cm"),
    
    axis.title.x = element_text(
      size = 20,
      color = "black",
      family = "Arial"
    ),
    axis.title.y = element_text(
      size = 20,
      color = "black",
      family = "Arial"
    ),
    axis.text.y = element_text(
      size = 15,
      color = "black",
      family = "Arial"
    ),
    axis.text.x = element_text(
      size = 15,
      color = "black",
      family = "Arial"
    ),
    
    aspect.ratio = 1,
    legend.background = element_rect(
      fill = "white",
      color = "white"
    )
  ) +
  labs(
    x = "distance between gHM and gCpG (kb)",
    y = "count of colocalizations"
  ) +
  scale_fill_manual(
    values = c("#aaaaff", "#33dd33"),
    labels = c(
      "opposite direction (32,426)",
      "same direction (20,135)"
    )
  ) +
  guides(
    fill = guide_legend(
      title = "effect directions of haQTL and mQTL",
      ncol = 1,
      byrow = TRUE
    )
  )

)
dev.off()



# Panel E: opposite-direction effect sizes.
## Restrict to opposite-direction colocalizations
data <- panelD2[panelD2$effect_direction == "opposite", ]

## Panel E beta values
data$betahaqtl <- abs(data$haqtl_beta)
data$beta <- -abs(data$mqtl_beta_harmonized)

## Check
cat("Number of opposite-direction pairs:", nrow(data), "\n")
summary(data$betahaqtl)
summary(data$beta)

tiff("PanelE_beta_corr.tiff", width = 6, height = 6, units = "in", res = 300)

print(ggplot(data, aes(x = betahaqtl, y = beta)) + 
  geom_point() +
  scale_x_continuous(limits = c(0.26, NA)) +
  geom_smooth(method = lm) +
  labs(
    x = "beta coefficient of haQTL",
    y = "beta coefficient of mQTL"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.title.x = element_text(
      size = 22, color = "black", family = "Arial"
    ),
    axis.title.y = element_text(
      size = 22, color = "black", family = "Arial"
    ),
    axis.text.y = element_text(
      size = 22, color = "black", family = "Arial"
    ),
    axis.text.x = element_text(
      size = 22, color = "black", family = "Arial"
    ),
    aspect.ratio = 1 / 1,
    legend.background = element_rect(fill = "white", color = "white")
  )

)
dev.off()

## Pearson correlation
cor.test(
  data$betahaqtl,
  data$beta,
  method = "pearson"
)
# Panel F: final measured CpG coverage input.
fdata <- fread("violin_plot_data.txt")
fdata[, coloc_cate := factor(coloc_cate, levels=status_levels)]
pf <- ggplot(fdata,aes(coloc_cate,count_cpg_within_10kbp,fill=coloc_cate))+geom_violin(trim=FALSE,color="black",alpha=.7)+geom_boxplot(width=.1,outlier.shape=NA,color="black",fatten=0)+stat_summary(fun=median,geom="point",size=5,color="blue",shape=18)+labs(x=NULL,y="nearby measured CpGs per gHM (within 10 kb)")+scale_fill_manual(values=colors,labels=status_labels)+scale_x_discrete(labels=c("strong coloc","suggestive coloc","no coloc","no coloc test"))+theme_classic(base_family="Arial")+theme(plot.margin=margin(2,2,2,2,"cm"),legend.position="none",axis.title.y=element_text(size=18,color="black",face="plain",margin=margin(r=10)),axis.text.y=element_text(size=18,color="black",face="plain"),axis.text.x=element_text(size=18,color="black",face="plain",angle=45,hjust=1,vjust=1))
ggsave("Panel_violin.tiff",pf,width=6,height=7,units="in",dpi=300,compression="lzw")
for (group in status_levels[1:3]) {
  test <- t.test(fdata[coloc_cate==group,count_cpg_within_10kbp], fdata[coloc_cate==status_levels[4],count_cpg_within_10kbp], alternative="two.sided",var.equal=TRUE)
  cat(group,"vs no test: P =",test$p.value,"\n")
}
# Panel G intentionally omitted.
