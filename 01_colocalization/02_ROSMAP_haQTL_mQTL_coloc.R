###############################################################################
# Colocalization of ROSMAP haQTLs and ROSMAP mQTLs
#
# This script performs coloc.abf() for 8,962 predefined haQTL peak–CpG pairs.
# The original analysis stored the pair definitions as thousands of individual
# objects (peak1, mqtl1, nh1, nm1, etc.). For readability, those definitions
# are stored in the companion file:
#
#   rosmap_haqtl_mqtl_pairs.tsv
#
# Required input files:
#   chr1MAF.txt ... chr22MAF.txt
#   haqtl_for_coloc_overlapping_snps.RDS
#   mqtl_for_coloc_overlapping_snps.RDS
#   rosmap_haqtl_mqtl_pairs.tsv
#
# Output:
#   rosmap_haqtl_mqtl_coloc_results.txt
###############################################################################

library(coloc)

###############################################################################
# 1. Read minor allele frequencies for autosomal variants
###############################################################################

maf_files <- paste0("chr", 1:22, "MAF.txt")

missing_maf_files <- maf_files[!file.exists(maf_files)]
if (length(missing_maf_files) > 0) {
  stop(
    "Missing MAF file(s): ",
    paste(missing_maf_files, collapse = ", ")
  )
}

maf_list <- lapply(
  maf_files,
  read.table,
  header = FALSE,
  stringsAsFactors = FALSE
)

maf <- do.call(rbind, maf_list)
colnames(maf)[1:2] <- c("snp", "MAF")
maf <- maf[, c("snp", "MAF")]

# Retain the first MAF record if a SNP occurs more than once.
maf <- maf[!duplicated(maf$snp), ]

###############################################################################
# 2. Prepare ROSMAP haQTL summary statistics
###############################################################################

ha <- readRDS("haqtl_for_coloc_overlapping_snps.RDS")

# Rename columns as required by coloc. These positions reproduce the original
# analysis code; verify them if the structure of the input RDS file changes.
colnames(ha)[1] <- "snp"
colnames(ha)[11] <- "pvalues"

ha$varbeta <- ha$se^2
ha$n <- 561
ha$type <- "quant"

# Remove columns not required for colocalization and add MAF.
ha2 <- ha[, !names(ha) %in% c(
  "chr", "pos", "peakPos", "se", "t", "A1", "A2"
)]

ha3 <- merge(ha2, maf, by = "snp")

###############################################################################
# 3. Prepare ROSMAP mQTL summary statistics
###############################################################################

ma <- readRDS("mqtl_for_coloc_overlapping_snps.RDS")

# Rename columns as required by coloc. These positions reproduce the original
# analysis code; verify them if the structure of the input RDS file changes.
colnames(ma)[1] <- "snp"
colnames(ma)[11] <- "pvalues"

ma$varbeta <- ma$se^2
ma$n <- 543
ma$type <- "quant"

# Remove columns not required for colocalization and add MAF.
ma2 <- ma[, !names(ma) %in% c(
  "chr", "pos", "CpGpos", "se", "t", "A1", "A2"
)]

ma3 <- merge(ma2, maf, by = "snp")

###############################################################################
# 4. Read the predefined haQTL peak–CpG pairs
###############################################################################

pair_file <- "rosmap_haqtl_mqtl_pairs.tsv"

if (!file.exists(pair_file)) {
  stop("Pair-definition file not found: ", pair_file)
}

pairs <- read.delim(
  pair_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_pair_columns <- c("pair_id", "peak", "CpG")
missing_pair_columns <- setdiff(required_pair_columns, names(pairs))

if (length(missing_pair_columns) > 0) {
  stop(
    "The pair-definition file is missing column(s): ",
    paste(missing_pair_columns, collapse = ", ")
  )
}

pairs <- pairs[order(pairs$pair_id), ]

if (anyDuplicated(pairs$pair_id)) {
  stop("pair_id must be unique in the pair-definition file.")
}

###############################################################################
# 5. Run coloc.abf() for each haQTL peak–CpG pair
###############################################################################

run_one_coloc <- function(pair_id, peak_id, cpg_id, ha_data, mqtl_data) {

  # Select all haQTL SNP associations for the specified acetylation peak.
  h_temp <- ha_data[ha_data$peak == peak_id, ]

  # Select all mQTL SNP associations for the specified CpG.
  m_temp <- mqtl_data[mqtl_data$CpG == cpg_id, ]

  # Retain the first occurrence of a duplicated SNP, matching the original
  # analysis code.
  h_temp <- h_temp[!duplicated(h_temp$snp), ]
  m_temp <- m_temp[!duplicated(m_temp$snp), ]

  if (nrow(h_temp) == 0 || nrow(m_temp) == 0) {
    warning(
      "No data for pair_id ", pair_id,
      " (peak = ", peak_id, ", CpG = ", cpg_id, ")"
    )
    return(NULL)
  }

  # Construct coloc input objects. Both traits are quantitative.
  h_dataset <- list(
    snp = h_temp$snp,
    beta = h_temp$beta,
    varbeta = h_temp$varbeta,
    pvalues = h_temp$pvalues,
    MAF = h_temp$MAF,
    N = h_temp$n[1],
    type = h_temp$type[1]
  )

  m_dataset <- list(
    snp = m_temp$snp,
    beta = m_temp$beta,
    varbeta = m_temp$varbeta,
    pvalues = m_temp$pvalues,
    MAF = m_temp$MAF,
    N = m_temp$n[1],
    type = m_temp$type[1]
  )

  result <- tryCatch(
    coloc.abf(h_dataset, m_dataset)$summary,
    error = function(e) {
      warning(
        "coloc.abf failed for pair_id ", pair_id,
        " (peak = ", peak_id, ", CpG = ", cpg_id, "): ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(result)) {
    return(NULL)
  }

  # Add pair identifiers to make the output traceable.
  result <- as.data.frame(t(result), check.names = FALSE)
  result$pair_id <- pair_id
  result$peak <- peak_id
  result$CpG <- cpg_id

  result[, c("pair_id", "peak", "CpG", setdiff(names(result), c(
    "pair_id", "peak", "CpG"
  )))]
}

results <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {

  results[[i]] <- run_one_coloc(
    pair_id = pairs$pair_id[i],
    peak_id = pairs$peak[i],
    cpg_id = pairs$CpG[i],
    ha_data = ha3,
    mqtl_data = ma3
  )

  if (i %% 500 == 0 || i == nrow(pairs)) {
    message("Completed ", i, " of ", nrow(pairs), " colocalization tests")
  }
}

###############################################################################
# 6. Combine and save results
###############################################################################

successful_results <- results[!vapply(results, is.null, logical(1))]

if (length(successful_results) == 0) {
  stop("No colocalization tests completed successfully.")
}

coloc_results <- do.call(rbind, successful_results)
rownames(coloc_results) <- NULL

write.table(
  coloc_results,
  file = "rosmap_haqtl_mqtl_coloc_results.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "Finished. ", nrow(coloc_results),
  " successful results were written to rosmap_haqtl_mqtl_coloc_results.txt"
)
