# ROSMAP DLPFC H3K9ac haQTL and NIMH WGBS mQTL colocalization
# Primary coloc.abf analysis using SSimp P-values across autosomes.
# Required packages: data.table, coloc.
# Required inputs: preprocessed haQTL RDS, hg19 peak-CpG pair table,
# and chromosome-specific SSimp TSV parts (generated separately).
# No imputation or genome-build conversion is performed here.
# Run from the input directory; outputs are written to that directory.

library(data.table)
library(coloc)

# ============================================================
# SETTINGS
# ============================================================

# Run from the directory containing the input files.
# Usage: Rscript 03_ROSMAP_haQTL_NIMH_mQTL_coloc.R [chromosome ...]
# No arguments: run all autosomes. Example argument: 17
args <- commandArgs(trailingOnly = TRUE)
chromosomes <- if (length(args)) as.integer(args) else 1:22
stopifnot(!anyNA(chromosomes), all(chromosomes %in% 1:22),
          !anyDuplicated(chromosomes))

cat("R version:", R.version.string, "\n")
cat("coloc version:", as.character(packageVersion("coloc")), "\n")

# ============================================================
# FILES
# ============================================================

ha_file   <- "haqtl_for_WGBS_coloc_wMAF.rds"
pair_file <- "ROSMAP_haQTL_WGBS_pairs_for_SSimp_coloc_hg19.txt"

stopifnot(file.exists(ha_file), file.exists(pair_file))

cat("Reading haQTL data...\n")

ha <- readRDS(ha_file)

# IMPORTANT:
# readRDS returns ha as a data.frame.
# Convert to data.table before using data.table syntax.
setDT(ha)

cat("ha class:", class(ha), "\n")
cat("ha rows:", nrow(ha), "\n")

cat("Reading peak-CpG pairs...\n")

pairs <- fread(pair_file)
ha <- ha[
    !is.na(snp) &
    snp != "" &
    !is.na(peak) &
    is.finite(beta) &
    is.finite(varbeta) &
    varbeta > 0 &
    is.finite(MAF) &
    MAF > 0 &
    MAF <= 0.5
]

# One SNP per haQTL peak
ha <- ha[!duplicated(ha[, .(peak, snp)])]

# Process one chromosome at a time to limit memory usage.
run_chromosome <- function(CHR) {
    pairs_chr <- pairs[chr == CHR]
    if (!nrow(pairs_chr)) {
        message("No pairs for chromosome ", CHR)
        return(invisible(NULL))
    }

    # Find chromosome-specific SSimp output files. An optional filename
    # label (e.g., a reference-panel label) may occur before "_part".
    # The delimiter after the chromosome prevents chr1 matching chr10-19.
    imp_files <- list.files(
        pattern = sprintf("^SSimp_WGBS_chr%d(_[^.]+)?_part[0-9]+[.]tsv$", CHR),
        full.names = TRUE
    )
    if (!length(imp_files)) stop("No SSimp input files for chromosome ", CHR)

    # Do not combine alternative imputation runs for the same chromosome.
    prefixes <- sub("_part[0-9]+[.]tsv$", "", basename(imp_files))
    if (length(unique(prefixes)) != 1L) {
        stop("Multiple SSimp file sets for chromosome ", CHR,
             ". Keep only the intended analysis files in the input directory.")
    }
    # Read in numeric part order to preserve first-record duplicate handling.
    part_ids <- as.integer(sub("^.*_part([0-9]+)[.]tsv$", "\\1", basename(imp_files)))
    if (anyDuplicated(part_ids) || !setequal(part_ids, seq_len(max(part_ids)))) {
        stop("Duplicate or missing intermediate SSimp parts for chromosome ", CHR)
    }
    imp_files <- imp_files[order(part_ids)]
    # Verify the final part count against the SSimp job outputs before running.
    # Filename discovery alone cannot detect a missing final part.
    print(imp_files)
    mq_chr <- rbindlist(lapply(imp_files, fread), use.names = TRUE, fill = TRUE)

    mq_chr <- mq_chr[
        !is.na(cpg) &
        !is.na(SNP) &
        SNP != "" &
        is.finite(P.imp) &
        P.imp > 0 &
        P.imp <= 1 &
        is.finite(maf) &
        maf > 0 &
        maf <= 0.5
    ]

    # One SNP per CpG
    mq_chr <- mq_chr[!duplicated(mq_chr[, .(cpg, SNP)])]

    cat("Clean haQTL rows:", nrow(ha), "\n")
    cat("Clean SSimp rows:", nrow(mq_chr), "\n")

    # ============================================================
    # RESULT CONTAINERS
    # ============================================================

    res_list  <- vector("list", nrow(pairs_chr))
    fail_list <- vector("list", nrow(pairs_chr))

    # ============================================================
    # COLOC LOOP
    # ============================================================

    for (i in seq_len(nrow(pairs_chr))) {

        peak_i <- pairs_chr$peak[i]
        cpg_i  <- pairs_chr$cpg[i]
        lead_i <- pairs_chr$lead_haSNP[i]

        if (i %% 100 == 0L || i == 1L) {
            cat(
                "Pair", i, "/", nrow(pairs_chr),
                "|", peak_i,
                "|", cpg_i,
                "\n"
            )
        }

        # --------------------------------------------------------
        # Extract haQTL peak and WGBS CpG
        # --------------------------------------------------------

        h <- ha[peak == peak_i]
        m <- mq_chr[cpg == cpg_i]

        if (nrow(h) == 0L) {

            fail_list[[i]] <- data.table(
                pair_id = pairs_chr$pair_id[i],
                peak = peak_i,
                cpg = cpg_i,
                reason = "no_haQTL_rows"
            )

            next
        }

        if (nrow(m) == 0L) {

            fail_list[[i]] <- data.table(
                pair_id = pairs_chr$pair_id[i],
                peak = peak_i,
                cpg = cpg_i,
                reason = "no_SSimp_rows"
            )

            next
        }

        # --------------------------------------------------------
        # Common SNPs
        # --------------------------------------------------------

        common <- intersect(h$snp, m$SNP)

        if (length(common) < 2L) {

            fail_list[[i]] <- data.table(
                pair_id = pairs_chr$pair_id[i],
                peak = peak_i,
                cpg = cpg_i,
                reason = "fewer_than_2_common_SNPs",
                n_common = length(common)
            )

            next
        }

        # Put both datasets in identical SNP order
        h2 <- h[match(common, snp)]
        m2 <- m[match(common, SNP)]

        stopifnot(all(h2$snp == m2$SNP))

        # --------------------------------------------------------
        # Construct coloc datasets
        #
        # haQTL:
        # actual beta + varbeta
        #
        # WGBS:
        # SSimp P.imp + MAF + original sample size N=165
        # --------------------------------------------------------

        d1 <- list(
            beta    = as.numeric(h2$beta),
            varbeta = as.numeric(h2$varbeta),
            snp     = as.character(h2$snp),
            MAF     = as.numeric(h2$MAF),
            N       = 561,
            type    = "quant"
        )

        d2 <- list(
            pvalues = as.numeric(m2$P.imp),
            snp     = as.character(m2$SNP),
            MAF     = as.numeric(m2$maf),
            N       = 165,
            type    = "quant"
        )

        # --------------------------------------------------------
        # Run coloc.abf
        # --------------------------------------------------------

        cc <- tryCatch(
            coloc.abf(
                dataset1 = d1,
                dataset2 = d2
            ),
            error = function(e) e
        )

        if (inherits(cc, "error")) {

            fail_list[[i]] <- data.table(
                pair_id = pairs_chr$pair_id[i],
                peak = peak_i,
                cpg = cpg_i,
                reason = paste0(
                    "coloc_error: ",
                    conditionMessage(cc)
                ),
                n_common = length(common)
            )

            next
        }

        # --------------------------------------------------------
        # Extract coloc summary
        # --------------------------------------------------------

        sm <- cc$summary

        # --------------------------------------------------------
        # SNP with highest SNP.PP.H4
        # --------------------------------------------------------

        top_snp <- NA_character_
        top_snp_pp_h4 <- NA_real_

        if (
            !is.null(cc$results) &&
            "SNP.PP.H4" %in% names(cc$results) &&
            nrow(cc$results) > 0L
        ) {

            jj <- which.max(cc$results$SNP.PP.H4)

            top_snp <- as.character(
                cc$results$snp[jj]
            )

            top_snp_pp_h4 <- as.numeric(
                cc$results$SNP.PP.H4[jj]
            )
        }

        # --------------------------------------------------------
        # QC information
        # --------------------------------------------------------

        n_observed <- sum(
            m2$snp_status == "observed",
            na.rm = TRUE
        )

        n_imputed <- sum(
            m2$snp_status == "imputed",
            na.rm = TRUE
        )

        min_mqtl_p <- min(
            m2$P.imp,
            na.rm = TRUE
        )

        min_ha_p <- min(
            h2$pvalues,
            na.rm = TRUE
        )

        lead_in_common <- lead_i %in% common

        # --------------------------------------------------------
        # Save result
        # --------------------------------------------------------

        res_list[[i]] <- data.table(
            pair_id = pairs_chr$pair_id[i],
            peak = peak_i,
            cpg = cpg_i,
            lead_haSNP = lead_i,

            n_common = length(common),
            n_observed = n_observed,
            n_imputed = n_imputed,

            lead_haSNP_in_common = lead_in_common,

            min_ha_p = min_ha_p,
            min_mqtl_p = min_mqtl_p,

            PP.H0 = as.numeric(sm["PP.H0.abf"]),
            PP.H1 = as.numeric(sm["PP.H1.abf"]),
            PP.H2 = as.numeric(sm["PP.H2.abf"]),
            PP.H3 = as.numeric(sm["PP.H3.abf"]),
            PP.H4 = as.numeric(sm["PP.H4.abf"]),

            top_H4_snp = top_snp,
            top_H4_snp_PP = top_snp_pp_h4
        )
    }

    # ============================================================
    # COMBINE RESULTS
    # ============================================================

    res <- rbindlist(
        res_list,
        use.names = TRUE,
        fill = TRUE
    )

    failed <- rbindlist(
        fail_list,
        use.names = TRUE,
        fill = TRUE
    )

    # ============================================================
    # ADD VERSION INFORMATION
    # ============================================================

    if (nrow(res) > 0L) {
        res[, R_version := R.version.string]
        res[, coloc_version := as.character(packageVersion("coloc"))]
    }

    if (nrow(failed) > 0L) {
        failed[, R_version := R.version.string]
        failed[, coloc_version := as.character(packageVersion("coloc"))]
    }

    # ============================================================
    # WRITE OUTPUT
    # ============================================================

    if (!ncol(failed)) failed <- data.table(pair_id = character(), peak = character(), cpg = character(), reason = character())
    if (!ncol(res)) res <- data.table(pair_id = character(), peak = character(), cpg = character(), n_common = integer(), PP.H4 = numeric())

    fwrite(
        res,
        sprintf("coloc_abf_chr%d_SSimp_results.tsv", CHR),
        sep = "\t",
        quote = FALSE,
        na = "NA"
    )

    fwrite(
        failed,
        sprintf("coloc_abf_chr%d_SSimp_failed.tsv", CHR),
        sep = "\t",
        quote = FALSE,
        na = "NA"
    )

    # ============================================================
    # SUMMARY
    # ============================================================

    cat("\n========================================\n")
    cat("chr", CHR, " coloc finished\n", sep = "")
    cat("========================================\n")

    cat("Total peak-CpG pairs:", nrow(pairs_chr), "\n")
    cat("Successful coloc:", nrow(res), "\n")
    cat("Failed/skipped:", nrow(failed), "\n")

    if (nrow(res) > 0L) {

        cat(
            "PP4 >= 0.8:",
            sum(res$PP.H4 >= 0.8, na.rm = TRUE),
            "\n"
        )

        cat(
            "PP4 >= 0.5:",
            sum(res$PP.H4 >= 0.5, na.rm = TRUE),
            "\n"
        )

        cat(
            "Median common SNPs:",
            median(res$n_common, na.rm = TRUE),
            "\n"
        )

        cat(
            "Median observed SNPs:",
            median(res$n_observed, na.rm = TRUE),
            "\n"
        )

        cat(
            "Median imputed SNPs:",
            median(res$n_imputed, na.rm = TRUE),
            "\n"
        )
    }

    if (nrow(failed) > 0L) {

        cat("\nFailure reasons:\n")
        print(failed[, .N, by = reason][order(-N)])
    }

    cat("\nOutput files:\n")
    cat(sprintf("coloc_abf_chr%d_SSimp_results.tsv\n", CHR))
    cat(sprintf("coloc_abf_chr%d_SSimp_failed.tsv\n", CHR))
}

# Run sequentially to avoid loading all chromosome SSimp files into memory.
for (CHR in chromosomes) {
    run_chromosome(CHR)
}
