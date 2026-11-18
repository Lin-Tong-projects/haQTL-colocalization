# NIMH WGBS mQTL summary statistics imputation (autosomes 1–22)
# Requires R/data.table, SSimp on PATH, and indexed reference VCFs.
# Keep the same N_TASKS and input pair ordering when restarting a run.
# Use a fresh output directory or remove old outputs/checkpoints when changing
# reference files, inputs, or task count. An interrupted append before checkpoint
# recording can leave duplicate rows; inspect outputs before downstream analysis.
# This script generalizes the chr21 implementation; use the original reference
# actually used for each chromosome via SSIMP_REF_MANIFEST if needed.

library(data.table)

# ============================================================
# SETTINGS
# ============================================================

# Run from the directory containing the input summary statistics.
# Usage: Rscript NIMH_WGBS_SSimp.R 21
# Optional: provide multiple chromosomes; defaults to all autosomes.
args <- commandArgs(trailingOnly = TRUE)
chromosomes <- if (length(args)) suppressWarnings(as.integer(args)) else 1:22
if (anyNA(chromosomes) || any(!chromosomes %in% 1:22)) {
    stop("Specify chromosome numbers between 1 and 22")
}
chromosomes <- unique(chromosomes)

pair_file <- "ROSMAP_haQTL_WGBS_pairs_for_SSimp_coloc_hg19.txt"
mqtl_file <- "WGBS_mQTL_selected_coloc_cpgs.csv.gz"
sample_file <- "ssimp/SSimp_ASW61_EUR29.txt"

# Set SSIMP_REF_DIR to the directory containing indexed Phase 3 VCFs.
ref_dir <- Sys.getenv("SSIMP_REF_DIR", unset = "reference_panel")

# Optional tab-separated manifest with columns chr and ref_vcf, for
# chromosome-specific reference filenames (e.g. a custom chr17 reference).
manifest_file <- Sys.getenv("SSIMP_REF_MANIFEST", unset = "")
manifest <- if (nzchar(manifest_file)) fread(manifest_file) else NULL
if (!is.null(manifest)) {
    stopifnot(all(c("chr", "ref_vcf") %in% names(manifest)))
    stopifnot(!anyDuplicated(manifest$chr))
}

ssimp_path <- Sys.which("ssimp")

if (ssimp_path == "") {
    stop("ssimp is not available in PATH")
}

# ============================================================
# SLURM ARRAY INFORMATION
# ============================================================

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
n_tasks <- as.integer(Sys.getenv("N_TASKS", unset = "1"))

if (is.na(task_id) || is.na(n_tasks) || n_tasks < 1L ||
    task_id < 1L || task_id > n_tasks) {
    stop("Invalid array task ID or N_TASKS")
}

# Read shared inputs once per R process, not once per chromosome.
stopifnot(file.exists(pair_file), file.exists(mqtl_file), file.exists(sample_file))
pairs <- fread(pair_file)
mqtl <- fread(mqtl_file)

# Each array task processes its own subset of CpGs for each chromosome.
run_chromosome <- function(CHR) {
    ref_vcf <- file.path(ref_dir, sprintf(
        "ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", CHR
    ))
    if (!is.null(manifest) && CHR %in% manifest$chr) {
        ref_vcf <- manifest$ref_vcf[match(CHR, manifest$chr)]
    }

    cat(
        "chr", CHR,
        " | task ", task_id,
        " / ", n_tasks,
        "\n",
        sep = ""
    )

    # ============================================================
    # OUTPUT / TEMP FILES FOR THIS TASK
    # ============================================================

    part_tag <- sprintf("%02d", task_id)

    tmp_dir <- file.path(
        sprintf("ssimp_chr%d_tmp", CHR),
        paste0("part_", part_tag)
    )

    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

    combined_file <- paste0(
        sprintf("SSimp_WGBS_chr%d_part", CHR),
        part_tag,
        ".tsv"
    )

    completed_file <- paste0(
        sprintf("SSimp_chr%d_completed_part", CHR),
        part_tag,
        ".txt"
    )

    failed_file <- paste0(
        sprintf("SSimp_chr%d_failed_part", CHR),
        part_tag,
        ".txt"
    )

    # ============================================================
    # CHECK REQUIRED FILES
    # ============================================================

    stopifnot(file.exists(pair_file))
    stopifnot(file.exists(mqtl_file))
    stopifnot(file.exists(sample_file))
    stopifnot(file.exists(ref_vcf))
    stopifnot(file.exists(paste0(ref_vcf, ".tbi")))

    # ============================================================
    # CREATE CHROMOSOME-SPECIFIC CpG DRIVER
    # ============================================================

    pairs_chr <- pairs[chr == CHR]

    cpg_chr <- unique(
        pairs_chr[, .(
            cpg,
            chr_hg19,
            pos_hg19
        )]
    )

    setorder(cpg_chr, pos_hg19)

    cat("chr", CHR, " peak-CpG pairs:", nrow(pairs_chr), "\n")
    cat("chr", CHR, " unique CpGs:", nrow(cpg_chr), "\n")

    # ============================================================
    # DIVIDE CpGs AMONG ARRAY TASKS
    #
    # Task 1 gets rows 1,9,17,...
    # Task 2 gets rows 2,10,18,...
    # etc. when n_tasks=8.
    # Every CpG belongs to exactly one task.
    # ============================================================

    idx <- if (task_id <= nrow(cpg_chr)) {
        seq(from = task_id, to = nrow(cpg_chr), by = n_tasks)
    } else integer()

    cpg_task <- cpg_chr[idx]

    cat(
        "CpGs assigned to this task:",
        nrow(cpg_task),
        "\n"
    )

    # ============================================================
    # RESTART SUPPORT FOR THIS TASK
    # ============================================================

    completed <- character()

    if (file.exists(completed_file)) {
        completed <- unique(
            scan(
                completed_file,
                what = character(),
                quiet = TRUE
            )
        )
    }

    cat(
        "Previously completed in this task:",
        length(completed),
        "\n"
    )

    # ============================================================
    # MAIN LOOP FOR THIS ARRAY TASK
    # ============================================================

    for (k in seq_len(nrow(cpg_task))) {

        cpg_i <- cpg_task$cpg[k]
        pos_i <- as.integer(cpg_task$pos_hg19[k])

        cat(
            "\n--------------------------------------------\n",
            "Task ", task_id,
            " | CpG ", k, " / ", nrow(cpg_task),
            " | ", cpg_i,
            "\n",
            sep = ""
        )

        # --------------------------------------------------------
        # Skip successfully completed CpGs on restart
        # --------------------------------------------------------

        if (cpg_i %in% completed) {
            cat("Already completed -- skip\n")
            next
        }

        # --------------------------------------------------------
        # Extract observed WGBS mQTL rows
        # --------------------------------------------------------

        m <- mqtl[cpg == cpg_i]

        if (nrow(m) == 0L) {

            cat(
                cpg_i, "\tno_WGBS_rows\n",
                file = failed_file,
                append = TRUE,
                sep = ""
            )

            next
        }

        # --------------------------------------------------------
        # Convert WGBS p-value + direction to signed normal Z
        # --------------------------------------------------------

        m[, Z := sign(statistic) *
              qnorm(pvalue / 2, lower.tail = FALSE)]

        # --------------------------------------------------------
        # Build SSimp input
        #
        # A1 = WGBS other allele
        # A2 = WGBS counted/effect allele
        # --------------------------------------------------------

        ssimp_in <- m[, .(
            SNP = snpRsNum,
            A1  = snpAlt,
            A2  = snpCounted,
            Z,
            N = 165L
        )]

        ssimp_in <- ssimp_in[
            !is.na(SNP) &
            SNP != "" &
            !is.na(A1) &
            !is.na(A2) &
            is.finite(Z)
        ]

        ssimp_in <- ssimp_in[!duplicated(SNP)]

        if (nrow(ssimp_in) == 0L) {

            cat(
                cpg_i, "\tzero_usable_input_SNPs\n",
                file = failed_file,
                append = TRUE,
                sep = ""
            )

            next
        }

        # --------------------------------------------------------
        # CpG +/-100 kb using hg19 coordinates
        # --------------------------------------------------------

        start_i <- max(1L, pos_i - 100000L)
        end_i   <- pos_i + 100000L

        range_i <- paste0(
            CHR, ":",
            start_i,
            "-", CHR, ":",
            end_i
        )

        # --------------------------------------------------------
        # Unique temporary filenames
        # --------------------------------------------------------

        safe_cpg <- gsub(
            "[^A-Za-z0-9_-]",
            "_",
            cpg_i
        )

        input_i <- file.path(
            tmp_dir,
            paste0(safe_cpg, "_input.txt")
        )

        output_i <- file.path(
            tmp_dir,
            paste0(safe_cpg, "_output.txt")
        )

        log_i <- file.path(
            tmp_dir,
            paste0(safe_cpg, "_ssimp.log")
        )

        fwrite(
            ssimp_in,
            input_i,
            sep = "\t",
            quote = FALSE,
            na = "NA"
        )

        # --------------------------------------------------------
        # Run SSimp
        # --------------------------------------------------------

        ssimp_args <- c(
            "--gwas", shQuote(input_i),
            "--ref", shQuote(ref_vcf),
            "--out", shQuote(output_i),
            "--sample.names", shQuote(sample_file),
            "--impute.range", range_i,
            "--missingness", "dep"
        )

        # Remove stale output from a previous failed attempt.
        if (file.exists(output_i)) unlink(output_i)

        status <- system2(
            ssimp_path,
            args = ssimp_args,
            stdout = log_i,
            stderr = log_i
        )

        # --------------------------------------------------------
        # Check SSimp success
        # --------------------------------------------------------

        if (status != 0L || !file.exists(output_i)) {

            cat(
                cpg_i, "\tSSimp_failed\n",
                file = failed_file,
                append = TRUE,
                sep = ""
            )

            cat("SSimp FAILED -- keeping temp files/log\n")
            next
        }

        # --------------------------------------------------------
        # Read SSimp output
        # --------------------------------------------------------

        imp <- fread(output_i)

        if (nrow(imp) == 0L) {

            cat(
                cpg_i, "\tempty_SSimp_output\n",
                file = failed_file,
                append = TRUE,
                sep = ""
            )

            next
        }

        # --------------------------------------------------------
        # Mark observed vs imputed
        #
        # We validated:
        # source == GWAS  -> observed WGBS SNP
        # source == SSIMP -> imputed SNP
        # --------------------------------------------------------

        imp[, snp_status := fifelse(
            source == "GWAS",
            "observed",
            "imputed"
        )]

        # --------------------------------------------------------
        # QUALITY FILTER
        #
        # Keep all observed SNPs returned by SSimp.
        # Keep imputed SNPs only if r2.pred > 0.30.
        # --------------------------------------------------------

        imp <- imp[
            source == "GWAS" |
            (
                source == "SSIMP" &
                !is.na(r2.pred) &
                r2.pred > 0.30
            )
        ]

        # --------------------------------------------------------
        # Add CpG identity
        # --------------------------------------------------------

        imp[, cpg := cpg_i]
        imp[, cpg_pos_hg19 := pos_i]

        # --------------------------------------------------------
        # Put important fields first
        # --------------------------------------------------------

        first_cols <- c(
            "cpg",
            "cpg_pos_hg19",
            "SNP",
            "snp_status",
            "source",
            "r2.pred",
            "P.imp",
            "z_imp",
            "maf",
            "N.imp"
        )

        first_cols <- first_cols[
            first_cols %in% names(imp)
        ]

        setcolorder(
            imp,
            c(
                first_cols,
                setdiff(names(imp), first_cols)
            )
        )

        # --------------------------------------------------------
        # Append once to THIS TASK'S part file.
        #
        # Different Slurm tasks never write to the same file.
        # --------------------------------------------------------

        first_write <- !file.exists(combined_file)

        fwrite(
            imp,
            combined_file,
            sep = "\t",
            append = !first_write,
            col.names = first_write,
            quote = FALSE,
            na = "NA"
        )

        # --------------------------------------------------------
        # Mark successful completion
        # --------------------------------------------------------

        cat(
            cpg_i,
            "\n",
            file = completed_file,
            append = TRUE
        )

        cat(
            "SUCCESS",
            " | input=", nrow(ssimp_in),
            " | observed=", sum(imp$source == "GWAS"),
            " | imputed_r2>0.3=", sum(imp$source == "SSIMP"),
            " | retained=", nrow(imp),
            "\n",
            sep = ""
        )

        # --------------------------------------------------------
        # Delete temporary files on success
        # --------------------------------------------------------

        unlink(input_i)
        unlink(output_i)
        unlink(log_i)
    }

    cat(
        "\nTask ", task_id,
        " finished.\n",
        sep = ""
    )
}

# Process the requested chromosomes sequentially within this array task.
for (CHR in chromosomes) {
    run_chromosome(CHR)
}
