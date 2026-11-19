# GTEx skeletal muscle H3K27ac haQTL and mQTL colocalization
# Reconstructed from original R history and manuscript Methods.
# Run from the directory containing the input files.
library(data.table)
library(coloc)

# 1. Parameters -------------------------------------------------------------
N_ha <- 108L
N_m <- 42L
p1 <- 1e-4; p2 <- 1e-4; p12 <- 1e-5
# Original P-value input and reconstructed variance input are run separately.
modes <- c("original", "estimated_varbeta")
# Optional recovered pair list, supplied alongside this script.
pair_file <- "GTEx_muscle_recovered_pairs.tsv"

# 2. Read and prepare summary statistics -----------------------------------
haqtl <- fread("muscle_haqtl_allch.txt")
maf <- fread("freq.txt")
mqtl <- fread("muscle_mqtl_for_coloc.txt")
stopifnot(!anyDuplicated(maf$snp))
stopifnot(all(c("snp", "gARE", "pvalues", "beta") %in% names(haqtl)),
          all(c("snp", "MAF") %in% names(maf)),
          all(c("snp", "cpg", "MAF", "pvalues", "beta", "SE") %in% names(mqtl)))
# Match GTEx GRCh38 SNP identifiers, as in the recovered R history.
haqtl2 <- merge(haqtl, maf, by = "snp")
mqtl2 <- merge(mqtl, maf[, .(snp)], by = "snp")
mqtl2[, varbeta := SE^2]
# Approximate variance assuming two-sided normal/Wald P-values.
haqtl2[, varbeta := NA_real_]
haqtl2[is.finite(beta) & beta != 0 & is.finite(pvalues) &
         pvalues > 0 & pvalues < 1,
       varbeta := (beta / qnorm(pvalues / 2, lower.tail = FALSE))^2]
haqtl3 <- haqtl2[, .(snp, gARE, pvalues, MAF, beta, varbeta)]
mqtl3 <- mqtl2[, .(snp, cpg, MAF, pvalues, beta, varbeta)]
stopifnot(!anyDuplicated(haqtl3[, .(gARE, snp)]),
          !anyDuplicated(mqtl3[, .(cpg, snp)]))

# 3. Select candidate pairs according to manuscript Methods ----------------
# High-confidence haQTL: lead SNP P < 1e-6.
# Link CpGs associated with that lead SNP at P < 1e-4.
# Retain tied lead SNPs; original tie-breaking rule is not recovered.
# Input mQTL file is assumed to contain the previously identified gCpGs.
leads <- haqtl[is.finite(pvalues), .SD[pvalues == min(pvalues)], by = gARE]
leads <- unique(leads[pvalues < 1e-6, .(gARE, snp)])
all_gHMs <- unique(leads$gARE)
linked <- merge(leads, unique(mqtl3[is.finite(pvalues) & pvalues < 1e-4,
                                    .(snp, cpg)]), by = "snp", allow.cartesian = TRUE)
pairs <- unique(linked[, .(gARE, cpg)])
fwrite(pairs, "GTEx_muscle_methods_pairs.tsv", sep = "\t")
# Compare reconstructed pairing against the actual recovered pair list.
if (file.exists(pair_file)) {
  recovered <- unique(fread(pair_file)[, .(gARE, cpg)])
  fwrite(fsetdiff(pairs, recovered), "GTEx_pairs_methods_only.tsv", sep = "\t")
  fwrite(fsetdiff(recovered, pairs), "GTEx_pairs_recovered_only.tsv", sep = "\t")
  message("Methods pairs: ", nrow(pairs), "; recovered pairs: ", nrow(recovered))
  # Reproduce the historical set by default when available.
  pairs <- recovered
}
message("High-confidence gHMs: ", length(all_gHMs), "; pairs: ", nrow(pairs))

# 4. Colocalization --------------------------------------------------------
setkey(haqtl3, gARE)
setkey(mqtl3, cpg)
valid_maf <- function(x) is.finite(x) & x > 0 & x <= 0.5
valid_var <- function(x) is.finite(x) & x > 0
run_pair <- function(g, c, mode) {
  h <- haqtl3[.(g), nomatch = 0L]
  m <- mqtl3[.(c), nomatch = 0L]
  # Preserve feature-specific SNP sets for coloc's sdY estimation, as in
  # original history. coloc.abf internally intersects SNP identifiers.
  h <- h[valid_maf(MAF)]
  m <- m[valid_maf(MAF) & is.finite(beta) & valid_var(varbeta)]
  if (mode == "original") {
    h <- h[is.finite(pvalues) & pvalues > 0 & pvalues <= 1]
    d1 <- list(snp=h$snp, pvalues=h$pvalues, MAF=h$MAF,
               N=N_ha, type="quant")
  } else {
    h <- h[is.finite(beta) & valid_var(varbeta)]
    d1 <- list(snp=h$snp, beta=h$beta, varbeta=h$varbeta,
               MAF=h$MAF, N=N_ha, type="quant")
  }
  d2 <- list(snp=m$snp, beta=m$beta, varbeta=m$varbeta,
             MAF=m$MAF, N=N_m, type="quant")
  n <- length(intersect(h$snp, m$snp))
  out <- data.table(gARE=g, cpg=c, nsnps=n, PP.H0.abf=NA_real_,
                    PP.H1.abf=NA_real_, PP.H2.abf=NA_real_,
                    PP.H3.abf=NA_real_, PP.H4.abf=NA_real_, status="ok")
  if (n == 0L) { out$status <- "no common SNPs"; return(out) }
  fit <- tryCatch(coloc.abf(d1, d2, p1=p1, p2=p2, p12=p12), error=identity)
  if (inherits(fit, "error")) { out$status <- conditionMessage(fit); return(out) }
  for (nm in paste0("PP.H", 0:4, ".abf")) out[[nm]] <- unname(fit$summary[nm])
  out$nsnps <- as.integer(fit$summary["nsnps"])
  out
}

# 5. Export pair-level and gHM-level results --------------------------------
results <- list()
for (mode in modes) {
  res <- rbindlist(lapply(seq_len(nrow(pairs)), function(i)
    run_pair(pairs$gARE[i], pairs$cpg[i], mode)))
  fwrite(res, paste0("GTEx_muscle_coloc_", mode, ".tsv"), sep="\t")
  results[[mode]] <- res
  best <- res[status == "ok" & is.finite(PP.H4.abf),
              .(PP4=max(PP.H4.abf)), by=gARE]
  ghm <- merge(data.table(gARE=unique(c(all_gHMs, pairs$gARE))), best,
               by="gARE", all.x=TRUE)
  ghm[, category := fcase(is.na(PP4), "No test", PP4 >= 0.8, "Strong",
                         PP4 >= 0.5, "Suggestive", default="No colocalization")]
  ghm[is.na(PP4) & gARE %in% pairs$gARE, category := "Test failed"]
  fwrite(ghm, paste0("GTEx_muscle_gHM_", mode, ".tsv"), sep="\t")
  print(ghm[, .N, by=category])
  print(res[status == "ok", .(pairs=.N, median_SNPs=median(nsnps),
                               mean_SNPs=mean(nsnps), IQR_SNPs=IQR(nsnps))])
  # Published targets: 169 strong, 54 suggestive, 37 no colocalization,
  # 225 no test; median 454 SNPs. Do not force these counts.
}

# 6. Compare reconstructed variances with original approach ----------------
if (all(c("original", "estimated_varbeta") %in% names(results))) {
  comparison <- merge(results$original, results$estimated_varbeta,
                      by=c("gARE", "cpg"), suffixes=c("_original", "_estimated"))
  comparison[, delta_PP4 := PP.H4.abf_estimated - PP.H4.abf_original]
  fwrite(comparison, "GTEx_muscle_coloc_method_comparison.tsv", sep="\t")
}
capture.output(sessionInfo(), file="GTEx_muscle_sessionInfo.txt")
