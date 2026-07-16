suppressPackageStartupMessages({
  library(data.table)
})

set.seed()

# 1. Load data
root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(root_dir, "data", "input", "immune")
out_dir <- file.path(root_dir, "results", "04_CIBERSORT_immune_infiltration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tpm_file <- file.path(input_dir, "training_TPM_gene_by_sample.csv")
meta_file <- file.path(input_dir, "training_metadata.csv")
lm22_file <- file.path(input_dir, "LM22.txt")
cibersort_file <- file.path(input_dir, "CIBERSORT.R")

required_files <- c(tpm_file, meta_file, lm22_file, cibersort_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) stop("Missing files: ", paste(missing_files, collapse = ", "))

read_expression <- function(path) {
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  genes <- trimws(as.character(x[[1]]))
  x[[1]] <- NULL
  mat <- as.matrix(x)
  rownames(mat) <- genes
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat) | mat < 0] <- 0
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  if (anyDuplicated(rownames(mat))) {
    dt <- as.data.table(mat, keep.rownames = "Gene")
    dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = Gene]
    mat <- as.matrix(dt[, setdiff(colnames(dt), "Gene"), with = FALSE])
    rownames(mat) <- dt$Gene
  }
  mat
}

expr <- read_expression(tpm_file)
meta <- fread(meta_file, data.table = FALSE, check.names = FALSE)
if (!all(c("Sample", "Group") %in% colnames(meta))) stop("Metadata must contain Sample and Group")
meta$Sample <- trimws(as.character(meta$Sample))
meta$Group <- factor(tolower(trimws(as.character(meta$Group))), levels = c("healthy", "sarcopenia"))
meta <- meta[match(colnames(expr), meta$Sample), , drop = FALSE]
if (anyNA(meta$Sample) || anyNA(meta$Group)) stop("Expression samples and metadata do not match")

# 2. Prepare CIBERSORT input
mixture_file <- file.path(out_dir, "CIBERSORT_mixture_TPM.txt")
write.table(
  cbind(GeneSymbol = rownames(expr), expr),
  file = mixture_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 3. Run CIBERSORT
source(cibersort_file)
if (!exists("CIBERSORT")) stop("CIBERSORT function was not found after sourcing CIBERSORT.R")
cibersort_result <- CIBERSORT(
  sig_matrix = lm22_file,
  mixture_file = mixture_file,
  perm = 1000,
  QN = FALSE
)
result <- as.data.frame(cibersort_result)
result$Sample <- rownames(result)

# 4. Extract cell fractions
metric_columns <- intersect(c("P-value", "Correlation", "RMSE", "Sample"), colnames(result))
cell_columns <- setdiff(colnames(result), metric_columns)
fractions <- result[, cell_columns, drop = FALSE]
fractions[] <- lapply(fractions, as.numeric)
rownames(fractions) <- result$Sample
meta <- meta[match(rownames(fractions), meta$Sample), , drop = FALSE]
if (anyNA(meta$Sample)) stop("CIBERSORT samples and metadata do not match")

# 5. Group comparisons
group_stats <- rbindlist(lapply(cell_columns, function(cell) {
  healthy <- fractions[meta$Group == "healthy", cell]
  sarcopenia <- fractions[meta$Group == "sarcopenia", cell]
  test <- wilcox.test(healthy, sarcopenia, exact = FALSE)
  data.frame(
    Cell_type = cell,
    Median_healthy = median(healthy, na.rm = TRUE),
    Median_sarcopenia = median(sarcopenia, na.rm = TRUE),
    P_value = test$p.value,
    stringsAsFactors = FALSE
  )
}))
group_stats$FDR <- p.adjust(group_stats$P_value, method = "BH")

# 6. Hub-gene correlations
safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) return(c(rho = NA_real_, p = NA_real_))
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(test$estimate), p = test$p.value)
}

hub_genes <- c("NTRK1", "TPT1", "TRAP1")
missing_hub <- setdiff(hub_genes, rownames(expr))
if (length(missing_hub) > 0) stop("Hub genes missing from TPM matrix: ", paste(missing_hub, collapse = ", "))
expr_log <- log2(expr + 1)
correlation_stats <- rbindlist(lapply(hub_genes, function(gene) {
  rbindlist(lapply(cell_columns, function(cell) {
    test <- safe_spearman(
      as.numeric(expr_log[gene, rownames(fractions)]),
      as.numeric(fractions[[cell]])
    )
    data.frame(
      Gene = gene,
      Cell_type = cell,
      Rho = as.numeric(test["rho"]),
      P_value = as.numeric(test["p"]),
      stringsAsFactors = FALSE
    )
  }))
}))
correlation_stats$FDR <- ave(correlation_stats$P_value, correlation_stats$Gene, FUN = function(x) p.adjust(x, method = "BH"))

# 7. Save output
fwrite(result, file.path(out_dir, "CIBERSORT_LM22_all_samples.csv"))
if ("P-value" %in% colnames(result)) {
  fwrite(result[result$`P-value` < 0.05, , drop = FALSE], file.path(out_dir, "CIBERSORT_LM22_P_lt_0.05.csv"))
}
fwrite(data.frame(Sample = rownames(fractions), fractions, check.names = FALSE), file.path(out_dir, "CIBERSORT_cell_fractions.csv"))
fwrite(group_stats, file.path(out_dir, "CIBERSORT_group_comparison_Wilcoxon.csv"))
fwrite(correlation_stats, file.path(out_dir, "CIBERSORT_hub_gene_Spearman_correlations.csv"))
fwrite(meta, file.path(out_dir, "CIBERSORT_metadata_used.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
