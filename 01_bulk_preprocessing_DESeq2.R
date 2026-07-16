suppressPackageStartupMessages({
  library(data.table)
  library(sva)
  library(DESeq2)
})

set.seed(1234)

# 1. Load data
root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(root_dir, "data", "input", "bulk")
gene_set_dir <- file.path(root_dir, "data", "input", "gene_sets")
out_dir <- file.path(root_dir, "results", "01_bulk_preprocessing_DESeq2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE111006", "GSE111010", "GSE226151", "GSE238215")
count_files <- setNames(file.path(input_dir, paste0(datasets, "_counts.csv")), datasets)
meta_file <- file.path(input_dir, "training_metadata.csv")
pan_file <- file.path(gene_set_dir, "PANoptosis_781.csv")

required_files <- c(count_files, meta_file, pan_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) stop("Missing files: ", paste(missing_files, collapse = ", "))

read_count_matrix <- function(path) {
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  genes <- trimws(as.character(x[[1]]))
  x[[1]] <- NULL
  mat <- as.matrix(x)
  rownames(mat) <- genes
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat) | mat < 0] <- 0
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  if (anyDuplicated(rownames(mat))) mat <- rowsum(mat, rownames(mat), reorder = FALSE)
  mat <- round(mat)
  storage.mode(mat) <- "integer"
  mat
}

counts_list <- lapply(count_files, read_count_matrix)
common_genes <- Reduce(intersect, lapply(counts_list, rownames))
counts_list <- lapply(counts_list, function(x) x[common_genes, , drop = FALSE])
merged_counts <- do.call(cbind, counts_list)
if (anyDuplicated(colnames(merged_counts))) stop("Duplicated sample identifiers in count matrices")

meta <- fread(meta_file, data.table = FALSE, check.names = FALSE)
required_meta <- c("Sample", "Group", "Dataset")
if (!all(required_meta %in% colnames(meta))) stop("training_metadata.csv must contain Sample, Group, and Dataset")
meta$Sample <- trimws(as.character(meta$Sample))
meta$Group <- factor(trimws(tolower(as.character(meta$Group))), levels = c("healthy", "sarcopenia"))
meta$Dataset <- factor(trimws(as.character(meta$Dataset)), levels = datasets)
if (anyNA(meta$Group) || anyNA(meta$Dataset)) stop("Invalid Group or Dataset values")
if (!all(meta$Sample %in% colnames(merged_counts))) stop("Metadata samples are missing from count matrices")
meta <- meta[match(colnames(merged_counts), meta$Sample), , drop = FALSE]
if (anyNA(meta$Sample)) stop("Count matrix samples are missing from metadata")
rownames(meta) <- meta$Sample

# 2. Filter low counts
keep <- rowSums(merged_counts >= 10) >= ceiling(ncol(merged_counts) * 0.50)
counts_filtered <- merged_counts[keep, , drop = FALSE]

# 3. ComBat-seq batch adjustment
counts_combat <- sva::ComBat_seq(
  counts = counts_filtered,
  batch = as.character(meta$Dataset),
  group = as.character(meta$Group)
)
counts_combat <- round(counts_combat)
storage.mode(counts_combat) <- "integer"

# 4. DESeq2 analysis
dds <- DESeqDataSetFromMatrix(
  countData = counts_combat,
  colData = meta,
  design = ~ Dataset + Group
)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("Group", "sarcopenia", "healthy"), alpha = 0.05)
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)
res_df <- res_df[, c("Gene", setdiff(colnames(res_df), "Gene"))]

# 5. Extract results
padj_cut <- 0.05
lfc_cut <- 0.25
deg <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < padj_cut &
    abs(res_df$log2FoldChange) > lfc_cut,
  , drop = FALSE
]
deg$Direction <- ifelse(deg$log2FoldChange > 0, "Up", "Down")
deg <- deg[order(deg$padj, -abs(deg$log2FoldChange)), , drop = FALSE]

pan <- fread(pan_file, data.table = FALSE, check.names = FALSE)
gene_col <- intersect(c("Gene", "gene", "SYMBOL", "Symbol", "symbol"), colnames(pan))[1]
if (is.na(gene_col)) gene_col <- colnames(pan)[1]
pan_genes <- unique(trimws(as.character(pan[[gene_col]])))
deprg <- deg[deg$Gene %in% pan_genes, , drop = FALSE]

# 6. Save output
dds_before <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = meta,
  design = ~ Dataset + Group
)
dds_before <- estimateSizeFactors(dds_before)
dds_before <- estimateDispersions(dds_before, quiet = TRUE)
vsd_before <- vst(dds_before, blind = FALSE)
vsd_after <- vst(dds, blind = FALSE)

write_matrix <- function(mat, path) {
  fwrite(data.table(Gene = rownames(mat), as.data.frame(mat, check.names = FALSE)), path)
}

pca_table <- function(mat, meta) {
  pca <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  out <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC1_variance_percent = variance[1],
    PC2_variance_percent = variance[2],
    stringsAsFactors = FALSE
  )
  cbind(out, meta[out$Sample, c("Group", "Dataset"), drop = FALSE])
}

fwrite(data.table(Gene = rownames(merged_counts), as.data.frame(merged_counts, check.names = FALSE)), file.path(out_dir, "merged_common_gene_counts.csv"))
fwrite(data.table(Gene = rownames(counts_filtered), as.data.frame(counts_filtered, check.names = FALSE)), file.path(out_dir, "filtered_counts_before_ComBat_seq.csv"))
fwrite(data.table(Gene = rownames(counts_combat), as.data.frame(counts_combat, check.names = FALSE)), file.path(out_dir, "ComBat_seq_adjusted_counts.csv"))
fwrite(meta, file.path(out_dir, "training_metadata_used.csv"), row.names = FALSE)
fwrite(res_df, file.path(out_dir, "DESeq2_all_results.csv"))
fwrite(deg, file.path(out_dir, "DESeq2_DEGs_padj0.05_absLFC0.25.csv"))
fwrite(deprg, file.path(out_dir, "DESeq2_PANoptosis_overlap_genes.csv"))
write_matrix(assay(vsd_before), file.path(out_dir, "VST_before_ComBat_seq.csv"))
write_matrix(assay(vsd_after), file.path(out_dir, "VST_after_ComBat_seq.csv"))
fwrite(pca_table(assay(vsd_before), meta), file.path(out_dir, "PCA_before_ComBat_seq.csv"))
fwrite(pca_table(assay(vsd_after), meta), file.path(out_dir, "PCA_after_ComBat_seq.csv"))
saveRDS(dds, file.path(out_dir, "DESeq2_dataset.rds"))
saveRDS(vsd_after, file.path(out_dir, "VST_object.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
