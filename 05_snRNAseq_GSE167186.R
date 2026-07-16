suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(Matrix)
})

set.seed()

# 1. Load data
root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(root_dir, "data", "input", "snrna")
out_dir <- file.path(root_dir, "results", "05_snRNAseq_GSE167186")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

object_file <- file.path(input_dir, "GSE167186_raw_seurat.rds")
doublet_file <- file.path(input_dir, "GSE167186_scrublet_doublet_calls.csv")
pan_file <- file.path(input_dir, "PANoptosis_top50.csv")

required_files <- c(object_file, pan_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) stop("Missing files: ", paste(missing_files, collapse = ", "))

obj <- readRDS(object_file)
DefaultAssay(obj) <- "RNA"

resolve_column <- function(meta, candidates) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) stop("Missing metadata column: ", paste(candidates, collapse = "/"))
  hit[1]
}

parse_logical <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(x))
  out[x %in% c("true", "t", "1", "yes", "y", "doublet")] <- TRUE
  out[x %in% c("false", "f", "0", "no", "n", "singlet")] <- FALSE
  out
}

sample_col <- resolve_column(obj@meta.data, c("Sample", "SampleID", "sample", "sample_id", "orig.ident"))
group_col <- resolve_column(obj@meta.data, c("Group", "Age_group", "age_group", "group"))
if (!"percent.mt" %in% colnames(obj@meta.data)) obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# 2. Quality control
meta <- obj@meta.data
meta$Cell <- rownames(meta)
meta$Sample_qc <- as.character(meta[[sample_col]])
meta$log_counts <- log1p(meta$nCount_RNA)
meta$log_features <- log1p(meta$nFeature_RNA)

qc_meta <- meta %>%
  group_by(Sample_qc) %>%
  mutate(
    counts_median = median(log_counts, na.rm = TRUE),
    counts_mad = mad(log_counts, na.rm = TRUE),
    features_median = median(log_features, na.rm = TRUE),
    features_mad = mad(log_features, na.rm = TRUE),
    mt_median = median(percent.mt, na.rm = TRUE),
    mt_mad = mad(percent.mt, na.rm = TRUE),
    keep_counts = abs(log_counts - counts_median) <= 3 * counts_mad,
    keep_features = abs(log_features - features_median) <= 3 * features_mad,
    keep_mt = percent.mt <= pmin(mt_median + 3 * mt_mad, 20)
  ) %>%
  ungroup()

if ("is_doublet" %in% colnames(qc_meta)) {
  qc_meta$doublet_flag <- parse_logical(qc_meta$is_doublet)
} else {
  if (!file.exists(doublet_file)) stop("Scrublet calls are required when is_doublet is absent from the Seurat metadata")
  doublet_calls <- fread(doublet_file, data.table = FALSE, check.names = FALSE)
  if (!all(c("Cell", "is_doublet") %in% colnames(doublet_calls))) stop("Scrublet file must contain Cell and is_doublet")
  qc_meta$doublet_flag <- parse_logical(doublet_calls$is_doublet[match(qc_meta$Cell, doublet_calls$Cell)])
  if (anyNA(qc_meta$doublet_flag)) stop("Scrublet calls are missing or invalid for some nuclei")
}

qc_meta$keep_final <- qc_meta$keep_counts & qc_meta$keep_features & qc_meta$keep_mt & !qc_meta$doublet_flag
obj <- subset(obj, cells = qc_meta$Cell[qc_meta$keep_final])

# 3. Normalize and cluster
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, seed.use = 1235, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 1235, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, seed.use = 1235, verbose = FALSE)

# 4. Identify cluster markers
markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# 5. Annotate cell types
cluster_to_celltype <- c(
  "0" = "Fast skeletal myonuclei",
  "2" = "Fast skeletal myonuclei",
  "9" = "Fast skeletal myonuclei",
  "1" = "Slow skeletal myonuclei",
  "3" = "Slow skeletal myonuclei",
  "5" = "NMJ-like myonuclei",
  "10" = "Satellite cells",
  "4" = "FAPs",
  "6" = "Endothelial cells",
  "8" = "Pericyte/smooth muscle cells",
  "7" = "Macrophages",
  "11" = "T/NK cells"
)
clusters_present <- sort(unique(as.character(obj$seurat_clusters)))
unmapped <- setdiff(clusters_present, names(cluster_to_celltype))
if (length(unmapped) > 0) stop("Unmapped clusters: ", paste(unmapped, collapse = ", "))
obj$Cell_type <- unname(cluster_to_celltype[as.character(obj$seurat_clusters)])
obj$Age_group <- tolower(trimws(as.character(obj@meta.data[[group_col]])))
obj$Age_group[obj$Age_group %in% c("aged", "older", "elderly")] <- "old"
obj$Age_group[obj$Age_group %in% c("younger", "adult")] <- "young"
obj$Age_group <- factor(obj$Age_group, levels = c("young", "old"))
obj$Sample <- as.character(obj@meta.data[[sample_col]])
Idents(obj) <- "Cell_type"

# 6. Calculate PANoptosis-associated module score
pan <- fread(pan_file, data.table = FALSE, check.names = FALSE)
gene_col <- intersect(c("Gene", "gene", "SYMBOL", "Symbol", "symbol"), colnames(pan))[1]
if (is.na(gene_col)) gene_col <- colnames(pan)[1]
pan_genes <- unique(trimws(as.character(pan[[gene_col]])))
pan_genes <- intersect(pan_genes, rownames(obj))
if (length(pan_genes) < 5) stop("Too few PANoptosis genes are present in the object")
obj <- AddModuleScore(obj, features = list(pan_genes), name = "PANoptosis_score", seed = 1235)
obj$PANoptosis_score <- obj$PANoptosis_score1
obj$PANoptosis_score1 <- NULL

# 7. Summarize composition and hub-gene expression
composition_sample <- as.data.frame(table(obj$Sample, obj$Age_group, obj$Cell_type), stringsAsFactors = FALSE)
colnames(composition_sample) <- c("Sample", "Age_group", "Cell_type", "Nuclei")
composition_sample <- composition_sample %>%
  group_by(Sample) %>%
  mutate(Fraction = Nuclei / sum(Nuclei)) %>%
  ungroup()

umap <- as.data.frame(Embeddings(obj, "umap"))
umap$Cell <- rownames(umap)
umap <- cbind(umap, obj@meta.data[umap$Cell, c("Sample", "Age_group", "seurat_clusters", "Cell_type", "PANoptosis_score"), drop = FALSE])

score_summary <- obj@meta.data %>%
  group_by(Cell_type, Age_group) %>%
  summarise(
    Nuclei = n(),
    Mean_score = mean(PANoptosis_score, na.rm = TRUE),
    Median_score = median(PANoptosis_score, na.rm = TRUE),
    .groups = "drop"
  )

get_normalized_data <- function(object) {
  tryCatch(
    GetAssayData(object, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(object, assay = "RNA", slot = "data")
  )
}

hub_genes <- intersect(c("NTRK1", "TPT1", "TRAP1"), rownames(obj))
expr <- get_normalized_data(obj)[hub_genes, , drop = FALSE]
hub_long <- rbindlist(lapply(hub_genes, function(gene) {
  data.frame(
    Cell = colnames(expr),
    Gene = gene,
    Expression = as.numeric(expr[gene, ]),
    stringsAsFactors = FALSE
  )
}))
hub_long <- merge(
  hub_long,
  data.frame(Cell = rownames(obj@meta.data), obj@meta.data[, c("Sample", "Age_group", "Cell_type"), drop = FALSE]),
  by = "Cell",
  all.x = TRUE
)

hub_celltype_summary <- hub_long %>%
  group_by(Gene, Cell_type, Age_group) %>%
  summarise(
    Nuclei = n(),
    Mean_expression = mean(Expression),
    Median_expression = median(Expression),
    Percent_expressing = 100 * mean(Expression > 0),
    .groups = "drop"
  )

hub_sample_summary <- hub_long %>%
  group_by(Gene, Sample, Age_group, Cell_type) %>%
  summarise(
    Nuclei = n(),
    Mean_expression = mean(Expression),
    Percent_expressing = 100 * mean(Expression > 0),
    .groups = "drop"
  )

hub_age_stats <- hub_sample_summary %>%
  group_by(Gene, Cell_type) %>%
  summarise(
    P_mean_expression = if (n_distinct(Age_group) == 2) wilcox.test(Mean_expression ~ Age_group, exact = FALSE)$p.value else NA_real_,
    P_percent_expressing = if (n_distinct(Age_group) == 2) wilcox.test(Percent_expressing ~ Age_group, exact = FALSE)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  group_by(Gene) %>%
  mutate(
    FDR_mean_expression = p.adjust(P_mean_expression, method = "BH"),
    FDR_percent_expressing = p.adjust(P_percent_expressing, method = "BH")
  ) %>%
  ungroup()

# 8. Save output
fwrite(qc_meta, file.path(out_dir, "QC_nucleus_level.csv"))
fwrite(markers, file.path(out_dir, "Seurat_cluster_markers.csv"))
fwrite(data.frame(Cluster = names(cluster_to_celltype), Cell_type = unname(cluster_to_celltype)), file.path(out_dir, "Cluster_annotation_map.csv"))
fwrite(composition_sample, file.path(out_dir, "Cell_type_composition_by_sample.csv"))
fwrite(umap, file.path(out_dir, "UMAP_coordinates_and_metadata.csv"))
fwrite(score_summary, file.path(out_dir, "PANoptosis_module_score_summary.csv"))
fwrite(hub_celltype_summary, file.path(out_dir, "Hub_gene_expression_by_cell_type_and_age.csv"))
fwrite(hub_sample_summary, file.path(out_dir, "Hub_gene_sample_level_summary.csv"))
fwrite(hub_age_stats, file.path(out_dir, "Hub_gene_age_group_statistics.csv"))
saveRDS(obj, file.path(out_dir, "GSE167186_processed_annotated.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
