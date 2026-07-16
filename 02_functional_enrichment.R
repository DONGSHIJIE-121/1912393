suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(fgsea)
})

set.seed(1234)

# 1. Load data
root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
bulk_dir <- file.path(root_dir, "results", "01_bulk_preprocessing_DESeq2")
gene_set_dir <- file.path(root_dir, "data", "input", "gene_sets")
out_dir <- file.path(root_dir, "results", "02_functional_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

deprg_file <- file.path(bulk_dir, "DESeq2_PANoptosis_overlap_genes.csv")
all_result_file <- file.path(bulk_dir, "DESeq2_all_results.csv")
vst_file <- file.path(bulk_dir, "VST_after_ComBat_seq.csv")
meta_file <- file.path(bulk_dir, "training_metadata_used.csv")
hallmark_gmt <- file.path(gene_set_dir, "h.all.v2024.1.Hs.symbols.gmt")
kegg_gmt <- file.path(gene_set_dir, "c2.cp.kegg.v2024.1.Hs.symbols.gmt")
reactome_gmt <- file.path(gene_set_dir, "c2.cp.reactome.v2024.1.Hs.symbols.gmt")

required_files <- c(deprg_file, all_result_file, vst_file, meta_file, hallmark_gmt, kegg_gmt, reactome_gmt)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) stop("Missing files: ", paste(missing_files, collapse = ", "))

deprg <- fread(deprg_file, data.table = FALSE)
all_result <- fread(all_result_file, data.table = FALSE)
vst <- fread(vst_file, data.table = FALSE, check.names = FALSE)
meta <- fread(meta_file, data.table = FALSE, check.names = FALSE)

genes <- unique(na.omit(trimws(as.character(deprg$Gene))))
universe_symbols <- unique(na.omit(trimws(as.character(all_result$Gene))))

# 2. Convert gene identifiers
gene_map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe_map <- bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_entrez <- unique(gene_map$ENTREZID)
universe_entrez <- unique(universe_map$ENTREZID)

# 3. GO and KEGG over-representation analysis
go_res <- enrichGO(
  gene = gene_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

kegg_res <- enrichKEGG(
  gene = gene_entrez,
  universe = universe_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1
)

go_df <- as.data.frame(go_res)
kegg_df <- as.data.frame(kegg_res)
go_sig <- go_df[!is.na(go_df$p.adjust) & go_df$p.adjust < 0.05, , drop = FALSE]
kegg_sig <- kegg_df[!is.na(kegg_df$p.adjust) & kegg_df$p.adjust < 0.05, , drop = FALSE]

# 4. Build ranked gene list
vst_mat <- as.matrix(vst[, -1, drop = FALSE])
rownames(vst_mat) <- vst[[1]]
storage.mode(vst_mat) <- "numeric"
meta$Sample <- trimws(as.character(meta$Sample))
meta <- meta[match(colnames(vst_mat), meta$Sample), , drop = FALSE]
if (anyNA(meta$Sample)) stop("VST samples are missing from metadata")
meta$Group <- factor(tolower(as.character(meta$Group)), levels = c("healthy", "sarcopenia"))
meta$Dataset <- factor(as.character(meta$Dataset))
design <- model.matrix(~ Dataset + Group, data = meta)
coef_name <- grep("^Group", colnames(design), value = TRUE)
if (length(coef_name) != 1) stop("Unable to identify the sarcopenia coefficient")
fit <- eBayes(lmFit(vst_mat, design))
rank_table <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
rank_table$Gene <- rownames(rank_table)
rank_table <- rank_table[is.finite(rank_table$t) & !is.na(rank_table$Gene), , drop = FALSE]
rank_table <- rank_table[order(rank_table$t, decreasing = TRUE), , drop = FALSE]
ranks <- rank_table$t
names(ranks) <- rank_table$Gene
ranks <- sort(ranks, decreasing = TRUE)

# 5. Gene set enrichment analysis
run_fgsea <- function(gmt_file, database_name) {
  pathways <- gmtPathways(gmt_file)
  overlap <- vapply(pathways, function(x) sum(x %in% names(ranks)), numeric(1))
  pathways <- pathways[overlap >= 10 & overlap <= 500]
  set.seed(1234)
  out <- fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0
  )
  out <- as.data.frame(out)
  out$Database <- database_name
  out$Direction <- ifelse(out$NES > 0, "Sarcopenia", "Healthy")
  out$leadingEdge <- vapply(out$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  out[order(out$padj, -abs(out$NES)), , drop = FALSE]
}

hallmark_res <- run_fgsea(hallmark_gmt, "Hallmark")
kegg_gsea <- run_fgsea(kegg_gmt, "KEGG")
reactome_res <- run_fgsea(reactome_gmt, "Reactome")
all_gsea <- bind_rows(hallmark_res, kegg_gsea, reactome_res)
all_gsea_sig <- all_gsea[!is.na(all_gsea$padj) & all_gsea$padj < 0.05, , drop = FALSE]

# 6. Save output
fwrite(gene_map, file.path(out_dir, "DEPRG_SYMBOL_to_ENTREZID.csv"))
fwrite(go_df, file.path(out_dir, "ORA_GO_all.csv"))
fwrite(go_sig, file.path(out_dir, "ORA_GO_padj_lt_0.05.csv"))
fwrite(kegg_df, file.path(out_dir, "ORA_KEGG_all.csv"))
fwrite(kegg_sig, file.path(out_dir, "ORA_KEGG_padj_lt_0.05.csv"))
fwrite(rank_table, file.path(out_dir, "GSEA_rank_limma_dataset_group.csv"))
fwrite(hallmark_res, file.path(out_dir, "GSEA_Hallmark_all.csv"))
fwrite(kegg_gsea, file.path(out_dir, "GSEA_KEGG_all.csv"))
fwrite(reactome_res, file.path(out_dir, "GSEA_Reactome_all.csv"))
fwrite(all_gsea_sig, file.path(out_dir, "GSEA_all_databases_padj_lt_0.05.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
