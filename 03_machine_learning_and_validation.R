suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(randomForest)
  library(xgboost)
  library(pROC)
  library(rms)
})

set.seed(1234)

# 1. Load data
root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
bulk_dir <- file.path(root_dir, "results", "01_bulk_preprocessing_DESeq2")
input_dir <- file.path(root_dir, "data", "input", "bulk")
out_dir <- file.path(root_dir, "results", "03_machine_learning_and_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

train_expr_file <- file.path(bulk_dir, "VST_after_ComBat_seq.csv")
train_meta_file <- file.path(bulk_dir, "training_metadata_used.csv")
feature_file <- file.path(bulk_dir, "DESeq2_PANoptosis_overlap_genes.csv")
valid_expr_file <- file.path(input_dir, "GSE111016_VST.csv")
valid_meta_file <- file.path(input_dir, "GSE111016_metadata.csv")

required_files <- c(train_expr_file, train_meta_file, feature_file, valid_expr_file, valid_meta_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) stop("Missing files: ", paste(missing_files, collapse = ", "))

read_expression <- function(path) {
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  genes <- trimws(as.character(x[[1]]))
  x[[1]] <- NULL
  mat <- as.matrix(x)
  rownames(mat) <- genes
  storage.mode(mat) <- "numeric"
  if (anyDuplicated(rownames(mat))) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }
  mat
}

align_metadata <- function(meta, samples) {
  if (!all(c("Sample", "Group") %in% colnames(meta))) stop("Metadata must contain Sample and Group")
  meta$Sample <- trimws(as.character(meta$Sample))
  meta$Group <- factor(tolower(trimws(as.character(meta$Group))), levels = c("healthy", "sarcopenia"))
  meta <- meta[match(samples, meta$Sample), , drop = FALSE]
  if (anyNA(meta$Sample) || anyNA(meta$Group)) stop("Expression samples and metadata do not match")
  rownames(meta) <- meta$Sample
  meta
}

make_stratified_foldid <- function(y, k, seed) {
  set.seed(seed)
  foldid <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- which(y == cls)
    foldid[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  foldid
}

roc_summary <- function(expr, group, genes, cohort) {
  out <- lapply(genes, function(gene) {
    direction <- if (
      median(expr[gene, group == "healthy"], na.rm = TRUE) >
      median(expr[gene, group == "sarcopenia"], na.rm = TRUE)
    ) ">" else "<"
    roc_obj <- roc(
      response = factor(group, levels = c("healthy", "sarcopenia")),
      predictor = as.numeric(expr[gene, ]),
      levels = c("healthy", "sarcopenia"),
      direction = direction,
      ci = TRUE,
      quiet = TRUE
    )
    ci_obj <- ci.auc(roc_obj, method = "delong")
    data.frame(
      Cohort = cohort,
      Gene = gene,
      Direction = direction,
      AUC = as.numeric(auc(roc_obj)),
      CI_lower = as.numeric(ci_obj[1]),
      CI_upper = as.numeric(ci_obj[3]),
      stringsAsFactors = FALSE
    )
  })
  rbindlist(out)
}

train_expr <- read_expression(train_expr_file)
train_meta <- align_metadata(fread(train_meta_file, data.table = FALSE), colnames(train_expr))
features_df <- fread(feature_file, data.table = FALSE)
if (!"Gene" %in% colnames(features_df)) stop("Feature file must contain a Gene column")
features <- unique(trimws(as.character(features_df$Gene)))
features <- intersect(features, rownames(train_expr))
if (length(features) < 2) stop("Fewer than two candidate genes are present in the training matrix")

x <- t(train_expr[features, , drop = FALSE])
storage.mode(x) <- "numeric"
y <- ifelse(train_meta$Group == "sarcopenia", 1, 0)

# 2. LASSO
lasso_foldid <- make_stratified_foldid(y, 10, 124)
set.seed(124)
lasso_cv <- cv.glmnet(
  x = x,
  y = y,
  family = "binomial",
  alpha = 1,
  foldid = lasso_foldid,
  standardize = TRUE,
  type.measure = "deviance"
)
lasso_coef <- as.matrix(coef(lasso_cv, s = "lambda.1se"))
lasso_table <- data.frame(
  Gene = rownames(lasso_coef),
  Coefficient = as.numeric(lasso_coef),
  stringsAsFactors = FALSE
)
lasso_table <- lasso_table[
  lasso_table$Gene != "(Intercept)" & lasso_table$Coefficient != 0,
  , drop = FALSE
]
lasso_genes <- lasso_table$Gene

# 3. Random forest
set.seed(123)
rf_fit <- randomForest(
  x = as.data.frame(x, check.names = FALSE),
  y = factor(train_meta$Group, levels = c("healthy", "sarcopenia")),
  ntree = 300,
  mtry = max(1, floor(sqrt(ncol(x)))),
  importance = TRUE
)
rf_error <- data.frame(Tree = seq_len(nrow(rf_fit$err.rate)), rf_fit$err.rate, check.names = FALSE)
rf_importance <- as.data.frame(importance(rf_fit))
rf_importance$Gene <- rownames(rf_importance)
rf_importance <- rf_importance[order(rf_importance$MeanDecreaseGini, decreasing = TRUE), , drop = FALSE]
rf_top10 <- head(rf_importance, 10)

# 4. XGBoost
xgb_foldid <- make_stratified_foldid(y, 5, 126)
xgb_folds <- split(seq_along(y), xgb_foldid)
dtrain <- xgb.DMatrix(data = x, label = y)
xgb_params <- list(
  booster = "gbtree",
  objective = "binary:logistic",
  eval_metric = "auc",
  eta = 0.05,
  max_depth = 3,
  min_child_weight = 1,
  subsample = 0.8,
  colsample_bytree = 0.8,
  gamma = 0,
  lambda = 1
)
set.seed(126)
xgb_cv <- xgb.cv(
  params = xgb_params,
  data = dtrain,
  nrounds = 300,
  folds = xgb_folds,
  early_stopping_rounds = 20,
  maximize = TRUE,
  verbose = 0
)
best_nrounds <- xgb_cv$best_iteration
if (is.null(best_nrounds) || is.na(best_nrounds)) {
  best_nrounds <- which.max(xgb_cv$evaluation_log$test_auc_mean)
}
set.seed(126)
xgb_fit <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)
xgb_importance <- as.data.frame(xgb.importance(feature_names = colnames(x), model = xgb_fit))
xgb_top10 <- head(xgb_importance, 10)

# 5. Consensus candidates
consensus_genes <- Reduce(intersect, list(lasso_genes, rf_top10$Gene, xgb_top10$Feature))
expected_consensus <- c("NTRK1", "TPT1", "TRAP1")
if (!setequal(consensus_genes, expected_consensus)) {
  stop(
    "Consensus genes differ from the manuscript. Observed: ",
    paste(consensus_genes, collapse = ", "),
    "; expected: ",
    paste(expected_consensus, collapse = ", ")
  )
}
consensus_genes <- expected_consensus
consensus_table <- data.frame(Gene = consensus_genes, stringsAsFactors = FALSE)

# 6. rms logistic model and nomogram
model_data <- data.frame(
  Outcome = y,
  x[, consensus_genes, drop = FALSE],
  check.names = FALSE
)
model_formula <- as.formula(
  paste("Outcome ~", paste(sprintf("`%s`", consensus_genes), collapse = " + "))
)

dd <- datadist(model_data)
old_datadist <- getOption("datadist")
options(datadist = "dd")

lrm_fit <- lrm(
  formula = model_formula,
  data = model_data,
  x = TRUE,
  y = TRUE
)

nomogram_fit <- nomogram(
  lrm_fit,
  fun = plogis,
  fun.at = seq(0.1, 0.9, by = 0.1),
  lp = FALSE,
  funlabel = "Sarcopenia probability"
)

pdf(file.path(out_dir, "Three_gene_nomogram.pdf"), width = 9, height = 6)
plot(nomogram_fit)
dev.off()

train_linear_predictor <- as.numeric(predict(lrm_fit, newdata = model_data, type = "lp"))
train_probability <- plogis(train_linear_predictor)

# 7. Bootstrap bias-corrected calibration
set.seed(127)
calibration_fit <- calibrate(
  lrm_fit,
  method = "boot",
  B = 1000
)
calibration_table <- as.data.frame(unclass(calibration_fit), check.names = FALSE)

pdf(file.path(out_dir, "Three_gene_calibration_bootstrap_B1000.pdf"), width = 6, height = 6)
plot(
  calibration_fit,
  subtitles = FALSE,
  xlab = "Predicted probability",
  ylab = "Observed probability"
)
dev.off()

# 8. Training-cohort model discrimination and decision curve
model_roc <- roc(
  response = factor(train_meta$Group, levels = c("healthy", "sarcopenia")),
  predictor = train_probability,
  levels = c("healthy", "sarcopenia"),
  direction = "<",
  ci = TRUE,
  quiet = TRUE
)
model_roc_ci <- ci.auc(model_roc, method = "delong")
model_roc_summary <- data.frame(
  Cohort = "Training",
  Model = "NTRK1_TPT1_TRAP1",
  AUC = as.numeric(auc(model_roc)),
  CI_lower = as.numeric(model_roc_ci[1]),
  CI_upper = as.numeric(model_roc_ci[3])
)

thresholds <- seq(0.01, 0.99, by = 0.01)
decision_curve <- rbindlist(lapply(thresholds, function(threshold) {
  predicted_positive <- train_probability >= threshold
  tp <- sum(predicted_positive & y == 1)
  fp <- sum(predicted_positive & y == 0)
  n <- length(y)
  prevalence <- mean(y)
  data.frame(
    Threshold = threshold,
    Model = tp / n - fp / n * threshold / (1 - threshold),
    Treat_all = prevalence - (1 - prevalence) * threshold / (1 - threshold),
    Treat_none = 0
  )
}))

pdf(file.path(out_dir, "Three_gene_decision_curve_training.pdf"), width = 6, height = 6)
matplot(
  decision_curve$Threshold,
  decision_curve[, c("Model", "Treat_all", "Treat_none")],
  type = "l",
  lty = c(1, 2, 3),
  xlab = "Threshold probability",
  ylab = "Net benefit"
)
legend(
  "topright",
  legend = c("Three-gene model", "Treat all", "Treat none"),
  lty = c(1, 2, 3),
  bty = "n"
)
dev.off()

# 9. Single-gene training and held-out validation
train_roc <- roc_summary(
  train_expr,
  as.character(train_meta$Group),
  consensus_genes,
  "Training"
)

valid_expr <- read_expression(valid_expr_file)
valid_meta <- align_metadata(fread(valid_meta_file, data.table = FALSE), colnames(valid_expr))
missing_valid <- setdiff(consensus_genes, rownames(valid_expr))
if (length(missing_valid) > 0) {
  stop("Consensus genes missing from validation matrix: ", paste(missing_valid, collapse = ", "))
}
valid_roc <- roc_summary(
  valid_expr,
  as.character(valid_meta$Group),
  consensus_genes,
  "GSE111016"
)

valid_expression_stats <- rbindlist(lapply(consensus_genes, function(gene) {
  healthy <- as.numeric(valid_expr[gene, valid_meta$Group == "healthy"])
  sarcopenia <- as.numeric(valid_expr[gene, valid_meta$Group == "sarcopenia"])
  test <- wilcox.test(healthy, sarcopenia, exact = FALSE)
  data.frame(
    Gene = gene,
    Median_healthy = median(healthy),
    Median_sarcopenia = median(sarcopenia),
    P_value = test$p.value,
    stringsAsFactors = FALSE
  )
}))
valid_expression_stats$FDR <- p.adjust(valid_expression_stats$P_value, method = "BH")

# 10. Save output
fwrite(lasso_table, file.path(out_dir, "LASSO_lambda_1se_genes.csv"))
fwrite(
  data.frame(lambda_min = lasso_cv$lambda.min, lambda_1se = lasso_cv$lambda.1se),
  file.path(out_dir, "LASSO_lambda_values.csv")
)
fwrite(rf_error, file.path(out_dir, "RandomForest_error_by_tree.csv"))
fwrite(rf_importance, file.path(out_dir, "RandomForest_importance_all.csv"))
fwrite(rf_top10, file.path(out_dir, "RandomForest_top10_MeanDecreaseGini.csv"))
fwrite(as.data.frame(xgb_cv$evaluation_log), file.path(out_dir, "XGBoost_cross_validation.csv"))
fwrite(xgb_importance, file.path(out_dir, "XGBoost_importance_all.csv"))
fwrite(xgb_top10, file.path(out_dir, "XGBoost_top10_gain.csv"))
fwrite(consensus_table, file.path(out_dir, "Consensus_candidate_genes.csv"))
fwrite(
  data.frame(Term = names(coef(lrm_fit)), Coefficient = as.numeric(coef(lrm_fit))),
  file.path(out_dir, "Three_gene_lrm_coefficients.csv")
)
fwrite(
  data.frame(Statistic = names(lrm_fit$stats), Value = as.numeric(lrm_fit$stats)),
  file.path(out_dir, "Three_gene_lrm_model_statistics.csv")
)
fwrite(
  data.frame(
    Sample = rownames(x),
    Group = train_meta$Group,
    Linear_predictor = train_linear_predictor,
    Probability = train_probability
  ),
  file.path(out_dir, "Training_model_probabilities.csv")
)
fwrite(calibration_table, file.path(out_dir, "Bootstrap_calibration_B1000.csv"))
fwrite(decision_curve, file.path(out_dir, "Training_decision_curve_table.csv"))
fwrite(model_roc_summary, file.path(out_dir, "Three_gene_training_ROC_summary.csv"))
fwrite(rbind(train_roc, valid_roc), file.path(out_dir, "Single_gene_ROC_summary.csv"))
fwrite(valid_expression_stats, file.path(out_dir, "GSE111016_expression_statistics.csv"))

saveRDS(lasso_cv, file.path(out_dir, "LASSO_cv_model.rds"))
saveRDS(rf_fit, file.path(out_dir, "RandomForest_model.rds"))
saveRDS(xgb_fit, file.path(out_dir, "XGBoost_model.rds"))
saveRDS(lrm_fit, file.path(out_dir, "Three_gene_lrm_model.rds"))
saveRDS(nomogram_fit, file.path(out_dir, "Three_gene_nomogram.rds"))
saveRDS(calibration_fit, file.path(out_dir, "Bootstrap_calibration_B1000.rds"))

options(datadist = old_datadist)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
