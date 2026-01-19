setwd("~/Desktop/StJudeJCO")
library(biomaRt)
library(DESeq2)
library(dplyr)
library(nestedcv)
library(randomForest)
library(glmnet)
library(e1071)
library(xgboost)
library(pROC)
library(caret)
library(dplyr)
library(Matrix)
library(mixOmics)
library(e1071)
library(ggplot2)
#library(TCGAbiolinks)
library(readxl)
library(SummarizedExperiment)

# Load expression and mutation data
exp <- readRDS("Data/vst_norm_exp.rds")
mut <- readRDS("Data/genome_alterations.rds")
sv <- readRDS("Data/Mut_SV_only.rds")
clin <- readRDS("Data/merged_clinical_data.rds")
clin <- unique(clin)
clean_ids <- sub("-[A-Z0-9]+-A$", "", colnames(exp))
colnames(exp) <- clin$PatientID[match(clean_ids, clin$case_submitter_id)]

# Load additional clinical and expression data
clin.follow.up <- read.delim("Data/follow_up.tsv", header = TRUE, sep = "\t")
clin.follow.up <- unique(clin.follow.up)
data <- readRDS("Data/tcga_summarized_exp.rds") # This 'data' object is a SummarizedExperiment
clin <- read_excel("Data/ds_jco.23.02238.xlsx", sheet = 4)
colnames(clin) <- clin[2,]
clin <- clin[-c(1:2),]

# Prepare expression_matrix from 'data' SummarizedExperiment object
expression_matrix <- assay(data) # Extract gene expression matrix

smallestGroupSize <- 10
keep <- rowSums(expression_matrix >= 10) >= smallestGroupSize
expression_matrix <- expression_matrix[keep,]

# Transpose the matrix to have samples as rows and genes as columns
expression_matrix <- t(expression_matrix)

exp.log <- readRDS("Data/vst_norm_exp.rds") # This seems to be the same as 'exp' loaded earlier, consider if one can be removed if they are identical.
exp.log <- t(exp.log)

ids.exp <- rownames(exp.log)
ids.exp <- sub("^[^-]+-([^-]+)-.*$", "\\1", ids.exp)
clin.exp <- clin[match(ids.exp, clin$USI),]
clin.follow.up$USI <- sub("^[^-]+-([^-]+).*$", "\\1", clin.follow.up$case_submitter_id)
meta_data <- as.data.frame(data@colData)

# Filter for coding genes
row_metadata <- as.data.frame(rowData(data))
row_metadata <- row_metadata %>% filter(gene_type == "protein_coding")
exp.log <- as.data.frame(exp.log)
exp.filt <- exp.log[, na.omit(match(row_metadata$gene_name, colnames(exp.log)))]
rownames(exp.filt) <- sub("^[^-]+-([^-]+)-.*", "\\1", rownames(exp.filt))
rownames(expression_matrix) <- sub("^[^-]+-([^-]+)-.*", "\\1", rownames(expression_matrix))

prepare_cox_data <- function(subtype, clin_exp, exp_log, clin_follow_up) {
  # Filter for the specific subtype
  ind_subtype <- grep(subtype, clin_exp$Subtype)
  exp_subtype <- exp_log[ind_subtype, ]
  # Extract patient IDs
  ids_subtype <- sub("^[^-]+-([^-]+)-.*$", "\\1", rownames(exp_subtype))
  # Match with clinical follow-up and expression data
  clin_follow_subtype <- clin_follow_up[match(ids_subtype, clin_follow_up$USI), ]
  clin_exp_subtype <- clin_exp[match(ids_subtype, clin_exp$USI), ]
  
  # Add patient IDs to clinical follow-up
  clin_follow_subtype$USI <- ids_subtype
  
  # Merge clinical follow-up and expression data
  clin_combined <- merge(clin_follow_subtype, clin_exp_subtype, by = "USI", all.x = TRUE)
  
  # Create "Time" column
  clin_combined$Time <- NA
  clin_combined$Time[clin_combined$days_to_recurrence == "'--"] <-
    as.numeric(as.character(clin_combined$days_to_progression_free[clin_combined$days_to_recurrence == "'--"]))
  clin_combined$Time[clin_combined$days_to_recurrence != "'--"] <-
    as.numeric(as.character(clin_combined$days_to_recurrence[clin_combined$days_to_recurrence != "'--"]))
  
  # Create "Status" column
  clin_combined$Status <- 0
  clin_combined$Status[clin_combined$Relapse == "Relapse"] <- 1
  # Scale expression data
  exp_scaled <- scale(exp_subtype)
  
  # Combine clinical and expression data
  cox_df <- cbind(clin_combined$Time, clin_combined$Status, exp_scaled)
  colnames(cox_df)[1:2] <- c("Time", "Status")
  
  # Remove NAs and transpose data
  cox_df <- na.omit(t(cox_df))
  cox_df <- t(cox_df)
  out <- list()
  out$Cox.df <- cox_df
  out$Clin.merged <- clin_combined
  return(out)
}

split_data_mut <- function(subtype, exp.raw, exp.vst, clin_df, mut.df, sv.df, clin.follow.up, split_ratio = 0.7) {
  set.seed(2025)
  
  # Subset and merge clinical data
  clin.sub <- clin_df %>% filter(Subtype == subtype)
  clin.sub <- merge(clin.sub, clin.follow.up, by = "USI", all.x = TRUE)
  clin.sub <- as.data.frame(clin.sub)
  
  # Clean recurrence/progression fields
  clin.sub$days_to_recurrence <- as.character(clin.sub$days_to_recurrence)
  clin.sub$days_to_progression_free <- as.character(clin.sub$days_to_progression_free)
  clin.sub$days_to_recurrence[clin.sub$days_to_recurrence == "--"] <- NA
  clin.sub$days_to_recurrence <- as.numeric(clin.sub$days_to_recurrence)
  clin.sub$days_to_progression_free <- as.numeric(clin.sub$days_to_progression_free)
  
  clin.sub$Time <- ifelse(is.na(clin.sub$days_to_recurrence),
                          clin.sub$days_to_progression_free,
                          clin.sub$days_to_recurrence)
  
  clin.sub$Status <- 0
  clin.sub$Status[clin.sub$Relapse == "Relapse"] <- 1
  
  # Remove patients with missing Time
  clin.sub <- clin.sub[!is.na(clin.sub$Time), ]
  
  # Match sample IDs
  raw.samples <- rownames(exp.raw)
  vst.samples <- rownames(exp.vst)
  clin.samples <- clin.sub$USI
  mut.samples <- as.character(na.omit(clin.sub$USI[match(rownames(mut.df), clin.sub$PatientID)]))
  sv.samples <- clin.sub$USI[match(rownames(sv.df), clin.sub$PatientID)]
  
  common.samples <- Reduce(intersect, list(raw.samples, vst.samples, clin.samples, mut.samples, sv.samples))
  
  # Subset and align all datasets to common samples
  sample.mut <- clin.sub$PatientID[match(common.samples, clin.sub$USI)]
  exp.raw <- exp.raw[common.samples, , drop = FALSE]
  exp.vst <- exp.vst[common.samples, , drop = FALSE]
  rownames(clin.sub) <- clin.sub$USI
  clin.sub <- clin.sub[common.samples, , drop = FALSE]
  mut.sub <- mut.df[sample.mut, , drop = FALSE]
  sv.sub <- sv.df[sample.mut, , drop = FALSE]
  
  mut.comb <- cbind(mut.sub, sv.sub)
  ind.rm <- which(is.na(colnames(mut.comb)))
  if (length(ind.rm) > 0) { # Only try to remove if there are NAs
    mut.comb <- mut.comb[, -ind.rm, drop = FALSE]
  }
  rownames(mut.comb) <- common.samples # Ensure USI rownames
  
  # Split
  num_samples <- length(common.samples)
  train_indices <- sample(seq_len(num_samples), size = floor(split_ratio * num_samples))
  test_indices <- setdiff(seq_len(num_samples), train_indices)
  
  # Assign split with names
  out <- list()
  out$train_exp_raw <- exp.raw[train_indices, , drop = FALSE]
  out$test_exp_raw <- exp.raw[test_indices, , drop = FALSE]
  
  out$train_exp_vst <- exp.vst[train_indices, , drop = FALSE]
  out$test_exp_vst <- exp.vst[test_indices, , drop = FALSE]
  
  out$train_clin <- clin.sub[train_indices, , drop = FALSE]
  out$test_clin <- clin.sub[test_indices, , drop = FALSE]
  
  out$train_mut <- mut.comb[train_indices, , drop = FALSE]
  out$test_mut <- mut.comb[test_indices, , drop = FALSE]
  
  # All outputs have rownames = USI
  stopifnot(all(sapply(out, function(x) is.null(rownames(x)) == FALSE)))
  
  return(out)
}

test <- split_data_mut("ETV6::RUNX1", expression_matrix, exp.filt, clin.exp, mut, sv, clin.follow.up, 0.8)

exp.train<-test$train_exp_vst
exp.train<-as.matrix(exp.train)
exp.train<-t(exp.train)
clin.train <-test$train_clin
clin.train$Risk<-as.factor(clin.train$Risk)
clin.train$Relapse<-as.factor(clin.train$Relapse)
class.relapse<-as.factor(test$train_clin$Relapse)
names(class.relapse)<-colnames(exp.train)

exp.test <- test$test_exp_vst
exp.test <- as.matrix(exp.test)
exp.test <- t(exp.test)

# -------------------------------------------------------------------------
library(fgsea)
library(GSVA)
pathway.biocarta <- gmtPathways("Data/c2.cp.biocarta.v2024.1.Hs.symbols.gmt")
gsvaPar.train <- gsvaParam(exp.train,pathway.biocarta)
gsva.es.train <- gsva(gsvaPar.train,verbose=FALSE)

gsvaPar.test <- gsvaParam(exp.test,pathway.biocarta)
gsva.es.test <- gsva(gsvaPar.test,verbose=FALSE)


gsva.es.train <- t(gsva.es.train)
relapse_status_for_gsva_train <- class.relapse[rownames(gsva.es.train)]
gsva.train <- as.data.frame(gsva.es.train) 
gsva.train$Relapse_Status <- relapse_status_for_gsva_train



#decoupleR
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(TCGAbiolinks)
exp.train
net <- readRDS("Data/net.RDS")

sample_acts.train <- decoupleR::run_ulm(
  mat = exp.train,
  net = net,
  .source = 'source',
  .target = 'target',
  .mor = 'mor',
  minsize = 5 # Default is 5, you can adjust this
)

decoupler_activities_train_wide <- sample_acts.train %>%
  dplyr::select(source, condition, score) %>% # Keep only the relevant columns
  tidyr::pivot_wider(names_from = source, values_from = score) %>% # Pivot 'source' names to new columns
  tibble::column_to_rownames(var = "condition") %>% # Set sample IDs as row names
  as.data.frame() # Convert to a standard data.frame, which randomForest prefers
relapse_status_for_decoupler_train <- class.relapse[rownames(decoupler_activities_train_wide)]
decoupler_and_relapse.train_df <- decoupler_activities_train_wide
decoupler_and_relapse.train_df$Relapse_Status <- relapse_status_for_decoupler_train
colnames(decoupler_and_relapse.train_df) <- gsub("-", "_", colnames(decoupler_and_relapse.train_df))

# Boruta ------------------------------------------------------------------
library(Boruta)

#Boruta with rna-seq
traindata <- t(exp.train)
relapse_status.train <- class.relapse[rownames(traindata)]
traindata <- as.data.frame(traindata)
traindata$Relapse_Status <- relapse_status.train

set.seed(12345)
exp.boruta.train <- Boruta(Relapse_Status ~ ., data = traindata, doTrace = 2)
exp.final.boruta <- TentativeRoughFix(exp.boruta.train)
exp.features <- getSelectedAttributes(exp.final.boruta, withTentative = FALSE)
exp.boruta.df <- attStats(exp.final.boruta)

data_for_exp_rf <- traindata %>%
  dplyr::select(dplyr::all_of(exp.features), Relapse_Status)

exp.rf_genes <- randomForest(Relapse_Status ~ ., data = data_for_exp_rf, proximity = TRUE,classwt=c("None" = 1, "Relapse" = 3.06))
print(exp.rf_genes)

testdata <- t(exp.test)
data_for_exp_rf_test <- testdata %>%
  as.data.frame() %>%
  dplyr::select(dplyr::all_of(exp.features))
test_sample_ids <- rownames(data_for_exp_rf_test)
actual_relapse_status_test <- test$test_clin$Relapse[match(test_sample_ids, test$test_clin$USI)]
test_data_for_rf <- cbind(data_for_exp_rf_test, Relapse_Status = as.factor(actual_relapse_status_test))
colnames(test_data_for_rf)[5] <- "Relapse_Status"
test_data_for_rf$Relapse_Status <- as.factor(test_data_for_rf$Relapse_Status)
predicted_relapse_status <- predict(exp.rf_genes, newdata = test_data_for_rf[, -which(colnames(test_data_for_rf) == "Relapse_Status")])
confusion_matrix <- table(Predicted = predicted_relapse_status, Actual = test_data_for_rf$Relapse_Status)
print(confusion_matrix)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
cat("Accuracy:", accuracy, "\n")

if (length(levels(test_data_for_rf$Relapse_Status)) == 2) {
  library(pROC)
  actual_numeric <- as.numeric(test_data_for_rf$Relapse_Status) - 1
  predicted_prob <- predict(exp.rf_genes, newdata = test_data_for_rf[, -which(colnames(test_data_for_rf) == "Relapse_Status")], type = "prob")[,2]
  roc_obj <- roc(actual_numeric, predicted_prob)
  plot(roc_obj, main = "ROC Curve (Gene expression)")
  cat("AUC:", auc(roc_obj), "\n")
}



#Boruta with pathways
set.seed(12345)
gsva.boruta.train <- Boruta(Relapse_Status~., data = gsva.train, doTrace = 2)
gsva.final.boruta <- TentativeRoughFix(gsva.boruta.train)
print(gsva.final.boruta)
gsva.features <- getSelectedAttributes(gsva.final.boruta, withTentative = F)
gsva.boruta.df <- attStats(gsva.final.boruta)

data_for_gsva_rf <- gsva.train %>%
  dplyr::select(all_of(gsva.features), Relapse_Status)
gsva.rf <- randomForest(Relapse_Status~.,data = data_for_gsva_rf, proximity = TRUE,classwt=c("None" = 1, "Relapse" = 3.06))
print(gsva.rf)

testdata.gsva <- t(gsva.es.test)
data_for_gsva_rf_test <- testdata.gsva %>%
  as.data.frame() %>%
  dplyr::select(dplyr::all_of(gsva.features))
test_sample_ids <- rownames(testdata.gsva)
actual_relapse_status_test <- test$test_clin$Relapse[match(test_sample_ids, test$test_clin$USI)]
test_data_for_rf <- cbind(data_for_gsva_rf_test, Relapse_Status = as.factor(actual_relapse_status_test))
colnames(test_data_for_rf)[7] <- "Relapse_Status"
test_data_for_rf$Relapse_Status <- as.factor(test_data_for_rf$Relapse_Status)
predicted_relapse_status <- predict(gsva.rf , newdata = test_data_for_rf[, -which(colnames(test_data_for_rf) == "Relapse_Status")])
confusion_matrix <- table(Predicted = predicted_relapse_status, Actual = test_data_for_rf$Relapse_Status)
print(confusion_matrix)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
cat("Accuracy:", accuracy, "\n")

if (length(levels(test_data_for_rf$Relapse_Status)) == 2) {
  library(pROC)
  actual_numeric <- as.numeric(test_data_for_rf$Relapse_Status) - 1
  predicted_prob <- predict(gsva.rf, newdata = test_data_for_rf[, -which(colnames(test_data_for_rf) == "Relapse_Status")], type = "prob")[,2]
  roc_obj <- roc(actual_numeric, predicted_prob)
  plot(roc_obj, main = "ROC Curve (GSVA Pathway)")
  cat("AUC:", auc(roc_obj), "\n")
}

#Boruta with TF
set.seed(12345)
TF.train <- decoupler_and_relapse.train_df
TF.boruta.train <- Boruta(Relapse_Status~., data = TF.train, doTrace = 2)
TF.final.boruta <- TentativeRoughFix(TF.boruta.train)
print(TF.final.boruta)
TF.features <- getSelectedAttributes(TF.final.boruta, withTentative = F)
TF.boruta.df <- attStats(TF.final.boruta)

data_for_TF_rf <- TF.train %>%
  dplyr::select(all_of(TF.features), Relapse_Status)
TF.rf <- randomForest(Relapse_Status~.,data = data_for_TF_rf, proximity = TRUE,classwt=c("None" = 1, "Relapse" = 3.06))
print(TF.rf)

sample_acts.test <- decoupleR::run_ulm(
  mat = exp.test,
  net = net,
  .source = 'source',
  .target = 'target',
  .mor = 'mor',
  minsize = 5 # Default is 5, you can adjust this
)
decoupler_activities_train_wide <- sample_acts.test%>%
  dplyr::select(source, condition, score) %>% 
  tidyr::pivot_wider(names_from = source, values_from = score) %>% 
  tibble::column_to_rownames(var = "condition") %>% 
  as.data.frame() 
data_for_TF_test <- decoupler_activities_train_wide %>%
  as.data.frame() %>%
  dplyr::select(dplyr::all_of(TF.features))
test_sample_ids <- rownames(data_for_TF_test)
actual_relapse_status_test <- test$test_clin$Relapse[match(test_sample_ids, test$test_clin$USI)]
test_data_for_rf <- cbind(data_for_TF_test, Relapse_Status = as.factor(actual_relapse_status_test))
predicted_relapse_status <- predict(TF.rf , newdata = test_data_for_rf[, -which(colnames(test_data_for_rf) == "Relapse_Status")])
confusion_matrix <- table(Predicted = predicted_relapse_status, Actual = test_data_for_rf$Relapse_Status)
print(confusion_matrix)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
cat("Accuracy:", accuracy, "\n")

if (length(levels(test_data_for_rf$Relapse_Status)) == 2) {
  library(pROC)
  actual_numeric <- as.numeric(test_data_for_rf$Relapse_Status) - 1
  predicted_prob <- predict(TF.rf, newdata = test_data_for_rf[, -which(colnames(test_data_for_rf) == "Relapse_Status")], type = "prob")[,2]
  roc_obj <- roc(actual_numeric, predicted_prob)
  plot(roc_obj, main = "ROC Curve (Transcription Factor)")
  cat("AUC:", auc(roc_obj), "\n")
}


#Combine Expression and Pathway for feature selection

gene.train <- traindata %>%
  dplyr::select(all_of(exp.features))
relapse_status.train <- traindata$Relapse_Status
gsva.train.sel <- gsva.train %>%
  dplyr::select(all_of(gsva.features))
combined.train <- cbind(gene.train, gsva.train.sel)
combined.train$Relapse_Status <- relapse_status.train

set.seed(12345)
combined.rf <- randomForest(Relapse_Status ~ ., data = combined.train, proximity = TRUE, classwt = c("None" = 1, "Relapse" = 3.06))
print(combined.rf)

gene.test <- testdata %>%
  as.data.frame() %>%
  dplyr::select(dplyr::all_of(exp.features))
gsva.test <- testdata.gsva %>%
  as.data.frame() %>%
  dplyr::select(dplyr::all_of(gsva.features))
combined.test <- cbind(gene.test, gsva.test)
test_sample_ids <- rownames(combined.test)
actual_relapse_status_test <- test$test_clin$Relapse[match(test_sample_ids, test$test_clin$USI)]
combined.test$Relapse_Status <- as.factor(actual_relapse_status_test)

predicted_relapse_combined <- predict(combined.rf,
                                      newdata = combined.test[, -which(colnames(combined.test) == "Relapse_Status")])
confusion_matrix_combined <- table(Predicted = predicted_relapse_combined,
                                   Actual = combined.test$Relapse_Status)
print(confusion_matrix_combined)
accuracy_combined <- sum(diag(confusion_matrix_combined)) / sum(confusion_matrix_combined)
cat("Combined Accuracy:", accuracy_combined, "\n")

#Correlation network
build_correlation_network <- function(data, feat_data, cor_threshold = 0.6, only_positive_edges = TRUE) {
  library(dplyr)
  library(igraph)
  
  expr <- t(data$train_exp_vst)
  clin <- data$train_clin
  
  # Extract genes from feat_data
  seed_genes <- feat_data 
  # Subset expression to sensitive ("R") samples
  sensitive_expr <- expr
  sensitive_expr <- t(sensitive_expr)
  # Compute correlation matrix for each seed gene
  cor_matrix <- sapply(seed_genes, function(seed) {
    if (seed %in% colnames(sensitive_expr)) {
      cor(sensitive_expr, sensitive_expr[, seed], method = "pearson", use = "pairwise.complete.obs")
    } else {
      rep(NA, ncol(sensitive_expr))
    }
  })
  rownames(cor_matrix) <- colnames(sensitive_expr)
  # Identify highly correlated genes
  high_cor_genes <- unique(rownames(cor_matrix)[apply(abs(cor_matrix), 1, max, na.rm = TRUE) > cor_threshold])
  
  network_genes <- union(seed_genes, high_cor_genes)
  sub_expr <- sensitive_expr[, network_genes, drop = FALSE]
  
  # Full correlation network on selected genes
  network_cor <- cor(sub_expr, method = "pearson", use = "pairwise.complete.obs")
  
  # Create edge list from upper triangle
  edge_list <- which(abs(network_cor) > 0.6 & lower.tri(network_cor), arr.ind = TRUE)
  edges <- data.frame(
    GENE1 = colnames(network_cor)[edge_list[,1]],
    GENE2 = colnames(network_cor)[edge_list[,2]],
    WEIGHT = network_cor[edge_list]
  )
  
  # Optionally keep only positive edges
  if (only_positive_edges) {
    edges <- edges[edges$WEIGHT > 0, ]
  }
  
  # Save to CSV
  
  # Create igraph object
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # Detect communities using Louvain algorithm
  modules <- cluster_louvain(g, weights = E(g)$WEIGHT)
  V(g)$module <- membership(modules)
  
  # Create module assignment table
  module_df <- data.frame(
    Gene = V(g)$name,
    Module = V(g)$module
  )
  
  # Optional plot
  plot(modules, g,
       vertex.label = V(g)$name,
       vertex.size = 5,
       vertex.label.cex = 0.7,
       main = "Gene Co-expression Network Modules (Louvain)")
  
  # Return outputs
  return(list(
    graph = g,
    modules = modules,
    module_table = module_df,
    module_list = split(module_df$Gene, module_df$Module)
  ))
}

genes.filt <- exp.features
etv6.cor<-build_correlation_network(test,genes.filt,cor_threshold = 0.6,only_positive_edges = TRUE)

#Module enrichment
run_module_enrichment <- function(result, expression_data, min_genes = 10, q_cutoff = 0.2) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  # Extract modules
  module_df <- result$module_table
  modules_list <- split(module_df$Gene, module_df$Module)
  
  # Filter modules by size
  module_sizes <- lengths(modules_list)
  module.sig <- which(module_sizes > min_genes)
  
  # Gene universe from expression matrix
  universe_genes <- rownames(expression_data)
  
  # Helper function to run GO enrichment
  run_enrichment <- function(gene_list, universe) {
    enrichGO(gene         = gene_list,
             OrgDb        = org.Hs.eg.db,
             keyType      = "SYMBOL",
             ont          = "BP",
             pAdjustMethod = "BH",
             minGSSize = 10,
             universe     = universe,
             qvalueCutoff = q_cutoff,
             readable     = TRUE)
  }
  
  # Apply enrichment to selected modules
  enrich_results <- lapply(modules_list[module.sig], run_enrichment, universe = universe_genes)
  
  # Name results by module
  names(enrich_results) <- paste0("Module_", names(modules_list[module.sig]))
  
  return(enrich_results)
}
pathway.etv6<-run_module_enrichment(etv6.cor,exp.train)

##### tell which cohort did the model predict correctly and which didn't


