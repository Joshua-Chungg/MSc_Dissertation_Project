#setwd("/home/regmjkc/Scratch/StJudeJCO")
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
library(FactoMineR) # For PCA computation, commented out PCA code so can be removed
library(factoextra) # For visualization, commented out PCA code so can be removed
library(SuperLearner) # Added to the list as it's used in feat_select_mixOmics_nestedCV
library(BiocParallel) # Added as it's used in tune.block.splsda


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
clin.clust <- readRDS("Data/clin.etv6.clust.rds")

# Filter for coding genes
row_metadata <- as.data.frame(rowData(data))
row_metadata <- row_metadata %>% filter(gene_type == "protein_coding")
exp.log <- as.data.frame(exp.log)
exp.filt <- exp.log[, na.omit(match(row_metadata$gene_name, colnames(exp.log)))]
rownames(exp.filt) <- sub("^[^-]+-([^-]+)-.*", "\\1", rownames(exp.filt))
rownames(expression_matrix) <- sub("^[^-]+-([^-]+)-.*", "\\1", rownames(expression_matrix))
expression_matrix<-readRDS("Data/exp_raw_counts.rds")
expression_matrix<-t(expression_matrix)
split_data_mut <- function(subtype, exp.raw, exp.vst, clin_df, mut.df, sv.df, clin.follow.up, split_ratio = 0.7) {
  set.seed(2025)
  
  # Subset and merge clinical data
  clin.sub <- clin_df %>% filter(Subtype == subtype)
  #clin.sub <- merge(clin.sub, clin.follow.up, by = "USI", all.x = TRUE)
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

  train_indices<-caret::createDataPartition(clin.sub$Cluster,p=split_ratio,list=FALSE)
  test_indices<-setdiff(1:nrow(clin.sub),train_indices)
  
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

test <- split_data_mut("ETV6::RUNX1", expression_matrix, exp.filt, clin.clust, mut, sv, clin.follow.up, 0.8)
# Filter for coding genes

exp.filt<-test$train_exp_raw
library(tidyr)
housekeeping<- c("18SRNA", "ACTB", "B2M", "GUSB", "GAPDH", "HPRT1", "MT-ATP6", "RPS17")
housekeeping_data <- exp.filt[, na.omit(match(housekeeping, colnames(exp.filt)))]

housekeeping_long <- as.data.frame(housekeeping_data) %>%
       tibble::rownames_to_column(var = "Sample") %>%
       pivot_longer(
            cols = -Sample,
             names_to = "Gene",
             values_to = "Expression"
         )
 ggplot(housekeeping_long, aes(x = Gene, y = Expression, fill = Gene)) +
       geom_boxplot() +
       labs(title = "Distribution of Housekeeping Gene Expression",
                       x = "Housekeeping Gene",
                       y = "Expression Value") +
       theme_minimal() +
       theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
       guides(fill = "none") # This removes the legend for a cleaner look
 library(matrixStats)
 hk_expr_pseudo <- housekeeping_data + 1
 
 # 4. Calculate geometric mean per sample (columns)
 geo_means <- exp(rowMeans(log(hk_expr_pseudo)))
 
 # 5. Normalize all gene expression values (add pseudocount here too!)
expr_matrix_pseudo <- test$train_exp_raw + 1  # very important!
expr_matrix_norm <- sweep(expr_matrix_pseudo, 1, geo_means, FUN = "/")
expr_matrix_log <- log2(expr_matrix_norm) 
mean_expr <- colMeans(expr_matrix_log)

exp_filtered<-expr_matrix_log[,which(mean_expr >= -3)]
exp_filtered <- t(exp_filtered)

gene_var <- apply(exp_filtered, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:5000]

# 3. Subset to top 5000
expr_top5000 <- exp_filtered[top_genes, ]
expr_scaled<-scale(t(expr_top5000))
exp.cox.train<-as.data.frame(expr_scaled)
exp.cox.train <-  cbind(test$train_clin %>% dplyr::select(c("Time","Status")),exp.cox.train)
library(survival)
cox_univariate_analysis <- function(cox.train.exp) {
  cox.train.exp<-as.data.frame(cox.train.exp)
  # Extract covariate names (excluding first two columns: Time, Status)
  covariates <- colnames(cox.train.exp)[-c(1,2)]
  
  # Rename gene columns systematically
  colnames(cox.train.exp)[-c(1,2)] <- paste0("Gene", seq_along(covariates))
  
  # Create a conversion table mapping original gene names to new IDs
  conv.table <- cbind(genes = covariates, ID = paste0("Gene", seq_along(covariates)))
  
  # Create univariate Cox regression formulas for each gene
  univ_formulas <- sapply(colnames(cox.train.exp)[-c(1,2)], 
                          function(x) as.formula(paste('Surv(Time, Status) ~', x)))
  cox.train.exp<-as.data.frame(cox.train.exp)
  # Fit Cox models for each gene
  univ_models <- lapply(univ_formulas, function(x) { 
    coxph(x, data = cox.train.exp) 
  })
  
  # Extract summary statistics from each model
  univ_results <- lapply(univ_models, function(x) {
    x_summary <- summary(x)
    p.value <- signif(x_summary$wald["pvalue"], digits = 2)
    wald.test <- signif(x_summary$wald["test"], digits = 2)
    beta <- signif(x_summary$coef[1], digits = 2)  # Coefficient beta
    HR <- signif(x_summary$coef[2], digits = 2)  # Hazard ratio (exp(beta))
    HR.confint.lower <- signif(x_summary$conf.int[,"lower .95"], 2)
    HR.confint.upper <- signif(x_summary$conf.int[,"upper .95"], 2)
    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
    
    res<-c(beta, HR, wald.test, HR.confint.lower, HR.confint.upper,
           p.value)
    names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", "HR_Lower","HR_Upper",
                  "p.value")
    return(res)
  })
  
  # Convert results to a data frame
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res <- as.data.frame(res)
  #rownames(res)<-conv.table[,1]
  # Convert p-values to numeric and apply Benjamini-Hochberg (BH) FDR correction
  res$p.value <- as.numeric(res$p.value)
  res$p.adjust <- p.adjust(res$p.value, method = "BH")
  res$Gene<-conv.table[,1][match(rownames(res),conv.table[,2])]
  # Return results
  return(list(results = res, conversion_table = conv.table))
}
cox.uni.etv6<-cox_univariate_analysis(exp.cox.train)
View(cox.uni.etv6[["results"]] %>% filter(p.value<0.01))

relapse <- data.frame(time=test$train_clin$Time,status=test$train_clin$Status)

cv_results <- list()
lambda_min_results <- numeric(5)
for (i in 1:5){
glm_model <- cv.glmnet(expr_scaled, clin.clust$Relapse, family = "cox", type.measure = "C",)
cv_results[[i]] <- cv_model
lambda_min_results[i] <- cv_model$lambda.min
}


#CVglmnet 



