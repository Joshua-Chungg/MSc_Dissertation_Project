library(dplyr)
library(ggplot2)
library(caret)
library(vip)
library(randomForest)
library(sva)
library(tidymodels)
library(edgeR)
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(survminer)
library(boot)
library(glmnet)
library(SGL)
library(rsample)


exp_raw_counts <- readRDS("/Users/JoshuaChung 1/Desktop/MSc Project/ALL DATA/exp_raw_counts.rds")
merged_clinical_data <- readRDS("/Users/JoshuaChung 1/Desktop/MSc Project/ALL DATA/merged_clinical_data.rds")

##Match colnames of the expression matrix to the USI column of the clinical data
#order the colnames
ind <- match(colnames(exp_raw_counts),merged_clinical_data$USI) 
#reorders
merged_clinical_data <- merged_clinical_data[na.omit(ind),] 
exp_raw_counts <- exp_raw_counts[,match(merged_clinical_data$USI,colnames(exp_raw_counts))]
##Filter for lowly expressed genes
dge <- DGEList(counts = exp_raw_counts)
keep <- filterByExpr(dge)
dge_filtered <- dge[keep, ]
exp_filtered <- dge_filtered$counts

##COMBAT-Seq batch effect correction
batch <- as.factor(merged_clinical_data$Batch)
count_matrix <- as.matrix(exp_filtered)
corrected_counts <- ComBat_seq(counts = count_matrix, batch = batch, group = NULL)
#Get the common samples
common_samples_combat <- intersect(colnames(exp_filtered), merged_clinical_data$USI)
# Subset and order the expression data
exp_combat_matched <- exp_filtered[, common_samples_combat]
# Subset and order the clinical data (and thus the batch variable) to match
clinical_combat_matched <- merged_clinical_data[match(common_samples_combat, merged_clinical_data$USI), ]
batch_combat_matched <- clinical_combat_matched$Batch
# Verify the order of sample IDs
identical(colnames(exp_combat_matched), as.character(clinical_combat_matched$USI))
# Now apply ComBat-seq with the matched data
count_matrix_combat <- as.matrix(exp_combat_matched)
corrected_counts <- ComBat_seq(counts = count_matrix_combat, batch = batch_combat_matched, group = NULL)

##Apply log(batch corrected counts +1) transformation
log_transformed_counts <- log10(corrected_counts + 1)

##Filtering data and scaling
total_counts <- rowSums(log_transformed_counts)
log_transformed_counts <- log_transformed_counts[total_counts >= 10, ]

#Take the top 2000 most variable genes
SDs=apply(log_transformed_counts,1,var )
topPreds=order(SDs,decreasing = TRUE)[1:2000]
nzv_log_transformed_counts=log_transformed_counts[topPreds,]

#Centering of data
nzv_log_transformed_counts <- scale(t(nzv_log_transformed_counts))


# ETV6 --------------------------------------------------------------------

mut <- readRDS("/Users/JoshuaChung 1/Desktop/MSc Project/Part1/ETV6 subtype/Mut_SV_only.rds")
gen_alt <- readRDS("/Users/JoshuaChung 1/Desktop/MSc Project/Part1/ETV6 subtype/genome_alterations.rds")
ind.rm<-which(is.na(match(rownames(gen_alt),clinical_combat_matched$PatientID)))
mut.filt<-gen_alt[-ind.rm,]
rownames(mut.filt)<-clinical_combat_matched$USI[match(rownames(mut.filt),clinical_combat_matched$PatientID)]
clin.etv6<- clinical_combat_matched %>% filter(Subtype=="ETV6::RUNX1")
exp.log.etv6<-log_transformed_counts[,match(clin.etv6$USI,colnames(log_transformed_counts))]
mut.etv6<-mut.filt[na.omit(match(clin.etv6$USI,rownames(mut.filt))),]
clin.etv6<-clin.etv6[match(rownames(mut.etv6),clin.etv6$USI),]  
exp.log.etv6<-exp.log.etv6[,match(rownames(mut.etv6),colnames(exp.log.etv6))]  
mut.etv6[which(mut.etv6==2,arr.ind = T)]<- 1

# scale
etv6.counts <- rowSums(exp.log.etv6)
exp.log.etv6 <- exp.log.etv6[etv6.counts >= 10, ]

# select top 2k variance genes
SDs=apply(exp.log.etv6,1,var )
topPreds=order(SDs,decreasing = TRUE)[1:2000]
top.genes.etv6=exp.log.etv6[topPreds,]

#Centering of data
top.genes.etv6 <- scale(t(top.genes.etv6))


# Machine Learning -----------------------------------------------------

mut.dat<- mut.etv6
clin.dat<- clin.etv6
exp.dat<- exp.log.etv6

#filter for 20 mutations
mut.dat<-mut.dat[,which(colSums(mut.dat)>20)]
#combine relapse status to mut
mut.df<-  mut.dat
mut.df[,"Relapse"]<- clin.dat$Relapse

#split data (mut)
set.seed(123)
mut.split<- initial_split(mut.df,prop=0.8,strata= "Relapse")
mut.train<- training(mut.split)
mut.test<- testing(mut.split)

#scale expression matrix
#exp.dat.scale<- scale(t(exp.dat))
#combine relapse status to mut
#exp.df<- exp.dat.scale
#exp.df[,"Relapse"]<- clin.dat$Relapse



# Logistic Regression -----------------------------------------------------








# Random Forest -----------------------------------------------------------
mut.train$Relapse <- as.factor(mut.train$Relapse)
mut.rf <- randomForest(Relapse~., data = mut.train, proximity = TRUE, mtry = sqrt(190), ntree = 500)
plot(mut.rf)

#mtry=6 has OOB error of 0.23
t <- tuneRF(mut.train[,-191], mut.train[,191],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)

#Normal distribution of tree nodes, with bulk of the trees with 100-110 nodes
hist(treesize(mut.rf),
     main = "No. of Nodes for the Trees",
     col = "green")

#Decision tree visualisation (NOT WORKING)
#tree_spec <- decision_tree() %>%
  #set_engine("rpart") %>%
  #set_mode("classification")
#tree_fit <- tree_spec %>%
  #fit(Relapse~.,data = mut.train)
#rpart.plot(tree_fit$fit,type = 4,extra = 101,under = TRUE,cex=0.8,box.palette = "auto")

#nonsynonymous_SNV_PAX5 3.68, good predictor of relapse
#nonsynonymous_SNV_ETV6 2.14
#HETDEL.5MB_VPREB1 1.90
#HETDEL_NR3C1 1.89
#HETDEL_ACTB 1.87
#nonsynonymous_SNV_NSD2 1.84
#HETDEL_CTCF 1.79
#GAIN_LINC00649 1.70
#HETDEL_CHD4 1.67
#HETDEL.5MB_CDKN1B 1.63
varImpPlot(mut.rf,
           sort=T,
           n.var = 40,
           main="Top 40-Variable Importance")
mut.real.importance<-importance(mut.rf)

partialPlot(mut.rf,mut.train,nonsynonymous_SNV_PAX5,"Relapse")

#shuffle relapse column 100 times and keep the rest of the mutation matrix fixed
pre.shuffle<- mut.train
perm<- list()
for (i in 1:190){
  set.seed(i)
  shuffled<- transform(pre.shuffle,"Relapse"=sample(Relapse))
  shuffled.rf<- randomForest(Relapse~., data = shuffled, proximity = FALSE, mtry = sqrt(190), ntree = 150)
  shuffled.imp<-importance(shuffled.rf)
  perm[[i]]<- shuffled.imp
}

randomforest.imp<- do.call(cbind,perm)
perm.percent<- list()
for (i in 1:nrow(mut.real.importance)){
  x<- mut.real.importance[i]-randomforest.imp[i,]
  num<- length(which(x >= 0))
  total<- 190
  percentage<- (num/total)*100
  perm.percent[[i]]<-percentage
  
}

sig.10<-colnames(mut.train)[which(signif(as.numeric(unlist(perm.percent)),2)<= 5)]
perc.10<-cbind(sig.10,unlist(perm.percent)[which(signif(unlist(perm.percent),2)<= 5)])
perc.10[,2]<-as.numeric(as.character(perc.10[,2]))


mut.filt.df<-mut.train[,sig.10]
mut_filt_df_USIs <- rownames(mut.filt.df)
matched_indices_in_clin_etv6 <- match(mut_filt_df_USIs, clin.etv6$USI)
mut.filt.df$Relapse <- clin.etv6$Relapse[matched_indices_in_clin_etv6]
mut.filt.df$Relapse <- as.factor(mut.filt.df$Relapse)
mut.filt.rf <- randomForest(Relapse~., data = mut.filt.df, proximity = TRUE, mtry = sqrt(9), ntree = 500)

reprtree:::plot.getTree(mut.filt.rf)

##combine mutation, age, MRD and risk and randomforest classification on training data
mut.filt.df$MRD<- clin.etv6$MRD[matched_indices_in_clin_etv6]
mut.filt.df$Risk<- clin.etv6$Risk[matched_indices_in_clin_etv6]
mut.filt.rf <- randomForest(Relapse~., data = mut.filt.df, proximity = TRUE, mtry = sqrt(9), ntree = 500)

reprtree:::plot.getTree(mut.filt.rf)

#also incorporate gene expression data
exp.dat<- t(exp.dat)#so that rows are USI and columns are genes
exp.dat<- as.data.frame(exp.dat)
exp.dat$USI <- rownames(exp.dat)
mut.filt.df_with_USI <- mut.filt.df
mut.filt.df_with_USI$USI <- rownames(mut.filt.df_with_USI)
mut.exp.df <- merge(mut.filt.df_with_USI, exp.dat, by = "USI", all.x = TRUE)
rownames(mut.exp.df) <- mut.exp.df$USI
mut.exp.df$USI <- NULL # Remove the USI column
names(mut.exp.df) <- make.names(names(mut.exp.df))
mut.exp.rf<- randomForest(Relapse~., data = mut.exp.df, proximity = TRUE, mtry = sqrt(13544), ntree = 500)

reprtree:::plot.getTree(mut.exp.rf)            


##plot variable important score (PCTP,LILRB3,INHBB, CCDC9B,EPPK1,RTKN,PLD1,CALN1,TAF8,PLD2)
varImpPlot(mut.exp.rf,
           sort=T,
           n.var = 10,
           main="Top 10-Variable Importance")
mut.real.importance<-importance(mut.exp.rf)



# -------------------------------------------------------------------------
set.seed(123)
top.genes.etv6 #expression matrix (filtered low-expressed genes, normalised and centered)
mut.etv6 #mutational profile
clin.etv6 #clinical data

#split data
clin.df <- as.data.frame(clin.etv6)
clin.df <- clin.df[,c("USI","Sex","Age","Karyotype","Risk","Relapse","MRD","days_to_follow_up","days_to_recurrence",
  "days_to_progression_free")]
clin.df$Karyotype <- substr(clin.df$Karyotype, start = 1, stop = 2) #take the first two numbers of the Karyotype column
clin.df$Relapse <- ifelse(clin.df$Relapse=="Relapse",1,0) #change the relapse column into binary

#combine days to progression-free and days to recurrence and filter out NA
clin.df$days_to_relapse <- ifelse((clin.df$days_to_progression_free=="'--"),clin.df$days_to_recurrence,clin.df$days_to_progression_free)
clin.df <- clin.df[clin.df$days_to_relapse != "'--", ]

clin.split <- initial_split(clin.df,prop=0.8,strata= "Relapse")
clin.train<- training(clin.split)
clin.test<- testing(clin.split)

#match the expression matrix and mutational profile to the training data
common_usis <- intersect(clin.train$USI, rownames(mut.etv6))
common_usis_ordered <- clin.train$USI[clin.train$USI %in% common_usis]
mut.train <- mut.etv6[match(common_usis_ordered, rownames(mut.etv6)), ]

common_usis <- intersect(clin.train$USI, rownames(top.genes.etv6))
common_usis_ordered <- clin.train$USI[clin.train$USI %in% common_usis]
exp.train <- top.genes.etv6[match(common_usis_ordered, rownames(top.genes.etv6)), ]

cox_dat <- data.frame(Relapse = clin.train$Relapse,days_to_relapse = clin.train$days_to_relapse,USI=clin.train$USI)
exp.train_USI <- as.data.frame(exp.train)
exp.train_USI$USI<- rownames(exp.train)
cox_dat <- merge(cox_dat,exp.train_USI,by="USI")
cox_dat$days_to_relapse <- as.numeric(cox_dat$days_to_relapse)

#cox regression to get coefficient and p value for genes
num <- ncol(cox_dat)-3
res.cox <- matrix(NA,nrow=num,ncol = 5)
for (i in 1:num){
  cox <- coxph(Surv(days_to_relapse,Relapse)~cox_dat[,(3+i)],data = cox_dat)
  cox_summary <- summary(cox)
  res.cox[i,]<- cox_summary$coefficients
  }
res.cox <- as.data.frame(res.cox)
colnames(res.cox)[1] <- "coef"
colnames(res.cox)[2] <- "exp(coef)"
colnames(res.cox)[3] <- "se(coef)"
colnames(res.cox)[4] <- "z"
colnames(res.cox)[5] <- "p.value"
res.cox$fdr.p.value <- p.adjust(res.cox$p.value, method = "BH") #do FDR correction
rownames(res.cox) <- colnames(cox_dat)[4:2003]
res.cox_sig <- res.cox[res.cox$p.value < 0.01, ] #filter for significant genes (0.01 of unadjusted p-value)

#identified 56 significant genes

#median of expression levels and then see if above or below (highly expressed or lowly expressed)
median <- matrix(NA,nrow(res.cox_sig),1)
median.df <- as.data.frame(median)
colnames(median.df)[1]<- "Median"
rownames(median.df)<- rownames(res.cox_sig)
for (i in 1:nrow(median.df)){
  gene <- rownames(res.cox_sig)[i]
  median.df[i,]<- median(exp.train[,gene])
}

exp.level <- matrix(NA,nrow = nrow(exp.train),ncol=nrow(res.cox_sig))
exp.level.df <- as.data.frame(exp.level)
colnames(exp.level.df) <- rownames(res.cox_sig)

for (i in 1:ncol(exp.level.df)){
  gene <- colnames(exp.level.df)[i]
  exp.level.df[,i] <- exp.train[,gene]
  exp.level.df[,i] <- ifelse(exp.level.df[,i]>median.df[i,],"high","low")
}

#cox for sig genes (high vs low)
exp.level.df$Relapse <- cox_dat$Relapse
exp.level.df$days.to.relapse <- cox_dat$days_to_relapse


num <- ncol(exp.level.df)-2
res.cox <- matrix(NA,nrow=num,ncol = 5)
for (i in 1:num){
  cox <- coxph(Surv(days.to.relapse,Relapse)~exp.level.df[,i],data = exp.level.df)
  cox_summary <- summary(cox)
  res.cox[i,]<- cox_summary$coefficients
}
res.cox <- as.data.frame(res.cox)
colnames(res.cox)[1] <- "coef"
colnames(res.cox)[2] <- "exp(coef)"
colnames(res.cox)[3] <- "se(coef)"
colnames(res.cox)[4] <- "z"
colnames(res.cox)[5] <- "p.value"
res.cox$fdr.p.value <- p.adjust(res.cox$p.value, method = "BH") #do FDR correction
rownames(res.cox) <- rownames(res.cox_sig)
#identified 2 significant genes (TREML1 and PLEKHA6) for feature selection using unadjusted p value

#Kaplan meier curve for TREML1
surv_object <- Surv(time = exp.level.df$days.to.relapse, event = exp.level.df$Relapse)
surv_object
TREML1_fit <- survfit(surv_object ~ TREML1, data = exp.level.df)
summary(TREML1_fit)
ggsurvplot(TREML1_fit, data = exp.level.df, pval = TRUE)

#Kaplan meier curve for PLEKHA6
surv_object <- Surv(time = exp.level.df$days.to.relapse, event = exp.level.df$Relapse)
surv_object
PLEKHA6_fit <- survfit(surv_object ~ PLEKHA6, data = exp.level.df)
summary(PLEKHA6_fit)
ggsurvplot(PLEKHA6_fit, data = exp.level.df, pval = TRUE)
##### however could be bias and more subjec to overfitting

#compare logistic and cox with glmnet for all genes

#glmnet for logistic
set.seed(123)
relapse <- cox_dat$Relapse
logistic <- cv.glmnet(x = exp.train,y = relapse, family = "binomial",type.measure = "class")
best_lambda <- logistic$lambda.min
best_lambda
plot(logistic) 
best_model <- glmnet(exp.train, relapse, family = "binomial", lambda = best_lambda)
selected_features <- which(coef(best_model)!=0)
logistic_features <- rownames(coef(best_model))[selected_features]
logistic_features <- logistic_features[logistic_features!= "(Intercept)"]
print(logistic_features)
# [1] "HLA-DQA2" "HLA-DQB2" "PTPRM"    "ZNF595"   "USP6"     "TIMP3"    "TBC1D3L" 
#[8] "KIR3DL2"  "GLB1L3"   "WDR49"    "ST14"     "CPAMD8"   "ADTRP"    "IL10"    
#[15] "NKG7"     "JAG1"     "TMEM52"   "MT-ATP8"  "DOCK6"    "TMCC3"    "PRKCZ"   
#[22] "FOXJ1"    "FAM240C"  "ABHD12B"  "RTN4RL2"  "CLCF1" 

surv_object <- Surv(time = exp.level.df$days.to.relapse, event = exp.level.df$Relapse)
cox_cv_model <- cv.glmnet(x = exp.train, y = surv_object, family = "cox", type.measure = "C")
best_lambda_cox <- cox_cv_model$lambda.min
best_cox_model <- glmnet(x = exp.train, y = surv_object, family = "cox", lambda = best_lambda_cox)
cox_features <- names(coef(best_cox_model)[coef(best_cox_model) != 0])
cox_features <- cox_features[cox_features != "(Intercept)"]
print(cox_features)


#nested cv with glmnet
relapse <- as.factor(as.character(relapse))
res.nestedcv <- nestcv.glmnet(y = relapse, x = exp.train,
                         family = "binomial", cv.cores = 8,
                         alphaSet = seq(0.7, 1, 0.05))
summary(res.nestedcv)

# Outer CV ROC
plot(res.nestedcv$roc, main = "Outer fold ROC", font.main = 1, col = 'blue')
legend("bottomright", legend = paste0("AUC = ", signif(pROC::auc(res.nestedcv$roc), 3)), bty = 'n')

# Inner CV ROC
rtx.inroc <- innercv_roc(res.nestedcv)
plot(rtx.inroc, main = "Inner fold ROC", font.main = 1, col = 'red')
legend("bottomright", legend = paste0("AUC = ", signif(pROC::auc(rtx.inroc), 3)), bty = 'n')


#SGL cross validation
mut.etv6
top.genes.etv6
mut.train.filt<-mut.train[,which(colSums(mut.train)>10)]
common_samples <- intersect(rownames(mut.train.filt), rownames(exp.train))
mut.exp <- cbind(exp.train[common_samples,],mut.train.filt[common_samples,])
n_genes <- ncol(exp.train)
n_mutations <- ncol(mut.train.filt)
relapse <- as.numeric(cox_dat$Relapse)
data <- list(x=mut.exp,y=relapse)
index <- c(rep(1, n_genes), rep(2, n_mutations))
cvFit <- cvSGL(data, index, type = "logit")

#top 10 mutations and 100 genes
topPreds <- order(SDs,decreasing=TRUE)[1:100]
exp.top100 <- exp.log.etv6[,topPreds]
exp.top100 <- scale(t(exp.top100))
exp.top100 <- as.matrix(exp.top100)
mut.top10<-mut.etv6[,which(colSums(mut.etv6)>127)]
mut.top100 <- as.matrix(mut.top100)
samples_in_exp_top100 <- rownames(exp.top100)
samples_in_mut_top10 <- rownames(mut.top10)
common_samples_final <- intersect(samples_in_exp_top100, samples_in_mut_top10)
exp.top100_matched_final <- exp.top100[common_samples_final, , drop = FALSE]
mut.top10_matched_final <- mut.top10[common_samples_final, , drop = FALSE]
top.exp.mut <- cbind(exp.top100_matched_final, mut.top10_matched_final)
top.exp.mut <- top.exp.mut[, !colnames(top.exp.mut) %in% "Relapse"]

relapse <- ifelse(clin.etv6$Relapse=="Relapse",1,0)
relapse <- as.numeric(relapse)
data <- list(x=top.exp.mut,y=relapse)
n_genes_sgl <- ncol(exp.top100_matched_final)
n_mutations_sgl <- ncol(mut.top10_matched_final)

data_for_cvSGL <- list(x = top.exp.mut, y = relapse)
index_for_cvSGL <- c(rep(1, n_genes_sgl), rep(2, n_mutations_sgl))
set.seed(123)

sgl_cv <- cvSGL(data = sgl_data, index = group, type = "logit", alpha = 1, nfold = 5, nlam = length(lambda.seq), lambdas = lambda.seq)
lambda.seq = exp(seq(log(1e-3), log(1), length.out = 100))
cvFit <- cvSGL(data_for_cvSGL, index_for_cvSGL, type = "logit", alpha = 1, nfold = 5, nlam = length(lambda.seq), lambdas = lambda.seq)
minlldiff <- min(cvFit$lldiff)
index_min <- which.min(cvFit$lldiff)

lambda.min <- cvFit$lambdas[index_min]


fit1 <- SGL(data=data_for_cvSGL, index=index_for_cvSGL, type = "logit", alpha = 1,lambdas = c(lambda.min))

lambda.min <- cvFit$lambda.min
lambda.1se <- sgl_cv$lambda.1se     
lambda_log_i <- lambda.min          # Fit final SGL on full training data at 1SE lambda     
sgl_fit <- SGL(data = sgl_data, index = group, type = "logit", alpha = 1, lambdas = lambda.min)

lambda_1se <- cvFit$lambda.1se
min_error_value <- cvFit$lldiff[index_min]

