library(edgeR)
library(sva)
library(caret)
library(ggplot2)
library(cluster)
library(factoextra)
library(mclust)
library(ConsensusClusterPlus)
library(dplyr)
library(pheatmap)
library(DESeq2)
library(fgsea)
library(umap)
library(discover)
library(readxl)
library(GSVA)

set.seed(123)
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

mut <- readRDS("/Users/JoshuaChung 1/Desktop/MSc Project/Part1//ETV6 subtype/Mut_SV_only.rds")
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


# Clustering --------------------------------------------------------------

set.seed(123)

#Kmeans
kmeans.etv6 <- kmeans(top.genes.etv6,centers = 13)
kmeans.config.etv6 <- kmeans.etv6$cluster
kmeans.cluster_df <- data.frame(USI= rownames(top.genes.etv6),Kmeans.cluster = as.factor(kmeans.config.etv6))
clin.etv6_1 <- merge(clin.etv6, kmeans.cluster_df, by = "USI")

#Partitioning around medoids (PAM)
pam.etv6 <- pam(top.genes.etv6,k = 20)
pam.config.etv6 <- pam.etv6$clustering
pam.cluster_df <- data.frame(USI= rownames(top.genes.etv6),PAM.cluster = as.factor(pam.config.etv6))
clin.etv6_2 <- merge(clin.etv6_1, pam.cluster_df, by = "USI")

#Hierarchical clustering
hcut.etv6 <- hcut(top.genes.etv6,method = "ward.D2",k=6)
hcut.config.etv6 <- hcut.etv6$cluster
hcut.cluster_df <- data.frame(USI= rownames(top.genes.etv6),HCUT.cluster = as.factor(hcut.config.etv6))
clin.etv6_3 <- merge(clin.etv6_2, hcut.cluster_df, by = "USI")
fviz_dend(hcut.etv6,
          k = 6, # Cut into 3 clusters (adjust this number based on your dendrogram)
          cex = 0.7, # Label size
          palette = "Set2", # Color palette for clusters
          color_labels_by_ci = TRUE, # Color labels by cluster
          rect = TRUE, # Add rectangles around clusters
          rect_border = "Set2",
          rect_fill = TRUE,
          main = "Agglomerative Hierarchical Clustering: Dendrogram with 3 Clusters"
)

#Consensus clustering
set.seed(123)
consensus_etv6 <- ConsensusClusterPlus(t(as.matrix(top.genes.etv6)),
                                       maxK = 10,
                                       reps = 50, # Number of resampling iterations
                                       pItem = 0.8, # Proportion of items to sample
                                       pFeature = 1, # Proportion of features to sample
                                       title = "ConsensusClustering",
                                       clusterAlg = "km") 

cons.config8.etv6 <- consensus_etv6[[6]]$consensusClass
cons.cluster_df <- data.frame(USI= rownames(top.genes.etv6),Cons.cluster.8 = as.factor(cons.config8.etv6),Cons.cluster.8 = as.factor(cons.config8.etv6))
clin.etv6_clust <- merge(clin.etv6_3, cons.cluster_df, by = "USI")

set.seed(123) # for reproducibility
umap_results_etv6 <- umap(top.genes.etv6)

# The UMAP coordinates are in umap_results_etv6$layout
umap_coords_df <- as.data.frame(umap_results_etv6$layout)
colnames(umap_coords_df) <- c("UMAP1", "UMAP2")

# 2. Add sample USIs (rownames of top.genes.etv6) to the UMAP coordinates
umap_coords_df$USI <- rownames(top.genes.etv6)

# 3. Add the cluster assignments (using the object you've named 'cons.config8.etv6' which now holds k=10 clusters)
umap_coords_df$Cluster <- as.factor(cons.config8.etv6[match(umap_coords_df$USI, names(cons.config8.etv6))])

# Verify that all samples from UMAP are in the cluster assignments (optional but good practice)
if(any(is.na(umap_coords_df$Cluster))) {
  warning("Some samples in UMAP results do not have corresponding cluster assignments.")
}

# 4. Plotting with ggplot2 (without sample names)
ggplot(umap_coords_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 3, alpha = 0.8) + # Adjust size and alpha as needed
  scale_color_discrete(name = "Consensus Cluster") + # Customize legend title
  labs(title = "UMAP of ETV6::RUNX1 Subtype Samples",
       subtitle = "Coloured by ConsensusClusterPlus (k=6)") + # <--- Updated subtitle for k=10
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), # Center plot title
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right")



# Mutation enrichment -----------------------------------------------------

x <- clin.etv6_clust
mut.etv6<-mut.etv6[,which(colSums(mut.etv6)>20)]
all.fisher.results.clust8 <- matrix(NA,nrow = ncol(mut.etv6),ncol = 6)
colnames(all.fisher.results.clust8)<- 1:6
rownames(all.fisher.results.clust8)<- colnames(mut.etv6)

for (c in 1:6){
  for (m in 1:ncol(mut.etv6)){
    x$Cons.cluster.8 <- ifelse(clin.etv6_clust$Cons.cluster.8 != c, "other", no = c)
    cont.tab <- cbind(Clust=x$Cons.cluster.8,Mut=mut.etv6[,m])
    cont.tab.df<- as.data.frame(cont.tab)
    cont.tab.df$Mut<- as.numeric(cont.tab.df$Mut)
    cont.tab.df<-table(cont.tab.df)
    print(cont.tab.df)
    test <-fisher.test(cont.tab.df)
    all.fisher.results.clust8[m,c] <- test$p.value
  }
}  

mutation_table <- read.csv("MSc Project/Fisher_clust_Mut.csv")
mutation_table <- as.data.frame(mutation_table)
rownames(mutation_table) <- mutation_table[,1]
mutation_table <- mutation_table[,-1]
colnames(mutation_table)[1] <- "Cluster 1"
colnames(mutation_table)[2] <- "Cluster 2"
colnames(mutation_table)[3] <- "Cluster 3"
colnames(mutation_table)[4] <- "Cluster 4"
colnames(mutation_table)[5] <- "Cluster 5"
colnames(mutation_table)[6] <- "Cluster 6"
cluster_pvalue_cols <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")
mutation_table <- mutation_table[
  rowSums(mutation_table[, cluster_pvalue_cols] < 0.05) > 0,
]
pheatmap(mutation_table,cluster_rows = FALSE,cluster_cols = FALSE,color = c("red","yellow"),breaks = c(0,0.05,1),
         main = "Gene P-values Across Clusters (Significant Genes Highlighted)",
         legend_breaks = c(0.025, 0.525),
         legend_labels =c("P < 0.05 (Red)", "P >= 0.05 (Yellow)"))
# Karyotype enrichment ----------------------------------------------------

karyotype <- readRDS("~/Desktop/Karyotype.rds")
kar<-karyotype
kar<-kar[na.omit(match(clin.etv6_clust$PatientID,rownames(kar))),] 
kar<-kar[-which(is.na(clin.etv6_clust$Karyotype)),]
kar<-kar[,which(colSums(kar)>20)]
kar.1<-kar
kar.1$PatientID <- rownames(kar)
clin.etv6_clust.kar <- merge(clin.etv6_clust,kar.1,by = "PatientID")
y<- clin.etv6_clust.kar

kar.fisher.results.clust8 <- matrix(NA,nrow = ncol(kar.1)-1,ncol = 6)
colnames(kar.fisher.results.clust8)<- 1:6
rownames(kar.fisher.results.clust8)<- colnames(kar.1[1:(ncol(kar.1)-1)])

for (c in 1:6){
  for (k in 1:ncol(kar)){
    y$Cons.cluster.8 <- ifelse(clin.etv6_clust.kar$Cons.cluster.8 != c, "other", no = c)
    cont.tab <- cbind(Clust=y$Cons.cluster.8,Kar=kar[,k])
    cont.tab.df<- as.data.frame(cont.tab)
    cont.tab.df$Kar<- as.numeric(cont.tab.df$Kar)
    cont.tab.df<-table(cont.tab.df)
    print(cont.tab.df)
    test <-fisher.test(cont.tab.df)
    kar.fisher.results.clust8[k,c] <- test$p.value
  }
}  

sig.kar.names <- rownames(kar.fisher.results.clust8)[unique(which(kar.fisher.results.clust8<0.05,arr.ind = T)[,1])]
sig.kar <- subset(kar.fisher.results.clust8,rownames(kar.fisher.results.clust8)%in%sig.kar.names)
pheatmap(sig.kar)

colors_kar <- colorRampPalette(c("#D73027", "yellow"))(length(which(sig.kar < 0.05)))
base_color_kar <- "#D73027" 

# Create a vector of colors for the heatmap matrix
heatmap_colors_kar <- matrix(base_color_kar, nrow = nrow(sig.kar), ncol = ncol(sig.kar))
heatmap_colors_kar[sig.kar < 0.05] <- colors_kar

# Define breaks to match the color scale
breaks_list_kar <- c(0, sort(unique(c(sig.kar[sig.kar < 0.05]))), 1)
pheatmap(sig.kar,
         color = c(base_color_kar, colors_kar), # Combine base color with the gradient
         breaks = c(0.05, breaks_list_kar), # Start break at 0.05
         na_col = "white", # Color for NA values
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Significant Karyotype Enrichment (p < 0.05 Highlighted)")

# DE genes etv6 ----------------------------------------------------------------
##cluster 1
all_genes <- list()
x <- clin.etv6_clust
x$Cons.cluster.8 <- ifelse(x$Cons.cluster.8 != "1", "other", "1")
x$Cons.cluster.8 <- as.factor(x$Cons.cluster.8)
exp.de.etv6<-exp_filtered[,match(x$USI,colnames(exp_filtered))]


x$Batch<-as.factor(x$Batch)
dds <- DESeqDataSetFromMatrix(countData = exp.de.etv6,
                              colData = x, 
                              design = ~ Batch + Cons.cluster.8)
dds <- DESeq(dds)
coefs<-resultsNames(dds)

res.clust1 <- results(dds,name="Cons.cluster.8_other_vs_1")
res.clust1<-as.data.frame(res.clust1)

res.clust1 <- res.clust1[order(res.clust1$padj), ]
head(res.clust1)
vol1<-EnhancedVolcano(res.clust1,
                lab = rownames(res.clust1),
                x = 'log2FoldChange',
                y = 'pvalue',
                title= "Cluster 1 DE genes")


# with(res.clust1, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot-Cluster1", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>2 and padj<0.01)
# with(subset(res.clust1, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(res.clust1, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

genes <- list(rownames(res.clust1 %>% filter(padj<.05 & abs(log2FoldChange)>1.5)))
all_genes <- append(all_genes,genes)



##Cluster2
x <- clin.etv6_clust
x$Cons.cluster.8 <- ifelse(x$Cons.cluster.8 != "2", "other", "2")
x$Cons.cluster.8 <- as.factor(x$Cons.cluster.8)
exp.de.etv6<-exp_filtered[,match(x$USI,colnames(exp_filtered))]


x$Batch<-as.factor(x$Batch)
dds <- DESeqDataSetFromMatrix(countData = exp.de.etv6,
                              colData = x, 
                              design = ~ Batch + Cons.cluster.8)
dds <- DESeq(dds)
coefs<-resultsNames(dds)

res.clust2 <- results(dds,name="Cons.cluster.8_other_vs_2")
res.clust2<-as.data.frame(res.clust2)

res.clust2 <- res.clust2[order(res.clust2$padj), ]
head(res.clust2)

# with(res.clust1, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot-Cluster2", xlim=c(-3,3)))
vol2<- EnhancedVolcano(res.clust2,
                lab = rownames(res.clust2),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Cluster 2 DE genes")
# Add colored points: blue if padj<0.01, red if log2FC>2 and padj<0.01)
# with(subset(res.clust2, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(res.clust2, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

genes <- list(rownames(res.clust2 %>% filter(padj<.05 & abs(log2FoldChange)>1.5)))
all_genes <- append(all_genes,genes)


##Cluster3
x <- clin.etv6_clust
x$Cons.cluster.8 <- ifelse(x$Cons.cluster.8 != "3", "other", "3")
x$Cons.cluster.8 <- as.factor(x$Cons.cluster.8)
exp.de.etv6<-exp_filtered[,match(x$USI,colnames(exp_filtered))]


x$Batch<-as.factor(x$Batch)
dds <- DESeqDataSetFromMatrix(countData = exp.de.etv6,
                              colData = x, 
                              design = ~ Batch + Cons.cluster.8)
dds <- DESeq(dds)
coefs<-resultsNames(dds)

res.clust3 <- results(dds,name="Cons.cluster.8_other_vs_3")
res.clust3<-as.data.frame(res.clust3)

res.clust3 <- res.clust3[order(res.clust3$padj), ]
head(res.clust3)
vol3<-EnhancedVolcano(res.clust3,
                lab = rownames(res.clust3),
                x = 'log2FoldChange',
                y = 'pvalue',
                title="Cluster 3 DE genes")

# with(res.clust1, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot-Cluster3", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>2 and padj<0.01)
# with(subset(res.clust3, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(res.clust3, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

genes <- list(rownames(res.clust3 %>% filter(padj<.05 & abs(log2FoldChange)>1.5)))
all_genes <- append(all_genes,genes)

##Cluster4
x <- clin.etv6_clust
x$Cons.cluster.8 <- ifelse(x$Cons.cluster.8 != "4", "other", "4")
x$Cons.cluster.8 <- as.factor(x$Cons.cluster.8)
exp.de.etv6<-exp_filtered[,match(x$USI,colnames(exp_filtered))]


x$Batch<-as.factor(x$Batch)
dds <- DESeqDataSetFromMatrix(countData = exp.de.etv6,
                              colData = x, 
                              design = ~ Batch + Cons.cluster.8)
dds <- DESeq(dds)
coefs<-resultsNames(dds)

res.clust4 <- results(dds,name="Cons.cluster.8_other_vs_4")
res.clust4<-as.data.frame(res.clust4)

res.clust4 <- res.clust4[order(res.clust4$padj), ]
head(res.clust4)
vol4<-EnhancedVolcano(res.clust4,
                lab = rownames(res.clust4),
                x = 'log2FoldChange',
                y = 'pvalue',
                title= "Cluster 4 DE genes")

# with(res.clust1, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot-Cluster4", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>2 and padj<0.01)
# with(subset(res.clust4, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(res.clust4, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

genes <- list(rownames(res.clust4 %>% filter(padj<.05 & abs(log2FoldChange)>1.5)))
all_genes <- append(all_genes,genes)

##Cluster5
x <- clin.etv6_clust
x$Cons.cluster.8 <- ifelse(x$Cons.cluster.8 != "5", "other", "5")
x$Cons.cluster.8 <- as.factor(x$Cons.cluster.8)
exp.de.etv6<-exp_filtered[,match(x$USI,colnames(exp_filtered))]


x$Batch<-as.factor(x$Batch)
dds <- DESeqDataSetFromMatrix(countData = exp.de.etv6,
                              colData = x, 
                              design = ~ Batch + Cons.cluster.8)
dds <- DESeq(dds)
coefs<-resultsNames(dds)

res.clust5 <- results(dds,name="Cons.cluster.8_other_vs_5")
res.clust5<-as.data.frame(res.clust5)

res.clust5 <- res.clust5[order(res.clust5$padj), ]
head(res.clust5)

vol5<-EnhancedVolcano(res.clust5,
                lab = rownames(res.clust5),
                x = 'log2FoldChange',
                y = 'pvalue',
                title="Cluster 5 DE genes")
# with(res.clust1, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot-Cluster5", xlim=c(-3,3)))
#  Add colored points: blue if padj<0.01, red if log2FC>2 and padj<0.01)
# with(subset(res.clust5, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(res.clust5, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

genes <- list(rownames(res.clust5 %>% filter(padj<.05 & abs(log2FoldChange)>1.5)))
all_genes <- append(all_genes,genes)

##Cluster6
x <- clin.etv6_clust
x$Cons.cluster.8 <- ifelse(x$Cons.cluster.8 != "6", "other", "6")
x$Cons.cluster.8 <- as.factor(x$Cons.cluster.8)
exp.de.etv6<-exp_filtered[,match(x$USI,colnames(exp_filtered))]


x$Batch<-as.factor(x$Batch)
dds <- DESeqDataSetFromMatrix(countData = exp.de.etv6,
                              colData = x, 
                              design = ~ Batch + Cons.cluster.8)
dds <- DESeq(dds)
coefs<-resultsNames(dds)

res.clust6 <- results(dds,name="Cons.cluster.8_other_vs_6")
res.clust6<-as.data.frame(res.clust6)

res.clust6 <- res.clust6[order(res.clust6$padj), ]
head(res.clust6)

vol6<-EnhancedVolcano(res.clust6,
                lab = rownames(res.clust6),
                x = 'log2FoldChange',
                y = 'pvalue',
                title="Cluster 6 DE genes")

# with(res.clust1, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot-Cluster6", xlim=c(-3,3)))
# # Add colored points: blue if padj<0.01, red if log2FC>2 and padj<0.01)
# with(subset(res.clust6, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(res.clust6, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

genes <- list(rownames(res.clust6 %>% filter(padj<.05 & abs(log2FoldChange)>1.5)))
all_genes <- append(all_genes,genes)



# Pathway enrichment ------------------------------------------------------

##subset gene symbol of each cluster to their padj
#cluster1
genes.clust1 <- all_genes[[1]]
sig.genes.clust1 <- as.data.frame(res.clust1[,"stat"])
sig.genes.clust1$Symbol <- rownames(res.clust1)
sig.genes.clust1 <- na.omit(sig.genes.clust1)
colnames(sig.genes.clust1)[1]<- "stat"
pathway.biocarta <- gmtPathways("~/Desktop/c2.cp.biocarta.v2024.1.Hs.symbols.gmt")
sig.genes.clust1 <- sig.genes.clust1[,c("Symbol","stat")]
a<- as.vector(sig.genes.clust1$stat)
names(a) <- sig.genes.clust1$Symbol
fgsea.clust1 <- fgsea(pathways = pathway.biocarta,
                      stats = a,
                      minSize = 15,
                      maxSize = 500)
fgsea.clust1.sig <-fgsea.clust1 %>% filter(padj<.05 )
topPathwaysUp <- fgsea.clust1[ES>0][head(order(pval),n=10),pathway]
topPathwaysDown <- fgsea.clust1[ES<0][head(order(pval),n=10),pathway]
topPathways <- c(topPathwaysUp,rev(topPathwaysDown))
plotGseaTable(pathway.biocarta[topPathways],a,fgsea.clust1,gseaParam=0.5)


ggplot(fgsea.clust1, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="BIOCARTA pathways NES from GSEA") + 
  theme_minimal()

##cluster2
genes.clust2 <- all_genes[[2]]
sig.genes.clust2 <- as.data.frame(res.clust2[,"stat"])
sig.genes.clust2$Symbol <- rownames(res.clust2)
sig.genes.clust2 <- na.omit(sig.genes.clust2)
colnames(sig.genes.clust2)[1]<- "stat"
pathway.biocarta <- gmtPathways("~/Desktop/c2.cp.biocarta.v2024.1.Hs.symbols.gmt")
sig.genes.clust2 <- sig.genes.clust2[,c("Symbol","stat")]
a<- as.vector(sig.genes.clust2$stat)
names(a) <- sig.genes.clust2$Symbol
fgsea.clust2 <- fgsea(pathways = pathway.biocarta,
                      stats = a,
                      minSize = 15,
                      maxSize = 500)
fgsea.clust2.sig <-fgsea.clust2 %>% filter(padj<.05 )
topPathwaysUp <- fgsea.clust2[ES>0][head(order(pval),n=10),pathway]
topPathwaysDown <- fgsea.clust2[ES<0][head(order(pval),n=10),pathway]
topPathways <- c(topPathwaysUp,rev(topPathwaysDown))
plotGseaTable(pathway.biocarta[topPathways],a,fgsea.clust2,gseaParam=0.5)


ggplot(fgsea.clust2, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="BIOCARTA pathways NES from GSEA") + 
  theme_minimal()


##cluster3
genes.clust3 <- all_genes[[3]]
sig.genes.clust3 <- as.data.frame(res.clust3[,"stat"])
sig.genes.clust3$Symbol <- rownames(res.clust3)
sig.genes.clust3 <- na.omit(sig.genes.clust3)
colnames(sig.genes.clust3)[1]<- "stat"
pathway.biocarta <- gmtPathways("~/Desktop/c2.cp.biocarta.v2024.1.Hs.symbols.gmt")
sig.genes.clust3 <- sig.genes.clust3[,c("Symbol","stat")]
a<- as.vector(sig.genes.clust3$stat)
names(a) <- sig.genes.clust3$Symbol
fgsea.clust3 <- fgsea(pathways = pathway.biocarta,
                      stats = a,
                      minSize = 15,
                      maxSize = 500)
fgsea.clust3.sig <-fgsea.clust3 %>% filter(padj<.05 )
topPathwaysUp <- fgsea.clust3[ES>0][head(order(pval),n=10),pathway]
topPathwaysDown <- fgsea.clust3[ES<0][head(order(pval),n=10),pathway]
topPathways <- c(topPathwaysUp,rev(topPathwaysDown))
plotGseaTable(pathway.biocarta[topPathways],a,fgsea.clust3,gseaParam=0.5)


ggplot(fgsea.clust3, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="BIOCARTA pathways NES from GSEA") + 
  theme_minimal()

##cluster4
genes.clust4 <- all_genes[[4]]
sig.genes.clust4 <- as.data.frame(res.clust4[,"stat"])
sig.genes.clust4$Symbol <- rownames(res.clust4)
sig.genes.clust4 <- na.omit(sig.genes.clust4)
colnames(sig.genes.clust4)[1]<- "stat"
pathway.biocarta <- gmtPathways("~/Desktop/c2.cp.biocarta.v2024.1.Hs.symbols.gmt")
sig.genes.clust4 <- sig.genes.clust4[,c("Symbol","stat")]
a<- as.vector(sig.genes.clust4$stat)
names(a) <- sig.genes.clust4$Symbol
fgsea.clust4 <- fgsea(pathways = pathway.biocarta,
                      stats = a,
                      minSize = 15,
                      maxSize = 500)
fgsea.clust4.sig <-fgsea.clust4 %>% filter(padj<.05 )
topPathwaysUp <- fgsea.clust4[ES>0][head(order(pval),n=10),pathway]
topPathwaysDown <- fgsea.clust4[ES<0][head(order(pval),n=10),pathway]
topPathways <- c(topPathwaysUp,rev(topPathwaysDown))
plotGseaTable(pathway.biocarta[topPathways],a,fgsea.clust4,gseaParam=0.5)


ggplot(fgsea.clust4, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="BIOCARTA pathways NES from GSEA") + 
  theme_minimal()

##cluster5
genes.clust5 <- all_genes[[5]]
sig.genes.clust5 <- as.data.frame(res.clust5[,"stat"])
sig.genes.clust5$Symbol <- rownames(res.clust5)
sig.genes.clust5 <- na.omit(sig.genes.clust5)
colnames(sig.genes.clust5)[1]<- "stat"
pathway.biocarta <- gmtPathways("~/Desktop/c2.cp.biocarta.v2024.1.Hs.symbols.gmt")
sig.genes.clust5 <- sig.genes.clust5[,c("Symbol","stat")]
a<- as.vector(sig.genes.clust5$stat)
names(a) <- sig.genes.clust5$Symbol
fgsea.clust5 <- fgsea(pathways = pathway.biocarta,
                      stats = a,
                      minSize = 15,
                      maxSize = 500)
fgsea.clust5.sig <-fgsea.clust5 %>% filter(padj<.05 )
topPathwaysUp <- fgsea.clust5[ES>0][head(order(pval),n=10),pathway]
topPathwaysDown <- fgsea.clust5[ES<0][head(order(pval),n=10),pathway]
topPathways <- c(topPathwaysUp,rev(topPathwaysDown))
plotGseaTable(pathway.biocarta[topPathways],a,fgsea.clust5,gseaParam=0.5)


ggplot(fgsea.clust5, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="BIOCARTA pathways NES from GSEA") + 
  theme_minimal()

##cluster6
genes.clust6 <- all_genes[[6]]
sig.genes.clust6 <- as.data.frame(res.clust6[,"stat"])
sig.genes.clust6$Symbol <- rownames(res.clust6)
sig.genes.clust6 <- na.omit(sig.genes.clust6)
colnames(sig.genes.clust6)[1]<- "stat"
pathway.biocarta <- gmtPathways("~/Desktop/c2.cp.biocarta.v2024.1.Hs.symbols.gmt")
sig.genes.clust6 <- sig.genes.clust6[,c("Symbol","stat")]
a<- as.vector(sig.genes.clust6$stat)
names(a) <- sig.genes.clust6$Symbol
fgsea.clust6 <- fgsea(pathways = pathway.biocarta,
                      stats = a,
                      minSize = 15,
                      maxSize = 500)
fgsea.clust6.sig <-fgsea.clust6 %>% filter(padj<.05 )
topPathwaysUp <- fgsea.clust6[ES>0][head(order(pval),n=10),pathway]
topPathwaysDown <- fgsea.clust6[ES<0][head(order(pval),n=10),pathway]
topPathways <- c(topPathwaysUp,rev(topPathwaysDown))
plotGseaTable(pathway.biocarta[topPathways],a,fgsea.clust6,gseaParam=0.5)


ggplot(fgsea.clust6, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="BIOCARTA pathways NES from GSEA") + 
  theme_minimal()



# Relapse Fisher test -----------------------------------------------------

relapse.fisher.results.clust8 <- matrix(NA, nrow = 1, ncol = 6) # Initialize as 1 row for p-values
colnames(relapse.fisher.results.clust8) <- paste0("Cluster_", 1:6)

for (c in 1:6) {
  x$Cons.cluster.8 <- ifelse(clin.etv6_clust$Cons.cluster.8 == c, as.character(c), "other")
  
  # Create the contingency table
  contingency_table <- table(x$Cons.cluster.8, clin.etv6_clust$Relapse)
  
  # Perform Fisher's exact test
  test <- fisher.test(contingency_table,alternative = "greater")
  
  # Store the p-value
  relapse.fisher.results.clust8[1, c] <- test$p.value
}

print(relapse.fisher.results.clust8)
#no significant difference :(

#Cluster_1 Cluster_2 Cluster_3 Cluster_4 Cluster_5 Cluster_6
#[1,] 0.4260607 0.3439954 0.4460403 0.8124316 0.6431479 0.7151727

# Within cluster gene alteration map --------------------------------------

#cluster1
cluster_samples <- clin.etv6_clust %>%
  filter(Cons.cluster.8 == 1) %>%
  pull(USI)
z <- as.data.frame(t(mut.etv6[rownames(mut.etv6) %in% cluster_samples, ]))
events <- discover.matrix(z)
subset <- rowSums(z) > 10
result.mutex <- pairwise.discover.test(events[subset,])
print(result.mutex, fdr.threshold=0.05)

old_par <- par(no.readonly = TRUE)#adjust plot parameters
par(mar = c(5.1, 10.1, 4.1, 2.1)) 
par(cex.axis = 0.7) 

plot(events[subset,])
  
#cluster2
cluster_samples <- clin.etv6_clust %>%
  filter(Cons.cluster.8 == 2) %>%
  pull(USI)
z <- as.data.frame(t(mut.etv6[rownames(mut.etv6) %in% cluster_samples, ]))
events <- discover.matrix(z)
subset <- rowSums(z) > 10
result.mutex <- pairwise.discover.test(events[subset,])
print(result.mutex, fdr.threshold=0.05)
plot(events[subset,])
  
#cluster3
cluster_samples <- clin.etv6_clust %>%
  filter(Cons.cluster.8 == 3) %>%
  pull(USI)
z <- as.data.frame(t(mut.etv6[rownames(mut.etv6) %in% cluster_samples, ]))
events <- discover.matrix(z)
subset <- rowSums(z) > 10
result.mutex <- pairwise.discover.test(events[subset,])
print(result.mutex, fdr.threshold=0.05)
plot(events[subset,])
 
#cluster4
cluster_samples <- clin.etv6_clust %>%
  filter(Cons.cluster.8 == 4) %>%
  pull(USI)
z <- as.data.frame(t(mut.etv6[rownames(mut.etv6) %in% cluster_samples, ]))
events <- discover.matrix(z)
subset <- rowSums(z) > 10
result.mutex <- pairwise.discover.test(events[subset,])
print(result.mutex, fdr.threshold=0.05)
plot(events[subset,])

#cluster5
cluster_samples <- clin.etv6_clust %>%
  filter(Cons.cluster.8 == 5) %>%
  pull(USI)
z <- as.data.frame(t(mut.etv6[rownames(mut.etv6) %in% cluster_samples, ]))
events <- discover.matrix(z)
subset <- rowSums(z) > 10
result.mutex <- pairwise.discover.test(events[subset,])
print(result.mutex, fdr.threshold=0.05)
plot(events[subset,])

#cluster6
cluster_samples <- clin.etv6_clust %>%
  filter(Cons.cluster.8 == 6) %>%
  pull(USI)
z <- as.data.frame(t(mut.etv6[rownames(mut.etv6) %in% cluster_samples, ]))
events <- discover.matrix(z)
subset <- rowSums(z) > 10
result.mutex <- pairwise.discover.test(events[subset,])
print(result.mutex, fdr.threshold=0.05)
plot(events[subset,])


par(old_par)


# t-test for presence or absence of mutation vs expression of gene --------
#across whole cohort

mut.etv6.filt<-mut.etv6[,which(colSums(mut.etv6)>10)] #filter for the mutation that is present in more than 10 samples
mut.etv6.filt[which(mut.etv6.filt>0,arr.ind = T)]<-1 #change whatever is more than 0 in the mutation matrix to 1

#extract character after "_"
alterations<-colnames(mut.etv6.filt)
genes.mut<-sub(".*_", "",alterations)

#so that all the mutations are the ones without NAs 
alterations<-alterations[-which(is.na(match(genes.mut,rownames(exp.log.etv6))))] 

alterations.t.test<-list() #creates an empty list to store the p-values from the t-test
for (i in 1:length(alterations)){
  gene<-sub(".*_", "",alterations[i]) #extract character after "_"
  exp.gene<-exp.log.etv6[gene,] #select for the rows that match the expression matrix with the gene symbols
  
  mut.gene<-mut.etv6.filt[,i] #select for the mutations
  exp.gene<-exp.gene[match(rownames(mut.etv6.filt),colnames(exp.log.etv6))] #reorders exp.gene so that it matches mut.gene
  mut.exp<- cbind(mut.gene,exp.gene) #combine mutational and expression matrix
  mut.exp<-as.data.frame(mut.exp)
  colnames(mut.exp)<-c("Mutation","Expression")
  t.test.stat<-t.test(Expression~Mutation,mut.exp)
  alterations.t.test[[alterations[[i]]]]<-t.test.stat$p.value #stores pvalue 
  #use batch corrected expression matrix (exp.log.etv6)
}

which(alterations.t.test<0.05)
#GAIN_ABI1       GAIN_ARID5B       GAIN_BCORL1          GAIN_ERG         GAIN_LCOR 
#GAIN_PTEN         GAIN_ETV6         GAIN_CKLF       GAIN_CREBBP         GAIN_TSC2 
#GAIN_UBE2I   HETDEL.5MB_ADD3  HETDEL.5MB_ARID2 HETDEL.5MB_CDKN2A  HETDEL.5MB_IKZF1
#HETDEL.5MB_RAG1   HETDEL.5MB_RAG2       HETDEL_ATRX      HETDEL_DDX3X      HETDEL_KDM6A 
#HETDEL_UBE2A     HETDEL_VPREB1     HETDEL_CDKN2A     HETDEL_CDKN2B      HETDEL_KMT2A 
#HETDEL_NR3C1     HETDEL_SLX4IP      HETDEL_PDS5B      HETDEL_INO80     HETDEL_ZNF260 


# t-test for combination of mutations vs gene expression-------------------
mut.etv6.filt
all_alterations <- list()

for (patients in 1:nrow(mut.etv6.filt)){
  mut.ind<- which(mut.etv6.filt[patients,]==1)
  names <- colnames(mut.etv6.filt)[mut.ind]
  all_alterations[[patients]] <- names
}

names(all_alterations) <- rownames(mut.etv6.filt)
unique_alteration_comb<- unique(all_alterations)






# Clustering based on gene alteration -------------------------------------
#consensus

clin.etv6_clust <- merge(clin.etv6_3, cons.cluster_df, by = "USI")

set.seed(123)
mut.etv6<-mut.etv6[,which(colSums(mut.etv6)>20)]
consensus_mut.etv6 <- ConsensusClusterPlus((as.matrix(mut.etv6)),
                                       maxK = 10,
                                       reps = 50, # Number of resampling iterations
                                       pItem = 0.8, # Proportion of items to sample
                                       pFeature = 1, # Proportion of features to sample
                                       title = "ConsensusClustering",
                                       clusterAlg = "hc")
pheatmap(t(mut.etv6),kmeans_k = 3)

ConsClust.relapse.t.test<-matrix(NA,nrow = ncol(mut.etv6),ncol = 3)
rownames(ConsClust.relapse.t.test)<- colnames(mut.etv6)
colnames(ConsClust.relapse.t.test)<- c("Cluster1","Cluster2","Cluster3")
ConsClust<- consensus_etv6[[3]]$consensusClass

for (m in 1:ncol(mut.etv6)){
  relapse.status.df <- as.data.frame(cbind(mut.etv6[,m], clin.etv6$Relapse)) #GAIN_ABI1
  colnames(relapse.status.df)[1] <- colnames(mut.etv6)[m] #GAIN_ABI1
  colnames(relapse.status.df)[2] <- "RelapseStatus"
  relapse.status.df$Cons.Clust.Asgmt <- ConsClust
  rownames(relapse.status.df)<- rownames(mut.etv6)
  relapse.status.df$Clust1 <- ifelse(relapse.status.df$Cons.Clust.Asgmt !="1","other","1")
  relapse.status.df$Clust2 <- ifelse(relapse.status.df$Cons.Clust.Asgmt !="2","other","2")
  relapse.status.df$Clust3 <- ifelse(relapse.status.df$Cons.Clust.Asgmt !="3","other","3")
  
  cont.table1<-relapse.status.df %>% filter(Cons.Clust.Asgmt=="1")
  cont.table1<-table(cont.table1[,1:2])
  test1<-fisher.test(cont.table1,alternative = "greater")
  cont.table2<-relapse.status.df %>% filter(Cons.Clust.Asgmt=="2")
  cont.table2<-table(cont.table2[,1:2])
  test2<-fisher.test(cont.table2,alternative = "greater")
  cont.table3<-relapse.status.df %>% filter(Cons.Clust.Asgmt=="3")
  cont.table3<-table(cont.table3[,1:2])
  test3<-fisher.test(cont.table3,alternative = "greater")
  
  ConsClust.relapse.t.test[m,1]<-test1$p.value#GAIN_ABI1
  ConsClust.relapse.t.test[m,2]<-test2$p.value#GAIN_ABI1
  ConsClust.relapse.t.test[m,3]<-test3$p.value#GAIN_ABI1
}

clus.sig<-list()
clus.sig[[1]]<-rownames(ConsClust.relapse.t.test)[which(ConsClust.relapse.t.test[,1]<0.05)]
clus.sig[[2]]<-rownames(ConsClust.relapse.t.test)[which(ConsClust.relapse.t.test[,2]<0.05)]

mut.etv6.clus<-mut.etv6[,clus.sig[[1]]]
relapse.status.df <- as.data.frame(cbind( ConsClust, clin.etv6$Relapse)) 
relapse.status.df <- relapse.status.df %>% filter(ConsClust=="1")
colnames(relapse.status.df)[2]<-"Relapse"

mut.etv6.clus<-mut.etv6.clus[rownames(relapse.status.df),]
z <- as.data.frame(t(mut.etv6.clus))
events <- discover.matrix(z)
subset <- rowSums(z) > 5
result.mutex <- pairwise.discover.test(events[subset,])
plot(events[subset,])
#heatmap showing relapse in cluster1
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
annotation_df <- data.frame(Relapse = relapse.status.df$Relapse)
rownames(annotation_df) <- rownames(relapse.status.df)

annotation_df$Relapse <- factor(annotation_df$Relapse, levels = c("Relapse", "None"))

ann_colors <- list(
  Relapse = c("Relapse" = "red", "None" = "yellow")
)

sorted_annotation_df <- annotation_df[order(annotation_df$Relapse), , drop = FALSE]

ordered_sample_names <- rownames(sorted_annotation_df)

heatmap_matrix_ordered_cols <- t(mut.etv6.clus)[ , ordered_sample_names, drop = FALSE]

pheatmap(heatmap_matrix_ordered_cols,
         annotation_col = sorted_annotation_df,
         annotation_colors = ann_colors,
         color = c("lightgray", "black"),
         legend = FALSE,
         show_colnames = FALSE,
         main = "Significant mutation mutations - Cluster1")

#Heatmap showing relapse in cluster2
mut.etv6.clus<-mut.etv6[,clus.sig[[2]]]
relapse.status.df <- as.data.frame(cbind( ConsClust, clin.etv6$Relapse)) 
relapse.status.df <- relapse.status.df %>% filter(ConsClust=="2")
colnames(relapse.status.df)[2]<-"Relapse"

mut.etv6.clus<-mut.etv6.clus[rownames(relapse.status.df),]
z <- as.data.frame(t(mut.etv6.clus))
events <- discover.matrix(z)
subset <- rowSums(z) > 5
result.mutex <- pairwise.discover.test(events[subset,])
plot(events[subset,])
#heatmap showing relapse in cluster1
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
annotation_df <- data.frame(Relapse = relapse.status.df$Relapse)
rownames(annotation_df) <- rownames(relapse.status.df)

annotation_df$Relapse <- factor(annotation_df$Relapse, levels = c("Relapse", "None"))

ann_colors <- list(
  Relapse = c("Relapse" = "red", "None" = "yellow")
)

sorted_annotation_df <- annotation_df[order(annotation_df$Relapse), , drop = FALSE]

ordered_sample_names <- rownames(sorted_annotation_df)

heatmap_matrix_ordered_cols <- t(mut.etv6.clus)[ , ordered_sample_names, drop = FALSE]

pheatmap(heatmap_matrix_ordered_cols,
         annotation_col = sorted_annotation_df,
         annotation_colors = ann_colors,
         color = c("lightgray", "black"),
         legend = FALSE,
         show_colnames = FALSE,
         main = "Significant mutation mutations - Cluster2")


#biclustering
library(biclust)
library(sparseBC)
library(BiBitR)

set.seed(123)
mut.matrix<- as.matrix(t(mut.etv6))
mut.matrix[mut.matrix>1]<- 1
mut.biclust_result <- biclust(mut.matrix,method = BCBimax(),minr= 3, minc= 5)
print(mut.biclust_result)

heatmapBC(mut.matrix, bicResult = mut.biclust_result, main= "Biclustering heatmap")



num_biclusters <- mut.biclust_result@Number
message(paste("Total number of biclusters found:", num_biclusters))

# Extract Mutation (Row) Memberships for each Bicluster ---

# Rows are original matrix rows (your mutations)
# Columns are the identified biclusters
mutation_membership_matrix <- mut.biclust_result@RowxNumber

bicluster_mutations <- list()

if (num_biclusters > 0) {
  for (i in 1:num_biclusters) {
    # Get the logical vector for the current bicluster
    is_mutation_in_bicluster <- mutation_membership_matrix[, i]
    
    # Get the names of the mutations (rows) that are TRUE for this bicluster
    mutations_in_this_bicluster <- rownames(mut.matrix)[is_mutation_in_bicluster]
    
    # Store them in the list
    bicluster_mutations[[paste0("Bicluster_", i, "_Mutations")]] <- mutations_in_this_bicluster
    
    message(paste0("Bicluster ", i, " has ", length(mutations_in_this_bicluster), " mutations."))
  }
} 

# --- 3. Extract Patient (Column) Memberships for each Bicluster ---

# mut.biclust_result@NumberxCol is a logical matrix:
# Rows are the identified biclusters
# Columns are original matrix columns (your patients)
patient_membership_matrix <- mut.biclust_result@NumberxCol

# Create a list to store the names of patients in each bicluster


if (num_biclusters > 0) {
  for (i in 1:num_biclusters) {
    # Get the logical vector for the current bicluster (row i in this matrix)
    is_patient_in_bicluster <- patient_membership_matrix[i, ]
    
    # Get the names of the patients (columns) that are TRUE for this bicluster
    patients_in_this_bicluster <- colnames(mut.matrix)[is_patient_in_bicluster]
    
    # Store them in the list
    bicluster_patients[[paste0("Bicluster_", i, "_Patients")]] <- patients_in_this_bicluster
    
    message(paste0("Bicluster ", i, " has ", length(patients_in_this_bicluster), " patients."))
    # Optionally print first few patient names
    # print(head(patients_in_this_bicluster))
  }
} 

mut.relapse.t.test <- list()
for (c in 1:100){
clust_patient_ids <- bicluster_patients[[c]]
relapse.status.df <- as.data.frame(cbind(clin.etv6$USI, clin.etv6$Relapse))
colnames(relapse.status.df)[1] <- "ID"
colnames(relapse.status.df)[2] <- "RelapseStatus"
relapse.status.df$ID <- as.character(relapse.status.df$ID)
relapse.status.df$Bicluster_Membership <- ifelse(relapse.status.df$ID %in% clust_patient_ids, "1", "0")

contingency_table <- table(relapse.status.df$Bicluster_Membership, relapse.status.df$RelapseStatus)

test <- fisher.test(contingency_table,alternative = "greater")
mut.relapse.t.test[[c]] <- test$p.value
}


# Differentiation state ---------------------------------------------------

HSC1_genes<- read_excel("~/Desktop/MSc Project/HSC1.xlsx")
HSC2_genes<- read_excel("~/Desktop/MSc Project/HSC2.xlsx")
CLP_genes<- read_excel("~/Desktop/MSc Project/CLP.xlsx")
proB_G2_M_S_genes<- read_excel("~/Desktop/MSc Project/pro-B G2_M_S.xlsx")
proB_G1_genes<- read_excel("~/Desktop/MSc Project/pro-B G1.xlsx")
preB_G2_M_S_genes<- read_excel("~/Desktop/MSc Project/pre-B I G1.xlsx")
preB_1_G1_genes<- read_excel("~/Desktop/MSc Project/HSC2.xlsx")
preB_2_G1_genes<- read_excel("~/Desktop/MSc Project/pre-B II G1.xlsx")
immature_B_1_genes<- read_excel("~/Desktop/MSc Project/immature B 1.xlsx")
immature_B_2_genes<- read_excel("~/Desktop/MSc Project/immature B 2.xlsx")
immature_B_3_genes<- read_excel("~/Desktop/MSc Project/immature B 3.xlsx")

n_genes<-100

HSC1_genes <- HSC1_genes %>%
  arrange(desc(score)) %>% 
  slice_head(n = n_genes) 
HSC2_genes <- HSC2_genes %>%
  arrange(desc(score)) %>% 
  slice_head(n = n_genes)
CLP_genes <- CLP_genes %>%
  arrange(desc(score)) %>% 
  slice_head(n = n_genes)
proB_G2_M_S_genes <- proB_G2_M_S_genes %>%
  arrange(desc(score)) %>% 
  slice_head(n = n_genes)
proB_G1_genes <- proB_G1_genes %>%
  arrange(desc(score)) %>% 
  slice_head(n = n_genes)
preB_G2_M_S_genes <- preB_G2_M_S_genes %>%
  arrange(desc(score)) %>% 
  slice_head(n = n_genes)
preB_1_G1_genes <- preB_1_G1_genes %>%
  arrange(desc(score)) %>% 
  slice_head(n = n_genes)
preB_2_G1_genes <- preB_2_G1_genes %>%
  arrange(desc(score)) %>% 
  slice_head(n = n_genes)
immature_B_1_genes <- immature_B_1_genes %>%
  arrange(desc(score)) %>% 
  slice_head(n = n_genes)
immature_B_2_genes <- immature_B_2_genes %>%
  arrange(desc(score)) %>% 
  slice_head(n = n_genes)
immature_B_3_genes <- immature_B_3_genes %>%
  arrange(desc(score)) %>% 
  slice_head(n = n_genes)

comb_gene_set<- list()
comb_gene_set[["HSC1 genes"]]<- as.list(HSC1_genes)$gene
comb_gene_set[["HSC2 genes"]]<- as.list(HSC2_genes)$gene
comb_gene_set[["CLP genes"]]<- as.list(CLP_genes)$gene
comb_gene_set[["pro-B G2/M/S genes"]]<- as.list(proB_G2_M_S_genes)$gene
comb_gene_set[["pro-B G1 genes"]]<- as.list(proB_G1_genes)$gene
comb_gene_set[["pre-B G1/M/S genes"]]<- as.list(preB_G2_M_S_genes)$gene
comb_gene_set[["pre-B I G1 genes"]]<- as.list(preB_1_G1_genes)$gene
comb_gene_set[["pre-B II G1"]]<- as.list(preB_2_G1_genes)$gene
comb_gene_set[["immature B 1 genes"]]<- as.list(immature_B_1_genes)$gene
comb_gene_set[["immature B 2 genes"]]<- as.list(immature_B_2_genes)$gene
comb_gene_set[["immature B 3 genes"]]<- as.list(immature_B_3_genes)$gene


expression_matrix_for_gsva <- as.matrix(exp.log.etv6)
expression_matrix_for_gsva <- expression_matrix_for_gsva[!duplicated(rownames(expression_matrix_for_gsva)), ]
gsvaPar <- gsvaParam(expression_matrix_for_gsva, comb_gene_set)
gsvaPar
gsva.es <- gsva(gsvaPar, verbose=FALSE)
dim(gsva.es)

cell_states_to_test<- c("HSC1 genes","HSC2 genes","CLP genes","pro-B G2/M/S genes","pro-B G1 genes","pre-B G1/M/S genes","pre-B I G1 genes","pre-B II G1","immature B 1 genes","immature B 2 genes","immature B 3 genes")
#Cluster vs others
cluster_state_t.test <- list()
cluster_gsva <- t(gsva.es)
cluster_gsva <- as.data.frame(cluster_gsva)
cluster_gsva$ClustAssignt <- clin.etv6_clust$Cons.cluster.8

for (k in 1:6){
  for (i in cell_states_to_test){
    cluster_state_gsva <- cluster_gsva
    cluster_state_gsva$ClustAssignt <- ifelse(cluster_state_gsva$ClustAssignt != k, "other",no = k)
    state <- cluster_state_gsva[,i]
    clust <- cluster_state_gsva$ClustAssignt
    state.clust <- cbind(state,clust)
    state.clust <- as.data.frame(state.clust)
    colnames(state.clust)<- c("State","Clust")
    state.clust$State<- as.numeric(state.clust$State)
    state.clust$Clust<- as.factor(state.clust$Clust)
    res <- t.test(State~Clust,data = state.clust)
    cluster_state_t.test[[paste0("Cluster_", k, "_", i)]]<- res$p.value
  }
}

#Cluster 3 bar plot showing average GSVA score for enriched cell type
cluster3 <- cluster_gsva[cluster_gsva$ClustAssignt == 3, ]
CLP_3 <- mean(cluster3$`CLP genes`)
all_CLP <- mean(cluster_gsva$`CLP genes`)
CLP_dat <- data.frame(ClusterAssignment= c("Cluster 3","Other Clusters"),MeanGSVA=c(-0.05202981,0.006313887))

ggplot(CLP_dat,aes(x=ClusterAssignment,y = MeanGSVA))+geom_bar(stat="identity")

cluster_state_gsva <- cluster_gsva
cluster_state_gsva$ClustAssignt <- ifelse(cluster_state_gsva$ClustAssignt != "3", "other",no = "3")
CLP <- cbind(cluster_state_gsva$`CLP genes`,cluster_state_gsva$ClustAssignt)
CLP <- as.data.frame(CLP)
colnames(CLP)[2] <- "Cluster"
colnames(CLP)[1] <- "GSVA score"
CLP$`GSVA score` <- as.numeric(CLP$`GSVA score`)
CLP$Cluster <- as.factor(CLP$Cluster)
ggplot(CLP,aes(x=Cluster,y=`GSVA score`))+
  geom_boxplot()

#Cluster 5
cluster_state_gsva <- cluster_gsva
cluster_state_gsva$ClustAssignt <- ifelse(cluster_state_gsva$ClustAssignt != "5", "other",no = "5")
preB2 <- cbind(cluster_state_gsva$`pre-B II G1`,cluster_state_gsva$ClustAssignt)
preB2 <- as.data.frame(preB2)
colnames(preB2)[2] <- "Cluster"
colnames(preB2)[1] <- "GSVA score"
preB2$`GSVA score` <- as.numeric(preB2$`GSVA score`)
preB2$Cluster <- as.factor(preB2$Cluster)
ggplot(preB2,aes(x=Cluster,y=`GSVA score`))+
  geom_boxplot()
cluster5 <- cluster_gsva[cluster_gsva$ClustAssignt == 5, ]
PreB_5 <- mean(cluster5$`pre-B II G1`)
all_Preb <- mean(cluster_gsva$`pre-B II G1`)
#Cluster 3 was significantly enriched for CLP, Cluster 5 for Pre B(G1)

gsva_cluster <- matrix(NA,nrow = 6,ncol = 12)
gsva_cluster <- as.data.frame(gsva_cluster)
colnames(gsva_cluster) <- colnames(cluster_state_gsva)
gsva_cluster <- gsva_cluster[,-ncol(gsva_cluster)]
for (i in 1:nrow(gsva_cluster)){
rownames(gsva_cluster)[i]<- paste0("Cluster", i)
}

#HSC1
for (i in 1:nrow(gsva_cluster)){
  cluster_state_gsva <- cluster_gsva
  cluster_state_gsva$ClustAssignt <- ifelse(cluster_state_gsva$ClustAssignt != i, "other",no = i)
  name <- colnames(cluster_state_gsva)[i]
  mean <- mean(cluster_state_gsva$name)
  gsva_cluster[,i] <- mean
}


gsva_cluster <- matrix(NA, nrow = 6, ncol = length(cell_states_to_test))
gsva_cluster <- as.data.frame(gsva_cluster)


colnames(gsva_cluster) <- cell_states_to_test

for (k_row_idx in 1:nrow(gsva_cluster)){ 
  rownames(gsva_cluster)[k_row_idx] <- paste0("Cluster", k_row_idx)
}

for (k_cluster_id in 1:6){ 
  
  data_for_current_cluster <- cluster_gsva[cluster_gsva$ClustAssignt == k_cluster_id, ]
  
  for (cell_state_name in cell_states_to_test){ 
   
    current_mean_gsva <- mean(data_for_current_cluster[[cell_state_name]])
    
   
    col_idx <- which(colnames(gsva_cluster) == cell_state_name)
    gsva_cluster[k_cluster_id, col_idx] <- current_mean_gsva
  }
}

print(gsva_cluster)





# Progression/Regression Site  --------------------------------------------
#subset relapsed patients then see if bone marrow or others then see how it relates to the enrichment scores

relapsed_patients_clin <- clin.etv6_clust %>% filter(Relapse=="Relapse")
relapsed_patients_clin <- relapsed_patients_clin %>%filter(progression_or_recurrence_anatomic_site!="'--")
bm_recurrence_patients_clin<- relapsed_patients_clin
bm_recurrence_patients_clin$progression_or_recurrence_anatomic_site <- ifelse(bm_recurrence_patients_clin$progression_or_recurrence_anatomic_site !="Bone marrow","other","Bone marrow")
relapsed_patients_gsva<- gsva.es[,match(relapsed_patients_clin$USI,colnames(gsva.es))]

#t-test
clin_t.test<- relapsed_patients_gsva
gsva.df <- as.data.frame(t(relapsed_patients_gsva))
gsva.df$Progression_site <- bm_recurrence_patients_clin$progression_or_recurrence_anatomic_site

recurrence_site_t.test <- list()
for (i in cell_states_to_test) {
  diff.state<- colnames(gsva.df)[i]
  state<-gsva.df[,i]
  site<- gsva.df$Progression_site
  state.site<- cbind(state,site)
  state.site.df<- as.data.frame(state.site)
  colnames(state.site.df)<- c("State","Site")
  state.site.df$State<- as.numeric(state.site.df$State)
  state.site.df$Site<- as.factor(state.site.df$Site)
  res <- t.test(State~Site,data = state.site.df)
  recurrence_site_t.test[[i]]<- res$p.value
}
  
#HSC2 and Pre-B genes are enriched for relapse in bone marrow
site<- gsva.df
HSC2<- cbind(site$`HSC2 genes`,site$Progression_site)
HSC2 <- as.data.frame(HSC2)
colnames(HSC2)[2] <- "Progression site"
colnames(HSC2)[1] <- "GSVA score for HSC2 genes"
HSC2$`GSVA score for HSC2 genes` <- as.numeric(HSC2$`GSVA score for HSC2 genes`)
HSC2$`Progression site` <- as.factor(HSC2$`Progression site`)
ggplot(HSC2,aes(x=`Progression site`,y=`GSVA score for HSC2 genes`))+
  geom_boxplot()

site<- gsva.df
preB<- cbind(site$`pre-B I G1 genes`,site$Progression_site)
preB <- as.data.frame(preB)
colnames(preB)[2] <- "Progression site"
colnames(preB)[1] <- "GSVA score for Pre-B genes"
preB$`GSVA score for Pre-B genes` <- as.numeric(preB$`GSVA score for Pre-B genes`)
preB$`Progression site` <- as.factor(preB$`Progression site`)
ggplot(preB,aes(x=`Progression site`,y=`GSVA score for Pre-B genes`))+
  geom_boxplot()









#  RandomForest -----------------------------------------------------------
set.seed(123)
diffstate.split<- initial_split(gsva.df,prop=0.8,strata= "Progression_site")
diffstate.train<- training(diffstate.split)
diffstate.test<- testing(diffstate.split)

colnames(diffstate.train)<-c("HSC1","HSC2","CLP","proB_G2_M_S","proB_G1","preB_G1_M_S","preBI_G1","preBII_G1","immatureB_1","immatureB_2","immatureB_3","Progression_site")
diffstate.train$Progression_site <- as.factor(diffstate.train$Progression_site)
diffstate.rf <- randomForest(Progression_site~., data = diffstate.train, proximity = TRUE, mtry = sqrt(11), ntree = 500)
plot(diffstate.rf)

varImpPlot(diffstate.rf,
           sort=T,
           n.var = 11,
           main="Top 5-Variable Importance")

reprtree:::plot.getTree(diffstate.rf)
#no significance

SE_celltype <- load("~/Desktop/SE_celltype.rda")


# -------------------------------------------------------------------------


#Top 5 significant pathways
# --- Prepare data for plotting (select top 5 pathways per cluster) ---
top_n_pathways_per_cluster <- 5 # Changed from 10 to 5

if (nrow(all_pathways_df) == 0) {
       message("No significant pathways found at padj < 0.05 across any cluster.")
       message("You may need to adjust the 'padj' threshold or inspect your data.")
    } else {
         # Get the top N pathways for each cluster
            plot_data_prep <- all_pathways_df %>%
                 group_by(Cluster, Direction) %>%
                 arrange(desc(abs(NES))) %>%
                 slice_head(n = top_n_pathways_per_cluster) %>%
                 ungroup() %>%
                 # Ensure Cluster factor levels for correct loop order
                 mutate(Cluster = factor(Cluster, levels = paste("Cluster", 1:6)))
           
             # --- Generate Separate Plots for Each Cluster ---
             unique_clusters <- levels(plot_data_prep$Cluster)
             
               for (current_cluster in unique_clusters) {
                     cluster_plot_data <- plot_data_prep %>%
                           filter(Cluster == current_cluster) %>%
                           mutate(pathway = factor(pathway, levels = unique(pathway[order(NES, decreasing = FALSE)])))
                     
                       p_single_cluster <- ggplot(cluster_plot_data, aes(x = NES, y = pathway, fill = Direction)) +
                             geom_bar(stat = "identity") +
                             scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
                             labs(
                                   title = paste0("Top ", top_n_pathways_per_cluster, " Significant Pathways for ", current_cluster),
                                   x = "Normalized Enrichment Score (NES)",
                                   y = "Pathway"
                               ) +
                             theme_minimal(base_size = 14) +
                             theme(
                                   plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                                   axis.text.y = element_text(size = 12, color = "black", hjust = 0),
                                   axis.title = element_text(face = "bold"),
                                   # --- NEW: Move Y-axis title further away ---
                                   axis.title.y = element_text(margin = margin(r = 20)), # Adjust 'r' (right) value as needed
                                   legend.position = "bottom",
                                   plot.margin = unit(c(0.5, 0.5, 0.5, 2), "cm") # Left margin set to 2cm, adjust if needed
                               )
                       
                         print(p_single_cluster)
                       
                         # Optional: Saving the plots with sufficient width
                         # file_name <- paste0("top", top_n_pathways_per_cluster, "_pathways_", gsub(" ", "_", tolower(current_cluster)), ".png")
                         # ggsave(file.path(base_path, file_name), plot = p_single_cluster, width = 10, height = 6, dpi = 300)
                     }
         }

relapse_table <- read_excel("cluster.xlsx")


# Assuming relapse_table is your data frame after reading
relapse_table <- read_excel("cluster.xlsx")
  ggplot(relapse_table, aes(fill=`Relapse status`, y=number, x=Cluster)) + 
  geom_bar(position="fill", stat="identity")+
    labs(y = "Proportion of Patients")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x = element_blank())



# Now set the new unique column as row names
rownames(relapse_table) <- relapse_table$unique_id

# You might then remove the original first column if it's no longer needed as a regular column
# relapse_table <- relapse_table[,-1]


# See if mutation correlates with expression level ------------------------

gen_alt <- mut.etv6
gen_alt <- gen_alt[,which(colSums(gen_alt)>10)]

exp <- t(exp.log.etv6)






