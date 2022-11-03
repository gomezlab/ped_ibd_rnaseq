# %%
## Ben's pipeline for RNA-seq instructions 
library("AnnotationDbi")
library("org.Hs.eg.db")
library(RUVSeq)
source("downstreamAnalysis_RNAseqFunctions.R")

# %%
# Read in data
cts = as.matrix(read.csv("../../data/processed_data/data_c_raw.csv", sep=',', row.names='gene'))

coldata = read.csv('../../data/metadata/coldata_cd_c.csv', sep=',', row.names='sample')
dim(coldata)
dim(cts)
#get coldata row names
coldata$sample = rownames(coldata)
sample = rownames(coldata)

head(coldata)


# %%
#set columns as factor
factor_cols = c('batch', 'sex', 'B2', 'B3', 'rem', 'surg')
for (col in factor_cols) {
  coldata[,col] = as.factor(coldata[,col])
}
batch = coldata$batch
sexBin = coldata$sex
age = coldata$dx_age

B2 = coldata$B2
B3 = coldata$B3
surg = coldata$surg
rem = coldata$rem

tin = coldata$TIN
#round counts to integers
cts <- round(cts);

#check if row names of coldata match colnames of cts
all(rownames(coldata) == colnames(cts))
#adjust for batch and sexBin
design <- as.formula(~ batch + sex)
#create the dds object
dds = DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = design)

vsd = vst(dds)
mat= assay(vsd)

# %%
pcaData <- plotPCA(vsd, intgroup=c("B2"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#save the plot
png("results/c/B2/pca_B2_tin.png", width=1024, height=1024)
ggplot(pcaData, aes(PC1, PC2, color=tin, shape=B2)) +
    geom_point(size=5) +
    #geom_text(aes(label=sample),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    theme(text = element_text(size = 24))
dev.off()

# %%
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, batch2 = sexBin, covariates=cbind(tin))

# %%
pcaData <- plotPCA(vsd, intgroup=c("B2"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#save the plot
png("results/c/B2/pca_B2_tin_corrected.png", width=1024, height=1024)
ggplot(pcaData, aes(PC1, PC2, color=tin, shape=B2)) +
    geom_point(size=5) +
    #geom_text(aes(label=sample),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    theme(text = element_text(size = 24))
dev.off()

# %%
#RUV normalization

# Calculate the average gene expression of each genes and take the top 5000 highest expressed
nc <- as.data.frame(assay(vsd))
nc$avg_exp <- rowSums(nc,na.rm=TRUE) / ncol(nc) # x should be your normalized matrix
nc <- arrange(nc, -avg_exp)%>%dplyr::select(-avg_exp)
nc <- nc[1:5000,]
# Calculate the variance of each genes, and choose the lowest 1000 genes as the negative control gene
nc$row_stv <- rowSds(as.matrix(nc))/(rowSums(nc,na.rm=TRUE) / ncol(nc))
nc <- arrange(nc, row_stv) %>% dplyr::select(-row_stv)
nc <- nc[1:1000,]
##The newSeqExpressionSet is given the raw count matrix, phenotypic data, and then the names of the samples 
# Create a new Expression Set and perform an upper quartile normalization
nc <- round(nc*1000)
set <- newSeqExpressionSet(as.matrix(nc),phenoData = data.frame(colData(dds),row.names=colnames(nc)))
#set <- betweenLaneNormalization(set, which="upper")
##Run RUV
##k refers to the number of "unwanted factors of variation" (kinda like PCs in PCA), with rows referring to the samples 
##and columns to the factors. You can play around with k, but 5 is a good starting number.
set_ruv<- RUVg(set, rownames(nc), k=1)

# Save the phenotype data
ruv_data <- pData(set_ruv)

# %%
png("results/c/B2/RUV_rle_plot_1.png")
plotRLE(set_ruv, outline=FALSE, ylim=c(-0.1, 0.1), col=as.numeric(as.factor(coldata$B2)))
dev.off()

# %%
png('results/c/B2/RUV_pca_plot_1.png')
plotPCA(set, col=as.numeric(as.factor(coldata$B2)), k=2, cex=1.2)
dev.off()

# %%
head(ruv_data)

# %%
write.csv(ruv_data, file = "results/c/ruv_data_cd.csv")

# %%
ruv_data$batch = as.factor(ruv_data$batch)

batch = ruv_data$batch
sexBin = as.numeric(as.factor(ruv_data$sex))
age = ruv_data$dx_age
ruv_data$B3 <- as.factor(ruv_data$B3)
ruv_data$B2 <- as.factor(ruv_data$B2)
ruv_data$rem <- as.factor(ruv_data$rem)
ruv_data$surg <- as.factor(ruv_data$surg)
B2 = ruv_data$B2
B3 = ruv_data$B3
surg = ruv_data$surg
rem = ruv_data$rem
tin = ruv_data$TIN

# %%
w_1 = ruv_data$W_1
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + TIN + W_1 + B2)

vsd = vst(dds)
mat= assay(vsd)
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, batch2=sexBin, covariates=cbind(tin + w_1))

pcaData <- plotPCA(vsd, intgroup=c("B2"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#save the plot
png("results/c/B2/pca_B2_tin_w1.png", width=1024, height=1024)
ggplot(pcaData, aes(PC1, PC2, color=B2)) +
    geom_point(size=5) +
    #geom_text(aes(label=sample),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    scale_color_manual(values = c("0" = "#0066cc", "1" = "#990000")) + 
    theme(text = element_text(size = 24))
dev.off()

# %%
deg_counts <- data.frame(matrix(ncol=2, nrow=0))
x <- c('complication', 'n_DEGs')
colnames(deg_counts) <- x
write.csv(deg_counts, file = "results/c/deg_counts_cd.csv", row.names = FALSE)

# %%
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + TIN + W_1 + B2)
y <- DGEList(counts=counts(dds1), group=coldata$B2)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]

dds1 <- DESeq(dds1)
res <- results(dds1)
resOrdered <- res[order(res$padj),]
#create a column named gene which is the gene name
resOrdered$gene = rownames(resOrdered)
#set gene as the first column by making the last column the first column
resOrdered = resOrdered[,c(7, 1, 2, 3, 4, 5, 6)]

deg_counts[nrow(deg_counts)+1,] <- c('B2', sum(resOrdered$padj < 0.05))

#plot a histogram of the p-values
png("results/c/B2/pvals_hist.png", width=512, height=512)
hist(resOrdered$pvalue, breaks=20, col="grey", main="p-value distribution", xlab="p-value")
dev.off()

write.csv(resOrdered, file = "results/c/B2/results_deseq.csv", sep = ",", row.names = FALSE, col.names = TRUE)

# %%
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + TIN + W_1 + B3)
y <- DGEList(counts=counts(dds1), group=coldata$B3)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]

dds1 <- DESeq(dds1)
res <- results(dds1)
resOrdered <- res[order(res$padj),]
#create a column named gene which is the gene name
resOrdered$gene = rownames(resOrdered)
#set gene as the first column by making the last column the first column
resOrdered = resOrdered[,c(7, 1, 2, 3, 4, 5, 6)]

deg_counts[nrow(deg_counts)+1,] <- c('B3', sum(resOrdered$padj < 0.05))

#plot a histogram of the p-values
png("results/c/B3/pvals_hist.png", width=512, height=512)
hist(resOrdered$pvalue, breaks=20, col="grey", main="p-value distribution", xlab="p-value")
dev.off()

write.csv(resOrdered, file = "results/c/B3/results_deseq.csv", sep = ",", row.names = FALSE, col.names = TRUE)

# %%
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + TIN + W_1 + rem)
y <- DGEList(counts=counts(dds1), group=coldata$rem)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]

dds1 <- DESeq(dds1)
res <- results(dds1)
resOrdered <- res[order(res$padj),]
#create a column named gene which is the gene name
resOrdered$gene = rownames(resOrdered)
#set gene as the first column by making the last column the first column
resOrdered = resOrdered[,c(7, 1, 2, 3, 4, 5, 6)]

deg_counts[nrow(deg_counts)+1,] <- c('rem', sum(resOrdered$padj < 0.05))

#plot a histogram of the p-values
png("results/c/rem/pvals_hist.png", width=512, height=512)
hist(resOrdered$pvalue, breaks=20, col="grey", main="p-value distribution", xlab="p-value")
dev.off()

write.csv(resOrdered, file = "results/c/rem/results_deseq.csv", sep = ",", row.names = FALSE, col.names = TRUE)

# %%
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + TIN + W_1 + surg)
y <- DGEList(counts=counts(dds1), group=coldata$surg)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]

dds1 <- DESeq(dds1)
res <- results(dds1)
resOrdered <- res[order(res$padj),]
#create a column named gene which is the gene name
resOrdered$gene = rownames(resOrdered)
#set gene as the first column by making the last column the first column
resOrdered = resOrdered[,c(7, 1, 2, 3, 4, 5, 6)]

deg_counts[nrow(deg_counts)+1,] <- c('surg', sum(resOrdered$padj < 0.05))

#plot a histogram of the p-values
png("results/c/surg/pvals_hist.png", width=512, height=512)
hist(resOrdered$pvalue, breaks=20, col="grey", main="p-value distribution", xlab="p-value")
dev.off()

write.csv(resOrdered, file = "results/c/surg/results_deseq.csv", sep = ",", row.names = FALSE, col.names = TRUE)

# %%
head(deg_counts)

# %%
write.csv(deg_counts, file = "results/c/deg_counts.csv", row.names = FALSE)

# %%
# Read in data
cts = as.matrix(read.csv("../data/processed_data/data_i_raw.csv", sep=',', row.names='gene'))

coldata = read.csv('../data/metadata/coldata_cd_i.csv', sep=',', row.names='sample')
dim(coldata)
dim(cts)
#get coldata row names
coldata$sample = rownames(coldata)
sample = rownames(coldata)

head(coldata)

#set columns as factor
factor_cols = c('batch', 'sex', 'B2', 'B3', 'rem', 'surg')
for (col in factor_cols) {
  coldata[,col] = as.factor(coldata[,col])
}
batch = coldata$batch
sexBin = coldata$sex
age = coldata$dx_age

B2 = coldata$B2
B3 = coldata$B3
surg = coldata$surg
rem = coldata$rem

tin = coldata$TIN
#round counts to integers
cts <- round(cts);

#check if row names of coldata match colnames of cts
all(rownames(coldata) == colnames(cts))
#adjust for batch and sexBin
design <- as.formula(~ batch + sex)
#create the dds object
dds = DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = design)

vsd = vst(dds)
mat= assay(vsd)
pcaData <- plotPCA(vsd, intgroup=c("B2"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#save the plot
png("results/i/B2/pca_B2_tin.png", width=1024, height=1024)
ggplot(pcaData, aes(PC1, PC2, color=tin, shape=B2)) +
    geom_point(size=5) +
    #geom_text(aes(label=sample),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    theme(text = element_text(size = 24))
dev.off()
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, batch2 = sexBin, covariates=cbind(tin))
pcaData <- plotPCA(vsd, intgroup=c("B2"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#save the plot
png("results/i/B2/pca_B2_tin_corrected.png", width=1024, height=1024)
ggplot(pcaData, aes(PC1, PC2, color=tin, shape=B2)) +
    geom_point(size=5) +
    #geom_text(aes(label=sample),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    theme(text = element_text(size = 24))
dev.off()
#RUV normalization

# Calculate the average gene expression of each genes and take the top 5000 highest expressed
nc <- as.data.frame(assay(vsd))
nc$avg_exp <- rowSums(nc,na.rm=TRUE) / ncol(nc) # x should be your normalized matrix
nc <- arrange(nc, -avg_exp)%>%dplyr::select(-avg_exp)
nc <- nc[1:5000,]
# Calculate the variance of each genes, and choose the lowest 1000 genes as the negative control gene
nc$row_stv <- rowSds(as.matrix(nc))/(rowSums(nc,na.rm=TRUE) / ncol(nc))
nc <- arrange(nc, row_stv) %>% dplyr::select(-row_stv)
nc <- nc[1:1000,]
##The newSeqExpressionSet is given the raw count matrix, phenotypic data, and then the names of the samples 
# Create a new Expression Set and perform an upper quartile normalization
nc <- round(nc*1000)
set <- newSeqExpressionSet(as.matrix(nc),phenoData = data.frame(colData(dds),row.names=colnames(nc)))
#set <- betweenLaneNormalization(set, which="upper")
##Run RUV
##k refers to the number of "unwanted factors of variation" (kinda like PCs in PCA), with rows referring to the samples 
##and columns to the factors. You can play around with k, but 5 is a good starting number.
set_ruv<- RUVg(set, rownames(nc), k=1)

# Save the phenotype data
ruv_data <- pData(set_ruv)
png("results/i/B2/RUV_rle_plot_1.png")
plotRLE(set_ruv, outline=FALSE, ylim=c(-0.1, 0.1), col=as.numeric(as.factor(coldata$B2)))
dev.off()
png('results/i/B2/RUV_pca_plot_1.png')
plotPCA(set, col=as.numeric(as.factor(coldata$B2)), k=2, cex=1.2)
dev.off()
head(ruv_data)
write.csv(ruv_data, file = "results/i/ruv_data_cd.csv")
ruv_data$batch = as.factor(ruv_data$batch)

batch = ruv_data$batch
sexBin = as.numeric(as.factor(ruv_data$sex))
age = ruv_data$dx_age
ruv_data$B3 <- as.factor(ruv_data$B3)
ruv_data$B2 <- as.factor(ruv_data$B2)
ruv_data$rem <- as.factor(ruv_data$rem)
ruv_data$surg <- as.factor(ruv_data$surg)
B2 = ruv_data$B2
B3 = ruv_data$B3
surg = ruv_data$surg
rem = ruv_data$rem
tin = ruv_data$TIN
w_1 = ruv_data$W_1
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + TIN + W_1 + B2)

vsd = vst(dds)
mat= assay(vsd)
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, batch2=sexBin, covariates=cbind(tin + w_1))

pcaData <- plotPCA(vsd, intgroup=c("B2"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#save the plot
png("results/i/B2/pca_B2_tin_w1.png", width=1024, height=1024)
ggplot(pcaData, aes(PC1, PC2, color=B2)) +
    geom_point(size=5) +
    #geom_text(aes(label=sample),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    scale_color_manual(values = c("0" = "#0066cc", "1" = "#990000")) + 
    theme(text = element_text(size = 24))
dev.off()
deg_counts <- data.frame(matrix(ncol=2, nrow=0))
x <- c('complication', 'n_DEGs')
colnames(deg_counts) <- x
write.csv(deg_counts, file = "results/i/deg_counts_cd.csv", row.names = FALSE)
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + TIN + W_1 + B2)
y <- DGEList(counts=counts(dds1), group=coldata$B2)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]

dds1 <- DESeq(dds1)
res <- results(dds1)
resOrdered <- res[order(res$padj),]
#create a column named gene which is the gene name
resOrdered$gene = rownames(resOrdered)
#set gene as the first column by making the last column the first column
resOrdered = resOrdered[,c(7, 1, 2, 3, 4, 5, 6)]

deg_counts[nrow(deg_counts)+1,] <- c('B2', sum(resOrdered$padj < 0.05))

#plot a histogram of the p-values
png("results/i/B2/pvals_hist.png", width=512, height=512)
hist(resOrdered$pvalue, breaks=20, col="grey", main="p-value distribution", xlab="p-value")
dev.off()

write.csv(resOrdered, file = "results/i/B2/results_deseq.csv", sep = ",", row.names = FALSE, col.names = TRUE)
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + TIN + W_1 + B3)
y <- DGEList(counts=counts(dds1), group=coldata$B3)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]

dds1 <- DESeq(dds1)
res <- results(dds1)
resOrdered <- res[order(res$padj),]
#create a column named gene which is the gene name
resOrdered$gene = rownames(resOrdered)
#set gene as the first column by making the last column the first column
resOrdered = resOrdered[,c(7, 1, 2, 3, 4, 5, 6)]

deg_counts[nrow(deg_counts)+1,] <- c('B3', sum(resOrdered$padj < 0.05))

#plot a histogram of the p-values
png("results/i/B3/pvals_hist.png", width=512, height=512)
hist(resOrdered$pvalue, breaks=20, col="grey", main="p-value distribution", xlab="p-value")
dev.off()

write.csv(resOrdered, file = "results/i/B3/results_deseq.csv", sep = ",", row.names = FALSE, col.names = TRUE)
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + TIN + W_1 + rem)
y <- DGEList(counts=counts(dds1), group=coldata$rem)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]

dds1 <- DESeq(dds1)
res <- results(dds1)
resOrdered <- res[order(res$padj),]
#create a column named gene which is the gene name
resOrdered$gene = rownames(resOrdered)
#set gene as the first column by making the last column the first column
resOrdered = resOrdered[,c(7, 1, 2, 3, 4, 5, 6)]

deg_counts[nrow(deg_counts)+1,] <- c('rem', sum(resOrdered$padj < 0.05))

#plot a histogram of the p-values
png("results/i/rem/pvals_hist.png", width=512, height=512)
hist(resOrdered$pvalue, breaks=20, col="grey", main="p-value distribution", xlab="p-value")
dev.off()

write.csv(resOrdered, file = "results/i/rem/results_deseq.csv", sep = ",", row.names = FALSE, col.names = TRUE)
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + TIN + W_1 + surg)
y <- DGEList(counts=counts(dds1), group=coldata$surg)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]

dds1 <- DESeq(dds1)
res <- results(dds1)
resOrdered <- res[order(res$padj),]
#create a column named gene which is the gene name
resOrdered$gene = rownames(resOrdered)
#set gene as the first column by making the last column the first column
resOrdered = resOrdered[,c(7, 1, 2, 3, 4, 5, 6)]

deg_counts[nrow(deg_counts)+1,] <- c('surg', sum(resOrdered$padj < 0.05))

#plot a histogram of the p-values
png("results/i/surg/pvals_hist.png", width=512, height=512)
hist(resOrdered$pvalue, breaks=20, col="grey", main="p-value distribution", xlab="p-value")
dev.off()

write.csv(resOrdered, file = "results/i/surg/results_deseq.csv", sep = ",", row.names = FALSE, col.names = TRUE)
head(deg_counts)

# %%
write.csv(deg_counts, file = "results/i/deg_counts.csv", row.names = FALSE)


