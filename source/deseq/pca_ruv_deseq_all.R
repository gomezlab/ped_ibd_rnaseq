# %%
## Ben's pipeline
library("AnnotationDbi")
library("org.Hs.eg.db")
library(RUVSeq)
library('openxlsx')
source("downstreamAnalysis_RNAseqFunctions.R")

# %%
# Read in data
cts = as.matrix(read.csv("../../data/processed_data/data_c_all_raw.csv", sep=',', row.names='gene'))
coldata = read.csv('../../data/metadata/coldata_all_c.csv', sep=',', row.names='sample')

#remove sample PNIBD26_C from cts and coldata
cts = cts[,-which(colnames(cts) == 'PNIBD26_C')]
coldata = coldata[-which(rownames(coldata) == 'PNIBD26_C'),]
dim(coldata)
dim(cts)

#create a new column for the sample name
coldata$sample = rownames(coldata)
sample = rownames(coldata)

#set columns as factor
coldata$status = as.factor(coldata$status)
coldata$batch = as.factor(coldata$batch)
coldata$sex = as.factor(coldata$sex)

#make variables based on co-variates for use in removebatcheffect
status = coldata$status
batch = coldata$batch
sexBin = coldata$sex
tin = coldata$TIN

#round counts to integers
cts <- round(cts);
#check if row names of coldata match colnames of cts
all(rownames(coldata) == colnames(cts))

#adjust for batch and sex
design <- as.formula(~ batch + sex)

#create the dds object
dds = DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = design)

#variance stabilizing transform
vsd = vst(dds)
mat= assay(vsd)


# %%
#PCA of raw data
pcaData <- plotPCA(vsd, intgroup=c("status"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#save the plot
png("results/c/cd/pca_raw.png", width=1024, height=1024)
ggplot(pcaData, aes(PC1, PC2, color=tin, shape=status)) +
    geom_point(size=5) +
    #geom_text(aes(label=sample),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    
    theme(text = element_text(size = 24))
dev.off()

# %%
#remove batch effect
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, batch2 = sexBin, covariates=cbind(tin))

# %%
#PCA of batch effect removed data
pcaData <- plotPCA(vsd, intgroup=c("status"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#save the plot
png("results/c/cd/pca_tin_corrected.png", width=1024, height=1024)
ggplot(pcaData, aes(PC1, PC2, color=tin, shape=status)) +
    geom_point(size=5) +
    #geom_text(aes(label=sample),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    
    theme(text = element_text(size = 24))
dev.off()

# %%
#RUVseq

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

##Run RUV
##k refers to the number of "unwanted factors of variation", set as 1 based on correlation with outcome and deseq results
set_ruv<- RUVg(set, rownames(nc), k=1)


# %%
#RLE plot
png("results/c/cd/RUV_rle_plot_1.png")
plotRLE(set_ruv, outline=FALSE, ylim=c(-0.1, 0.1), col=as.numeric(as.factor(coldata$status)))
dev.off()

# %%
# Creates a PCA plot based on negative control genes. Should not see much of a pattern.
png("results/c/cd/RUV_pca_1.png")
plotPCA(set, col=as.numeric(as.factor(coldata$status)), k=2, cex=1.2)
dev.off()

# %%
#save the coldata + RUV factor
ruv_data<-pData(set_ruv)
write.xlsx(ruv_data, file="results/c/cd/ruv_data_1.xlsx")

# %%
#redo data-prep steps with the ruv_data dataframe
#set columns as factor
ruv_data$status = as.factor(ruv_data$status)
ruv_data$batch = as.factor(ruv_data$batch)
ruv_data$sex = as.factor(ruv_data$Sex)

status = ruv_data$status
batch = ruv_data$batch
sexBin = ruv_data$Sex

tin = ruv_data$TIN

w_1 = ruv_data$W_1

# %%
# Set up new DESeq object to perform differential analyses
dds1 <- DESeqDataSetFromMatrix(countData=cts, colData=ruv_data, design = ~batch + Sex + TIN + W_1 + status)

# %%
vsd = vst(dds)
mat= assay(vsd)
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, batch2=sexBin, covariates=cbind(tin + w_1))

#save the PCA, accounting for batch effect, sex, tin score, and the RUV factor (W_1)
pcaData <- plotPCA(vsd, intgroup=c("status"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("results/c/cd/pca_cd_tin_w1_new.png", width=1024, height=1024)
ggplot(pcaData, aes(PC1, PC2, color=status)) +
    geom_point(size=5) +
    # geom_text(aes(label=sample),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    scale_color_manual(values = c("nonIBD" = "#0066cc", "CD" = "#990000")) + 
    theme(text = element_text(size = 24))
dev.off()


# %%
#now the DESeq analysis
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + Sex + TIN + W_1 + status)
#filter out lowly expressed genes
y <- DGEList(counts=counts(dds1), group=coldata$status)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]
#run DESeq
dds1 <- DESeq(dds1)
#get the results, with nonIBD as the reference
res <- results(dds1, contrast=c('status', 'CD', 'nonIBD'))
#order the results by adjusted p-value
resOrdered <- res[order(res$padj),]
#create a column named gene which is the gene name
resOrdered$gene = rownames(resOrdered)
#set gene as the first column by making the last column the first column
resOrdered = resOrdered[,c(7, 1, 2, 3, 4, 5, 6)]

#plot a histogram of the p-values
png("results/c/cd/pvals_hist.png", width=512, height=512)
hist(resOrdered$pvalue, breaks=20, col="grey", main="p-value distribution", xlab="p-value")
dev.off()

#save the results
write.csv(resOrdered, file = "results/c/cd/results_deseq.csv", sep = ",", row.names = FALSE, col.names = TRUE)

#save the number of genes that are significant
deg_counts <- data.frame('colon' = c(sum(resOrdered$padj < 0.05)))
write.xlsx(deg_counts, file = "results/c/cd/deg_counts.xlsx", sheetName='DEGs', row.names = TRUE, col.names = TRUE)


# %%


# %%
#same thing for the ileum data
# Read in data
cts = as.matrix(read.csv("../data/processed_data/data_i_all_raw.csv", sep=',', row.names='gene'))
coldata = read.csv('../data/metadata/coldata_all_i.csv', sep=',', row.names='sample')
dim(coldata)
dim(cts)
coldata$sample = rownames(coldata)
sample = rownames(coldata)
#set columns as factor
coldata$status = as.factor(coldata$status)
coldata$batch = as.factor(coldata$batch)
coldata$sex = as.factor(coldata$sex)

status = coldata$status
batch = coldata$batch
sexBin = coldata$sex

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
design0 <- model.matrix(~coldata$status, data=data.frame(assay(vsd)))

pcaData <- plotPCA(vsd, intgroup=c("status"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#save the plot
png("results/i/cd/pca_raw.png", width=1024, height=1024)
ggplot(pcaData, aes(PC1, PC2, color=tin, shape=status)) +
    geom_point(size=5) +
    #geom_text(aes(label=sample),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    
    theme(text = element_text(size = 24))
dev.off()

assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, batch2 = sexBin, covariates=cbind(tin))
pcaData <- plotPCA(vsd, intgroup=c("status"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#save the plot
png("results/i/cd/pca_tin_corrected.png", width=1024, height=1024)
ggplot(pcaData, aes(PC1, PC2, color=tin, shape=status)) +
    geom_point(size=5) +
    #geom_text(aes(label=sample),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    
    theme(text = element_text(size = 24))
dev.off()

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

png("results/i/cd/RUV_rle_plot_1.png")
plotRLE(set_ruv, outline=FALSE, ylim=c(-0.1, 0.1), col=as.numeric(as.factor(coldata$status)))
dev.off()
# Creates a PCA plot of the negative control genes. Should not see much of a pattern. Write out learned RUVseq factors 
png("results/i/cd/RUV_pca_1.png")
plotPCA(set, col=as.numeric(as.factor(coldata$status)), k=2, cex=1.2)
dev.off()
ruv_data<-pData(set_ruv)
write.csv(ruv_data, file="results/i/cd/ruv_data_1.csv")

#set columns as factor
ruv_data$status = as.factor(ruv_data$status)
ruv_data$batch = as.factor(ruv_data$batch)
ruv_data$sex = as.factor(ruv_data$sex)

status = ruv_data$status
batch = ruv_data$batch
sexBin = ruv_data$sex

tin = ruv_data$TIN

w_1 = ruv_data$W_1

# Set up new DESeq object to perform differential analyses
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + TIN + W_1 + status)
vsd = vst(dds)
mat= assay(vsd)
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, batch2=sexBin, covariates=cbind(tin + w_1))

pcaData <- plotPCA(vsd, intgroup=c("status"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#save the plot
png("results/i/cd/pca_cd_tin_w1.png", width=1024, height=1024)
ggplot(pcaData, aes(PC1, PC2, color=status)) +
    geom_point(size=5) +
    #geom_text(aes(label=sample),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    scale_color_manual(values = c("nonIBD" = "#0066cc", "CD" = "#990000")) + 
    theme(text = element_text(size = 24))
dev.off()
# Set up new DESeq object to perform differential analyses
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + TIN + W_1 + status)
y <- DGEList(counts=counts(dds1), group=coldata$status)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]

dds1 <- DESeq(dds1)
res <- results(dds1, contrast=c('status', 'CD', 'nonIBD'))
resOrdered <- res[order(res$padj),]
#create a column named gene which is the gene name
resOrdered$gene = rownames(resOrdered)
#set gene as the first column by making the last column the first column
resOrdered = resOrdered[,c(7, 1, 2, 3, 4, 5, 6)]



#plot a histogram of the p-values
png("results/i/cd/pvals_hist.png", width=512, height=512)
hist(resOrdered$pvalue, breaks=20, col="grey", main="p-value distribution", xlab="p-value")
dev.off()

write.csv(resOrdered, file = "results/i/cd/results_deseq.csv", sep = ",", row.names = FALSE, col.names = TRUE)

deg_counts = data.frame(matrix(ncol=2, nrow=0))
col_names = c('tissue', 'n_DEGs')
colnames(deg_counts) = col_names

sum(resOrdered$padj < 0.05)
deg_counts[nrow(deg_counts)+1,] <- c('ileum', sum(resOrdered$padj < 0.05))

write.xlsx(deg_counts, file = "results/i/cd/deg_counts.xlsx", sheetName='DEGs', row.names = TRUE, col.names = TRUE)