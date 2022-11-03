# %%
## Ben's pipeline for RNA-seq instructions 
library("AnnotationDbi")
library("org.Hs.eg.db")
library(edgeR)
library(RUVSeq)

# This is useful for when you want to run in a different directory
# Should be updated based on where you are running this
source("downstreamAnalysis_RNAseqFunctions.R")

# %%
# Read in data
cts = as.matrix(read.csv("../../data/processed_data/data_c_all_raw.csv", sep=',', row.names='gene'))

coldata = read.csv('../../data/metadata/coldata_all_c.csv', sep=',', row.names='sample')
dim(coldata)
dim(cts)
#get coldata row names
coldata$sample = rownames(coldata)
sample = rownames(coldata)
#set columns as factor
status = as.factor(coldata$status)
batch = coldata$batch
sexBin = as.numeric(as.factor(coldata$sex))
age = coldata$DiagnosisAge
status = as.factor(coldata$status)

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
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, covariates=cbind(sexBin))

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
# png("../results/B2/i/RUV_plot_1.png")
# plotRLE(set_ruv, outline=FALSE, ylim=c(-0.1, 0.1), col=as.numeric(as.factor(coldata$subtype)))
# dev.off()
# png("../results/B2/i/RUV_pca_1.png")
# Creates a PCA plot of the negative control genes. Should not see much of a pattern. Write out learned RUVseq factors 
# plotPCA(set, col=as.numeric(as.factor(coldata$status)), k=2, cex=1.2)
# dev.off()

# Save the phenotype data
ruv_data<-pData(set_ruv)
#write.table(ruv_data, file = "../data/ruvdata.csv", sep = ",", row.names = TRUE, col.names = TRUE)

coldata$status
dds$status <- factor(dds$status, levels = c("nonIBD","CD"))
dds$status <- relevel(dds$status, ref='nonIBD')
levels(dds$status)
# Set up new DESeq object to perform differential analyses
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + W_1 + status)

# remove genes with low counts for most samples
y <- DGEList(counts=counts(dds), group=coldata$status)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]

# %%
d <- plotCounts(dds1, gene='CADM1', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/ibd_vs_non/c/CADM1_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="CADM1", fill='Subtype') +
            scale_color_manual(values=c("#00BFC4", "#F8766D")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='HCAR3', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/ibd_vs_non/c/HCAR3_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="HCAR3", fill='Subtype') +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='IL24', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
# png(file="results/ibd_vs_non/c/IL24_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="IL24", fill='Subtype') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
# dev.off()

# %%
d <- plotCounts(dds1, gene='AQP9', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/ibd_vs_non/c/AQP9_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="AQP9", fill='Subtype') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='UCN2', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/ibd_vs_non/c/UCN2_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="UCN2", fill='Subtype') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='MMP3', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/ibd_vs_non/c/MMP3_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="MMP3", fill='Subtype') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()


# %%
d <- plotCounts(dds1, gene='MUC6', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/ibd_vs_non/c/MUC6_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="MUC6", fill='Subtype') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
# Read in data
cts = as.matrix(read.csv("../../data/processed_data/data_c_raw.csv", sep=',', row.names='gene'))

coldata = read.csv('../data/metadata/coldata_cd_c.csv', sep=',', row.names='sample')

dim(coldata)
dim(cts)


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
design <- as.formula(~ batch + sex + B2)
#create the dds object
dds = DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = design)

vsd = vst(dds)
mat= assay(vsd)

# %%
#remove batch effects
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, batch2=sexBin, covariates=cbind(tin))
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


# %%
# Save the phenotype data
ruv_data <- pData(set_ruv)

design1 <- as.formula(~batch + sex + W_1 + B2)
dds1 <- DESeqDataSetFromMatrix(countData=cts, colData=ruv_data, design = design1)

#remove genes with low counts for most samples
y <- DGEList(counts=counts(dds1), group=B2)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]

# %%
d <- plotCounts(dds1, gene='FGL2', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/c/B2/FGL2_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="FGL2", fill='Subtype') +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()


# %%
d <- plotCounts(dds1, gene='REG1A', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/c/B2/REG1A_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="REG1A", fill='Subtype') +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='MMP3', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/c/B2/MMP3_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="MMP3", fill='Subtype') +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='DUOX2', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/c/B2/DUOX2_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="DUOX2", fill='Subtype') +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='MMP1', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/c/B2/MMP1_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(10, 100, 1000, 10000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="MMP1", fill=B2) +
            # scale_fill_brewer(palette="Set2") +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%


# %%
d <- plotCounts(dds1, gene='UCN2', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/c/B2/UCN2_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(10, 100, 1000, 10000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="UCN2", fill=B2) +
            # scale_fill_brewer(palette="Set2") +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='CXCL5', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/c/B2/CXCL5_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(10, 100, 1000, 10000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="CXCL5", fill=B2) +
            # scale_fill_brewer(palette="Set2") +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='CHI3L1', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/c/B2/CHI3L1_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(10, 100, 1000, 10000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="CHI3L1", fill=B2) +
            # scale_fill_brewer(palette="Set2") +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='AQP9', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/c/B2/AQP9_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(10, 100, 1000, 10000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="AQP9", fill=B2) +
            # scale_fill_brewer(palette="Set2") +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='HCAR3', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/c/B2/HCAR3_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(10, 100, 1000, 10000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="HCAR3", fill=B2) +
            # scale_fill_brewer(palette="Set2") +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='S100A8', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/c/B2/S100A8_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(10, 100, 1000, 10000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="S100A8", fill=B2) +
            # scale_fill_brewer(palette="Set2") +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='TIMP1', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/c/B2/TIMP1_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(10, 100, 1000, 10000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="TIMP1", fill=B2) +
            # scale_fill_brewer(palette="Set2") +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
#adjust for batch and sexBin
design <- as.formula(~ batch + sex + B2)
#create the dds object
dds = DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = design)
# #remove genes with low counts for most samples
# library(edgeR)
# y <- DGEList(counts=counts(dds), group=coldata$status)
# keep <- filterByExpr(y)
# table(keep)
# y <- y[keep,]
# dds <- dds[keep,]
vsd = vst(dds)
mat= assay(vsd)


#remove batch effects
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, covariates=cbind(sexBin))
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

design1 <- as.formula(~batch + sex + W_1 + B2)
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = design1)

#remove genes with low counts for most samples
y <- DGEList(counts=counts(dds1), group=B2)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]

# %%
d <- plotCounts(dds1, gene='HCAR3', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/B2/c/HCAR3_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(10, 100, 1000, 10000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="HCAR3", fill=B2) +
            # scale_fill_brewer(palette="Set2") +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()


# %%
d <- plotCounts(dds1, gene='MMP3', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/B2/c/MMP3_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(10, 100, 1000, 10000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="MMP3", fill='B2') +
            # scale_fill_brewer(palette="Set2") +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            # theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                # axis.text=element_text(size=14, color="black"),
                # panel.background = element_rect(fill="white",color="black"),
                # legend.title=element_blank(), legend.key=element_blank())
dev.off()


# %%
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/B2/c/AQP9_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="AQP9", fill='B2') +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()
d <- plotCounts(dds1, gene='UCN2', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/B2/c/UCN2_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="UCN2", fill='B2') +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()
d <- plotCounts(dds1, gene='CXCL5', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/B2/c/CXCL5_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="CXCL5", fill='B2') +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()
d <- plotCounts(dds1, gene='CHI3L1', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/B2/c/CHI3L1_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="CHI3L1", fill='B2') +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()
d <- plotCounts(dds1, gene='IL24', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/B2/c/IL24_c.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="IL24", fill='B2') +
            scale_color_manual(values=c("#0066cc", "#990000")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
# Read in data
cts = as.matrix(read.csv("../../data/processed_data/data_i_all_raw.csv", sep=',', row.names='gene'))

coldata = read.csv('../../data/metadata/coldata_all_i.csv', sep=',', row.names='sample')
dim(coldata)
dim(cts)
#get coldata row names
coldata$sample = rownames(coldata)
sample = rownames(coldata)
#set columns as factor
status = as.factor(coldata$status)
batch = coldata$batch
sexBin = as.numeric(as.factor(coldata$sex))
age = coldata$DiagnosisAge
status = as.factor(coldata$status)

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
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, covariates=cbind(sexBin))

# Calculate the average gene expression of each genes and take the top 5000 highest expressed
nc <- as.data.frame(assay(vsd))
nc$avg_exp <- rowSums(nc,na.rm=TRUE) / ncol(nc) # x should be your normalized matrix
nc <- arrange(nc, -avg_exp)%>%dplyr::select(-avg_exp)
nc <- nc[1:1000,]

# Calculate the variance of each genes, and choose the lowest 1000 genes as the negative control gene
nc$row_stv <- rowSds(as.matrix(nc))/(rowSums(nc,na.rm=TRUE) / ncol(nc))
nc <- arrange(nc, row_stv) %>% dplyr::select(-row_stv)
nc <- nc[1:250,]
##The newSeqExpressionSet is given the raw count matrix, phenotypic data, and then the names of the samples 
# Create a new Expression Set and perform an upper quartile normalization
nc <- round(nc*1000)
set <- newSeqExpressionSet(as.matrix(nc),phenoData = data.frame(colData(dds),row.names=colnames(nc)))
#set <- betweenLaneNormalization(set, which="upper")
##Run RUV
##k refers to the number of "unwanted factors of variation" (kinda like PCs in PCA), with rows referring to the samples 
##and columns to the factors. You can play around with k, but 5 is a good starting number.
set_ruv<- RUVg(set, rownames(nc), k=1)
# png("../results/B2/i/RUV_plot_1.png")
# plotRLE(set_ruv, outline=FALSE, ylim=c(-0.1, 0.1), col=as.numeric(as.factor(coldata$subtype)))
# dev.off()
# png("../results/B2/i/RUV_pca_1.png")
# Creates a PCA plot of the negative control genes. Should not see much of a pattern. Write out learned RUVseq factors 
# plotPCA(set, col=as.numeric(as.factor(coldata$status)), k=2, cex=1.2)
# dev.off()

# Save the phenotype data
ruv_data<-pData(set_ruv)
#write.table(ruv_data, file = "../data/ruvdata.csv", sep = ",", row.names = TRUE, col.names = TRUE)

coldata$status
dds$status <- factor(dds$status, levels = c("nonIBD","CD"))
dds$status <- relevel(dds$status, ref='nonIBD')
levels(dds$status)
# Set up new DESeq object to perform differential analyses
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = ~batch + sex + W_1 + status)

# remove genes with low counts for most samples
y <- DGEList(counts=counts(dds), group=coldata$status)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]


# %%
d <- plotCounts(dds1, gene='CADM1', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/ibd_vs_non/i/CADM1_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="CADM1", fill='Subtype') +
            scale_color_manual(values=c("#00BFC4", "#F8766D")) +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%
d <- plotCounts(dds1, gene='HCAR3', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/ibd_vs_non/i/HCAR3_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="HCAR3", fill='Subtype') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()
d <- plotCounts(dds1, gene='IL24', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
# png(file="results/ibd_vs_non/i/IL24_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="IL24", fill='Subtype') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
# dev.off()
d <- plotCounts(dds1, gene='AQP9', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/ibd_vs_non/i/AQP9_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="AQP9", fill='Subtype') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()
d <- plotCounts(dds1, gene='UCN2', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/ibd_vs_non/i/UCN2_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="UCN2", fill='Subtype') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()
d <- plotCounts(dds1, gene='MMP3', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/ibd_vs_non/i/MMP3_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="MMP3", fill='Subtype') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

d <- plotCounts(dds1, gene='MUC6', intgroup='status', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/ibd_vs_non/i/MUC6_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=status, y=count, color=status)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="MUC6", fill='Subtype') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()


# %%

# Read in data
cts = as.matrix(read.csv("../../data/processed_data/data_i_raw.csv", sep=',', row.names='gene'))

coldata = read.csv('../../data/metadata/coldata_cd_i.csv', sep=',', row.names='sample')

dim(coldata)
dim(cts)

#get coldata row names
coldata$sample = rownames(coldata)
sample = rownames(coldata)
#set columns as factor
#coldata$status = as.factor(coldata$status)
#status = as.factor(coldata$status)
batch = coldata$batch
sexBin = as.numeric(as.factor(coldata$sex))
tissue = coldata$Tissue
age = coldata$DiagnosisAge
coldata$B3 <- as.factor(coldata$B3)
coldata$B2 <- as.factor(coldata$B2)
coldata$rem <- as.factor(coldata$rem)
B2 = coldata$B2
B3 = coldata$B3
coldata$surg = as.factor(coldata$surg)
surg = coldata$surg
rem = coldata$rem
#round counts to integers
cts <- round(cts);

#check if row names of coldata match colnames of cts
all(rownames(coldata) == colnames(cts))
#adjust for batch and sexBin
design <- as.formula(~ batch + sex)
#create the dds object
dds = DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = design)
# #remove genes with low counts for most samples
# library(edgeR)
# y <- DGEList(counts=counts(dds), group=coldata$status)
# keep <- filterByExpr(y)
# table(keep)
# y <- y[keep,]
# dds <- dds[keep,]
vsd = vst(dds)
mat= assay(vsd)


#remove batch effects
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, covariates=cbind(sexBin))
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

design1 <- as.formula(~batch + sex + W_1 + B2)
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=ruv_data,design = design1)

#remove genes with low counts for most samples
y <- DGEList(counts=counts(dds1), group=B2)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds1 <- dds1[keep,]
d <- plotCounts(dds1, gene='HCAR3', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/B2/i/HCAR3_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="Subtype", y="Normalized/Transformed Expression", title="HCAR3", fill='Subtype') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()
d <- plotCounts(dds1, gene='AQP9', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/B2/i/AQP9_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="AQP9", fill='B2') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()
d <- plotCounts(dds1, gene='UCN2', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/B2/i/UCN2_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="UCN2", fill='B2') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()
d <- plotCounts(dds1, gene='MMP3', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/B2/i/MMP3_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="MMP3", fill='B2') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

d <- plotCounts(dds1, gene='MUC6', intgroup='B2', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/B2/i/MUC6_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=B2, y=count, color=B2)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="B2", y="Normalized/Transformed Expression", title="MUC6", fill='B2') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()
d <- plotCounts(dds1, gene='MUC6', intgroup='surg', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/surg/i/MUC6_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=surg, y=count, color=surg)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="surg", y="Normalized/Transformed Expression", title="MUC6", fill='surg') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()
d <- plotCounts(dds1, gene='UCN2', intgroup='surg', returnData=TRUE)
#plot the scatter plot with group 0 renamed as 'non-IBD' and group 1 as 'CD'
png(file="results/surg/i/UCN2_i.png", width=2048, height=2048, res=300)
ggplot(d, aes(x=surg, y=count, color=surg)) +
            geom_point(position=position_jitter(w=0.1,h=0)) +
            geom_boxplot(alpha = 0) +
    # stat_summary(fun.y=mean, geom="line", aes(y=count, group=1), colour='black', group=1) +
    # stat_summary(fun.y=mean, geom="point", colour='black', size=3, shape=8) +
            scale_y_log10(breaks=c(100, 200, 500, 1000)) +
            labs(x="surg", y="Normalized/Transformed Expression", title="UCN2", fill='surg') +
            theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
dev.off()

# %%



