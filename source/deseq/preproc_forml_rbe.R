# %%
## Ben's pipeline for RNA-seq instructions 
library("AnnotationDbi")
library("org.Hs.eg.db")
library(RUVSeq)

# This is useful for when you want to run in a different directory
# Should be updated based on where you are running this
source("/home/kchen/Documents/ped_ibd/data/raw/downstreamAnalysis_RNAseqFunctions.R")
setwd("/home/kchen/Documents/ped_ibd/deseq2")
# Read in data
cts = as.matrix(read.csv("../data/data_c_raw.csv", sep=',', row.names='gene'))

coldata = read.csv('../data/clin_c.csv', sep=',', row.names='sample')
dim(coldata)
dim(cts)
#get coldata row names
coldata$sample = rownames(coldata)
sample = rownames(coldata)

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

#remove batch effect
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, batch2 = sexBin, covariates=cbind(tin))

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

ruv_data <- pData(set_ruv)
#set columns as factor
factor_cols = c('batch', 'sex', 'B2', 'B3', 'rem', 'surg')
for (col in factor_cols) {
  ruv_data[,col] = as.factor(ruv_data[,col])
}
batch = ruv_data$batch
sexBin = ruv_data$sex
age = ruv_data$dx_age

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
mat= assay(vsd)

write.csv(mat, file='../data/data_vst_ruv1_c.csv')

# %%
# Read in data
cts = as.matrix(read.csv("../data/data_i_raw.csv", sep=',', row.names='gene'))

coldata = read.csv('../data/clin_i.csv', sep=',', row.names='sample')
dim(coldata)
dim(cts)
#get coldata row names
coldata$sample = rownames(coldata)
sample = rownames(coldata)

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

#remove batch effect
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, batch2 = sexBin, covariates=cbind(tin))

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

ruv_data <- pData(set_ruv)
#set columns as factor
factor_cols = c('batch', 'sex', 'B2', 'B3', 'rem', 'surg')
for (col in factor_cols) {
  ruv_data[,col] = as.factor(ruv_data[,col])
}
batch = ruv_data$batch
sexBin = ruv_data$sex
age = ruv_data$dx_age

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
mat= assay(vsd)
write.csv(mat, file='../data/data_vst_ruv1_i.csv')

# %%



