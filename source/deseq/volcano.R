# %%
library(EnhancedVolcano)
library(tidyverse)
library(dplyr)

# %%
res <- read.csv('results/i/cd/results_deseq.csv')
#if res$'gene' starts with 'ENSG' then drop the row and save to new dataframe res2
res <- res %>% 
  filter(!grepl('ENSG', gene))

res <- res[res$log2FoldChange > -6.5,]
res <- res[res$log2FoldChange < 6.5,]

png('results/i/cd/volcano.png', width=1440, height=1440)
EnhancedVolcano(res,
    titleLabSize = 40,
    subtitleLabSize = 30,
    captionLabSize = 30,
    axisLabSize = 30,
    legendLabSize = 30,
    title = "CD vs non-IBD",
    subtitle = bquote(Ileum),
    lab = res$gene,
    selectLab = c('HCAR3', 'MUC1', 'CXCL8', 'AQP9', 'MMP1', 'CXCL5', 'DUOX2', 'CHI3L1', 'MUC6', 'IL1B', 'TM4SF4', 'GSTA2'),
    x = 'log2FoldChange',
    y = 'padj',
    legendLabels = c("NS", "p-adj", "L2FC", "p < 0.05 and L2FC > 1.5"),
    pCutoff = 0.05,
    FCcutoff = 1.5,
    pointSize = 3.0,
    labSize = 16.0,
    col = c("#e6f2ff", "#b3d9ff", "#b3d9ff", "#0066cc"))
    # colGradient = c('#990000', '#0066cc'))
dev.off()

# %%
head(res)

# %%
#find genes with log2FoldChange is < -1
head(res[res$log2FoldChange < -1,], 5)

# %%
#find genes with log2FoldChange is between 3 and 4
head(res[(res$log2FoldChange > 4) & (res$padj > 1e-10),], 5)

# %%
res <- read.csv('results/c/cd/results_deseq.csv')
#if res$'gene' starts with 'ENSG' then drop the row and save to new dataframe res2
res <- res %>% 
  filter(!grepl('ENSG', gene))

res <- res[res$log2FoldChange > -6.5,]
res <- res[res$log2FoldChange < 6.5,]

png('results/c/cd/volcano.png', width=1440, height=1440)
EnhancedVolcano(res,
    titleLabSize = 40,
    subtitleLabSize = 30,
    captionLabSize = 30,
    axisLabSize = 30,
    legendLabSize = 30,
    title = "CD vs non-IBD",
    subtitle = bquote(Colon),
    lab = res$gene,
    selectLab = c('HCAR3', 'CXCL8', 'AQP9', 'CXCL5', 'DUOX2', 'CHI3L1', 'MUC6', 'IL1B', 'IL1A', 'REG1B', 'AKR1B10', 'UGDH'),
    x = 'log2FoldChange',
    y = 'padj',
    legendLabels = c("NS", "p-adj", "L2FC", "p < 0.05 and L2FC > 1.5"),
    pCutoff = 0.05,
    FCcutoff = 1.5,
    pointSize = 3.0,
    labSize = 16.0,
    col = c("#e6f2ff", "#b3d9ff", "#b3d9ff", "#0066cc"))
    # colGradient = c('#990000', '#0066cc'))
dev.off()

# %%
res <- read.csv('results/c/B2/results_deseq.csv')
# if res$'gene' starts with 'ENSG' then drop the row and save to new dataframe res2
res <- res %>% 
  filter(!grepl('ENSG', gene))
res$X = res$gene
head(res)


# %%
#drop rows with padj is NA
res <- res[!is.na(res$padj),]

# %%
res_pos <- res[res$log2FoldChange > 2,]
head(res_pos)

# %%
# res <- res[res$log2FoldChange > -7.5,]
# res <- res[res$log2FoldChange < 7.5,]
png('results/c/B2/volcano.png', width=1440, height=1440)
EnhancedVolcano(res,
    lab = res$X,
    selectLab = c('MMP1', 'CXCL5', 'ZNF527', 'CXCL8', 'FOXP3', 'INHBA', 'IL1B', 'TIMP1', 'DUOXA2', 'CXCL1', 'MTCO1P12', 'CLDN8'),
    # selectLab = c('MMP3', 'UCN2', 'MMP1', 'CXCL5', 'CXCL8', 'IL24', 'CHI3L1', 'TGIF2-RAB5IF', 'FOXP3', 'AQP9', 'HCAR2', 'REG1A', 'CXCL6', 'MUC6', 'IL1B', 'MTCO1P12', 'CLDN8'),
    # selectLab = c('FOXP3', 'CHI3L1',  'REG1A', 'MTCO1P12'),
    x = 'log2FoldChange',
    y = 'padj',
    legendLabels = c("NS", "p-adj", "L2FC", "p < 0.05 and L2FC > 1.5"),
    titleLabSize = 40,
    subtitleLabSize = 30,
    captionLabSize = 30,
    axisLabSize = 30,
    legendLabSize = 30,
    title = 'B2 - Stricturing',
    subtitle = bquote(Colon),
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 4.0,
    labSize = 16.0,
    ylim = c(0, 4.8),
    xlim = c(-5.5, 5),
    # drawConnectors = TRUE,
    col = c("#e6f2ff", "#b3d9ff", "#b3d9ff", "#0066cc"))
    
dev.off()


