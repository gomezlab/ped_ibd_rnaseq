# %%
library(EnhancedVolcano)
library(fgsea)
library(tidyverse)
library(dplyr)

# %%
res <- read.csv('results/i/cd/results_deseq.csv')
res2 <- res %>%
    dplyr::select(gene, stat) %>%
    na.omit() %>%
    group_by(gene) %>%
    summarize(stat)
library(fgsea)
ranks <- deframe(res2)
head(ranks, 15)

pathways.hallmark <- gmtPathways("gsea/h.all.v7.5.1.symbols.gmt")
pathways.hallmark %>%
    head() %>%
    lapply(head)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
head(fgseaResTidy)

#change leadingEdge data type to string
fgseaResTidy$leadingEdge <- as.character(fgseaResTidy$leadingEdge)
#save to csv
write.csv(fgseaResTidy, 'results/i/cd/fgsea_hallmark.csv')

top10 = fgseaResTidy %>%
    arrange(padj) %>%
    head(15)
#replace '_' in pathway to ' '
top10$pathway <- gsub('_', ' ', top10$pathway)
#strip 'HALLMARK ' from pathway
top10$pathway <- gsub('HALLMARK ', '', top10$pathway)
#make pathway all lowercase
top10$pathway <- tolower(top10$pathway)

head(top10)
png('results/i/cd/fgsea_hallmark.png', width=1440, height=1440)
ggplot(top10, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj < 0.05)) +
  scale_fill_manual(values=c("#0066cc", "#990000")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways", subtitle='Ileum') + 
  theme_minimal(base_size=45)
dev.off()

# %%
res <- read.csv('results/c/cd/results_deseq.csv')

head(res)
res2 <- res %>%
    dplyr::select(gene, stat) %>%
    na.omit() %>%
    group_by(gene) %>%
    summarize(stat)
library(fgsea)
ranks <- deframe(res2)
head(ranks, 15)

pathways.hallmark <- gmtPathways("gsea/h.all.v7.5.1.symbols.gmt")
pathways.hallmark %>%
    head() %>%
    lapply(head)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
head(fgseaResTidy)

#change leadingEdge data type to string
fgseaResTidy$leadingEdge <- as.character(fgseaResTidy$leadingEdge)
#save to csv
write.csv(fgseaResTidy, 'results/c/cd/fgsea_hallmark.csv')

top10 = fgseaResTidy %>%
    arrange(padj) %>%
    head(15)
#replace '_' in pathway to ' '
top10$pathway <- gsub('_', ' ', top10$pathway)
#strip 'HALLMARK ' from pathway
top10$pathway <- gsub('HALLMARK ', '', top10$pathway)
#make pathway all lowercase
top10$pathway <- tolower(top10$pathway)

head(top10)
png('results/c/cd/fgsea_hallmark.png', width=1440, height=1440)
ggplot(top10, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj < 0.05)) +
  scale_fill_manual(values=c("#0066cc", "#990000")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways", subtitle='Colon') + 
  theme_minimal(base_size=45)
dev.off()

# %%
res <- read.csv('results/c/B2/results_deseq.csv')
res$X = res$gene
head(res)
res2 <- res %>%
    dplyr::select(gene, stat) %>%
    na.omit() %>%
    group_by(gene) %>%
    summarize(stat)
library(fgsea)
ranks <- deframe(res2)
head(ranks, 15)

# %%
pathways.hallmark <- gmtPathways("/home/kchen/Documents/ped_ibd/deseq/gsea/h.all.v7.5.1.symbols.gmt")
pathways.hallmark %>%
    head() %>%
    lapply(head)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
head(fgseaResTidy)

#change leadingEdge data type to string
fgseaResTidy$leadingEdge <- as.character(fgseaResTidy$leadingEdge)
#save to csv
write.csv(fgseaResTidy, 'results/c/B2/fgsea_hallmark.csv')

top10 = fgseaResTidy %>%
    arrange(NES) %>%
    head(15)
#replace '_' in pathway to ' '
top10$pathway <- gsub('_', ' ', top10$pathway)
#strip 'HALLMARK ' from pathway
top10$pathway <- gsub('HALLMARK ', '', top10$pathway)
#make pathway all lowercase
top10$pathway <- tolower(top10$pathway)

head(top10)
png('results/c/B2/fgsea_hallmark.png', width=1440, height=1440)
ggplot(top10, aes(reorder(pathway, NES), NES)) +
  #fill color based on padj < 0.05, where padj < 0.05 is blue, otherwise red
  geom_col(aes(fill=padj < 0.05)) +
  scale_fill_manual(values=c("#0066cc", "#990000")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark Pathways", subtitle='Colon') + 
  #increase font size
  theme_minimal(base_size=45)
dev.off()
# %%


