# Linking Gene Expression to Clinical Outcomes in Pediatric Crohn’s Disease Using Machine Learning 

This repository contains code to accompany the manuscript "Linking Gene Expression to Clinical Outcomes in Pediatric Crohn’s Disease Using Machine Learning" and reproduce the results.

There are two main sections to the code. The first, in the folder "deseq" performs differential expression analysis and prepares data for machine learning. The second, in the "ml" folder, contains scripts that perform hyperparameter searches and then train models for the outcomes.

##Data
The clinical and RNA-seq data for this project is not publicly available, but can be accessed at the following repository: '...'. The clinical data should then be placed in 'data/metadata' and RNA-seq data in 'data/processed_data'.

##Differential expression
Analysis starts with differential expression. The scripts 'pca_ruv_deseq_all.R' and 'pca_ruv_deseq_cd.R' in the folder 'source/deseq/' create PCA plots, and perform RUV/DESeq2 analysis for the CD vs nonIBD and specific clinical outcomes respectively. They produce results files ('results_deseq.csv') which are filed by tissue type 'c' - colon and 'i' - ileum and outcome 'cd' - CD vs nonIBD, B2, B3, remission, and surgery. The results files can then be fed into 'hallmark.R' and 'volcano.R' to perform pathway analysis and plot volcano plots. The script 'boxplot.R' creates boxplots for specific genes between the groups of interest.

Finally, 'preproc_forml_rbe.R' can be used to transform the raw counts data for downstream use in machine learning, creating the files 'data_vst_ruv1_c.csv' and 'data_vst_ruv1_i.csv'. 'preproc_comb.ipynb' then combines the RNA-seq data with the clinical data in the products 'comb_vst_i.csv' and 'comb_vst_c.csv'.

##Machine learning
