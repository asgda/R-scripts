#automate the script
data_mutation<- read.csv("D:\\COSMIC_DB_FAT4_PC-cell-lines\\V95_38_CLP_MUTANT.csv", header = T, as.is = T, check.names = F)
#data_mutation<- read.delim("D:\\paad_qcmg_uq_2016\\paad_qcmg_uq_2016\\data_mutations_mskcc.txt", as.is = TRUE, check.names = F)
library(dplyr)
library(survival)
library(survminer)
library(tidyr)
#sort(table(data_mutation$Hugo_Symbol), decreasing = TRUE)[1:25] -get the most recurrent genes
#my_list <- list("KRAS", "TP53", "MUC16", "SMAD4", "OBSCN", "PLEC", "SYNE1", "RYR1", "RYR3", "FAT4", "COL5A1", "DST") 
#my_list
#for (i in 1:length(my_list)) {
#gene<-my_list[i]
#print(gene)
#names(data_patient)<- NULL #remove existing header
#data_patient<- data_patient[-c(1:3),] #remove the first 3 rows not needed, only keep the needed header as a row value
#names(data_patient)<- data_patient[1,] #replace the header with value in the 1st row
#data_patient<- data_patient[-c(1),] #after renaming the header, remove the 1st row that had column names
data_mutation<- data_mutation[, c(1,5,6,22)]
data_mutation_total<- data.frame(data_mutation$ID_SAMPLE)
data_mutation$STATUS<- ifelse(data_mutation$GENE_NAME=="FAT4", "M", "WT")
mut_gene<- filter(data_mutation, GENE_NAME %in% "FAT4")
mut_gene_uni<- mut_gene %>% distinct(SAMPLE_NAME,.keep_all = TRUE)
colnames(data_mutation_total)<- gsub("data_mutation.ID_SAMPLE", "ID_SAMPLE", colnames(data_mutation_total))
total_gene_uni<- data_mutation_total %>% distinct(ID_SAMPLE,.keep_all = TRUE)
#abc_status<- data_patient[, c(1, 10)]
#pfs_status<- data_patient[, c(1, 36, 37)]
#mut_gene_uni$Tumor_Sample_Barcode<- gsub(".{3}$", "", mut_gene_uni$Tumor_Sample_Barcode)
#os_status$OS_STATUS<- substr(os_status$OS_STATUS, 0, 1)
#pfs_status$PFS_STATUS<- substr(pfs_status$PFS_STATUS, 0, 1)
#class(os_status$OS_STATUS)= "numeric"
#class(os_status$OS_MONTHS)= "numeric"
#class(pfs_status$PFS_STATUS)= "numeric"
#class(pfs_status$PFS_MONTHS)= "numeric"
#mut_gene_uni<- mutate(mut_gene_uni, Status=case_when(Variant_Classification=="Silent" ~ "WT"
#                                                     , Variant_Classification=="Intron" ~ "WT"
#                                                     , Variant_Classification=="UTR" ~ "WT"
#                                                     , TRUE ~ "M")
#)

#colnames(mut_gene_uni)<- gsub("Tumor_Sample_Barcode", "PATIENT_ID", colnames(mut_gene_uni))
table_final<- full_join(total_gene_uni, mut_gene_uni, "ID_SAMPLE")
#table_final1<- table_final[, -c(1,4:8)]
#table_final<- table_final[, -c(2,5)]
#pfs_final<- full_join(pfs_status, mut_gene_uni, "PATIENT_ID")
#pfs_final<- pfs_final[, -c(4,6)]
table_final$STATUS[is.na(table_final$STATUS)]<- "WT"
#pfs_final$Status[is.na(pfs_final$Variant_Classification)]<- "WT"
#fit_os<- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Variant_Classification, data=os_final)
#p_os<-surv_pvalue(fit_os)$pval
#fit_pfs<- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ Variant_Classification, data=pfs_final)
#p_pfs<-surv_pvalue(fit_pfs)$pval

#pick only significant values and plot their graphs
#if(p_os<0.05) {
#png(file= paste0("OS_MUT_SC_", "FAT4", ".png"))
#print(ggsurvplot(
#fit_os, 
#pval = TRUE, 
#risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE)
#)
#dev.off()#}
#if(p_pfs<0.05) {
#png(file= paste0("PFS_MUT_SC_", "FAT4", ".png"))
#print(ggsurvplot(
#fit_pfs, 
#pval = TRUE, 
#risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE)
#)
#dev.off()#}
#}

#creating rna_expr file for .gct and .cls file formats
rna_expr<- read.csv("D:\\COSMIC_DB_FAT4_PC-cell-lines\\CosmicCLP_CompleteGeneExpression.tsv\\CosmicCLP_CompleteGeneExpression.tsv", sep = "\t", header = TRUE, check.names = FALSE)
#names(rna_expr)<- gsub(".{3}$", "", names(rna_expr)) #remove last 3 strings of headers
rna_expr<- pivot_wider(rna_expr, id_cols = "GENE_NAME", names_from = "SAMPLE_ID", values_from = "Z_SCORE")
rna_expr<- as.data.frame(rna_expr)
wt_mut<- table_final[, -c(2:4)]
wt_expr<- filter(wt_mut, STATUS %in% "WT")
mut_expr<- filter(wt_mut, STATUS %in% "M")
wt_vec<- wt_expr$ID_SAMPLE
mut_vec<- mut_expr$ID_SAMPLE
rna_cols<- names(rna_expr)
wt_vec1<- wt_vec[wt_vec %in% rna_cols]
wt_vec1<- as.character(wt_vec1)
mut_vec1<- mut_vec[mut_vec %in% rna_cols]
mut_vec1<- as.character(mut_vec1)
wt_expr_1 <- rna_expr[, wt_vec1]
mut_expr_1 <- rna_expr[ , mut_vec1]
a<- cbind(mut_expr_1, wt_expr_1)
expr_cols<- rna_expr[,c(1)]
expr_data_pc<- cbind(expr_cols, a)
expr_data_pc$Hugo_Symbol[expr_data_pc$Hugo_Symbol==""] <- NA

write.csv(expr_data_pc, file = "expr_data_pc_cosmic_FAT4.csv")


#HEATMAP GENERATION:

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#library(ComplexHeatmap)
#library(dplyr)
#library(tidyr)
#b<- read.table("D:\\paad_tcga_pan_can_atlas_2018\\mol212490-sup-0003-tables3.csv", sep= ",", header = T, check.names = FALSE)
#expr_data_pc<- read.csv("D:\\Expression\\expr_data_pc.csv", header = TRUE, check.names = FALSE)
#names(b)<- b[1,]
#b<- b[-c(1),]
#class(b$`WNT Pathway`)= "numeric"
#b_mod<- select(b, 1:2)
#b_mod<- filter(b_mod, `WNT Pathway` %in% 1)
#gene_vec<- b_mod$`Gene Name`
#gene_expr<- expr_data_pc$Hugo_Sym
#gene_expr<- gene_expr[gene_expr != ""]
#gene_expr1<- gene_expr[gene_expr %in% gene_vec]
#gene_expr_wnt <- expr_data_pc[expr_data_pc$Hugo_Sym %in% gene_vec, ]
#write.csv(gene_expr_wnt, file = "rna_expr_wnt.csv")


#BOX PLOT
library(tidyr)
library(dplyr)
library(ggplot2)
box_plot_table<- read.csv("expr_data_pc_cosmic_FAT4.csv", header = T, check.names = F)
box_plot_table<- box_plot_table[,-c(1)]
df_status <- data.frame(Mutation_Status = c(rep("MUT", 6), rep("WT", 26)))
emt_genes <- read.table("emt_gene_list.csv", sep=",", header=T, check.names=F)
genes<- emt_genes$common_genes
genes<- genes[genes %in% box_plot_table$expr_cols]
box_plot_list <- as.list(genes)
p_val<- data.frame()
#creating empty table with column names
for (i in 1:length(box_plot_list)) {
  emt_gene<-box_plot_list[i]
  print(emt_gene)
  box_plot_table_gene<- box_plot_table[box_plot_table$expr_cols %in% emt_gene,]
  box_plot_gene_longer<- pivot_longer(box_plot_table_gene, cols = c(2:ncol(box_plot_table_gene)), names_to = "GENE_ID", values_to = "Z_SCORE")
  box_plot_gene_longer_1<- cbind(box_plot_gene_longer, df_status)
  #subsetting data for t and wilcox test 
  exp_mut <- subset(box_plot_gene_longer_1, box_plot_gene_longer_1$Mutation_Status =="MUT")$Z_SCORE
  exp_wt <- subset(box_plot_gene_longer_1, box_plot_gene_longer_1$Mutation_Status=="WT")$Z_SCORE
  #p value testing by wilcox and t test
  p<-wilcox.test(exp_mut, exp_wt)$p.val
  p_t<- t.test(exp_mut, exp_wt)$p.val
  if(p_t<0.05){
    png(paste0("EMT_PC_FAT4_BoxPlot_Cosmic_p_t_", emt_gene ,".png"))
    print(bp<- ggplot(box_plot_gene_longer_1, aes(x=Mutation_Status, y=log2(Z_SCORE))) + 
            geom_boxplot(aes(fill=Mutation_Status)) +
            geom_jitter())
    dev.off()}
  
  #creating csv file with storing the values in iteration one after another in the loop
  y<- rbind.data.frame(p_val, cbind.data.frame(emt_gene, p))
  write.table(y, file = "genes_emt_p-value_geneexpr_wt_Cosmic.csv", append = T, col.names = F, row.names = F, sep = ",")
  z<- rbind.data.frame(p_val, cbind.data.frame(wnt_gene, p_t))
  write.table(z, file = "genes_emt_p-value_geneexpr_tt_Cosmic.csv", append = T, col.names = F, row.names = F, sep = ",")
}