#automate the script
data_patient<- read.delim("D:\\paad_qcmg_uq_2016\\paad_qcmg_uq_2016\\data_clinical_patient.txt", as.is = TRUE, check.names = F)
data_mutation<- read.delim("D:\\paad_qcmg_uq_2016\\paad_qcmg_uq_2016\\data_mutations_mskcc.txt", as.is = TRUE, check.names = F)
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
data_patient<- data_patient[-c(1:3),] #remove the first 3 rows not needed, only keep the needed header as a row value
names(data_patient)<- data_patient[1,] #replace the header with value in the 1st row
data_patient<- data_patient[-c(1),] #after renaming the header, remove the 1st row that had column names
mut<- select(data_mutation, Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)
mut$mutation_status<- ifelse(mut$Hugo_Symbol=="FAT4", "M", "WT")
mut_gene<- filter(mut, Hugo_Symbol %in% "FAT4")
mut_gene_uni<- mut_gene %>% distinct(Tumor_Sample_Barcode,.keep_all = TRUE)
abc_status<- data_patient[, c(1, 10)]
#pfs_status<- data_patient[, c(1, 36, 37)]
#mut_gene_uni$Tumor_Sample_Barcode<- gsub(".{3}$", "", mut_gene_uni$Tumor_Sample_Barcode)
#os_status$OS_STATUS<- substr(os_status$OS_STATUS, 0, 1)
#pfs_status$PFS_STATUS<- substr(pfs_status$PFS_STATUS, 0, 1)
#class(os_status$OS_STATUS)= "numeric"
#class(os_status$OS_MONTHS)= "numeric"
#class(pfs_status$PFS_STATUS)= "numeric"
#class(pfs_status$PFS_MONTHS)= "numeric"
mut_gene_uni<- mutate(mut_gene_uni, Status=case_when(Variant_Classification=="Silent" ~ "WT"
                                                                     , Variant_Classification=="Intron" ~ "WT"
                                                                     , Variant_Classification=="UTR" ~ "WT"
                                                                     , TRUE ~ "M")
)

colnames(mut_gene_uni)<- gsub("Tumor_Sample_Barcode", "PATIENT_ID", colnames(mut_gene_uni))
table_final<- full_join(abc_status, mut_gene_uni, "PATIENT_ID")
table_final<- table_final[, -c(2,5)]
#pfs_final<- full_join(pfs_status, mut_gene_uni, "PATIENT_ID")
#pfs_final<- pfs_final[, -c(4,6)]
table_final$Status[is.na(table_final$Status)]<- "WT"
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
rna_expr<- read.delim("D:\\paad_qcmg_uq_2016\\paad_qcmg_uq_2016\\data_RNA_Seq_v2_expression_median.txt", header = TRUE, check.names = FALSE)
#names(rna_expr)<- gsub(".{3}$", "", names(rna_expr)) #remove last 3 strings of headers

wt_mut<- table_final[, -c(2,3)]
wt_expr<- filter(wt_mut, Status %in% "WT")
mut_expr<- filter(wt_mut, Status %in% "M")
wt_vec<- wt_expr$PATIENT_ID
mut_vec<- mut_expr$PATIENT_ID
rna_cols<- names(rna_expr)
wt_vec1<- wt_vec[wt_vec %in% rna_cols]
mut_vec1<- mut_vec[mut_vec %in% rna_cols]
wt_expr_1 <- rna_expr[, wt_vec1]
mut_expr_1 <- rna_expr[, mut_vec1]
a<- cbind(mut_expr_1, wt_expr_1)
expr_cols<- rna_expr[,c(1,2)]
expr_data_pc<- cbind(expr_cols, a)
expr_data_pc$Hugo_Symbol[expr_data_pc$Hugo_Symbol==""] <- NA

write.csv(expr_data_pc, file = "expr_data_pc_bailey_FAT4.csv")


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
box_plot_table<- read.csv("expr_data_pc_bailey_FAT4.csv", header = T, check.names = F)
box_plot_table<- box_plot_table[,-c(1,3)]
df_status <- data.frame(Mutation_Status = c(rep("MUT", 5), rep("WT", 91)))
wnt_genes <- read.table("rna_expr_wnt.csv", sep=",", header=T, check.names=F)
genes<- wnt_genes$Hugo_Sym
genes<- genes[genes %in% box_plot_table$Hugo_Symbol]
box_plot_list <- as.list(genes)
p_val<- data.frame()
#creating empty table with column names
for (i in 1:length(box_plot_list)) {
  wnt_gene<-box_plot_list[i]
  print(wnt_gene)
  box_plot_table_gene<- box_plot_table[box_plot_table$Hugo_Symbol %in% wnt_gene,]
  box_plot_gene_longer<- pivot_longer(box_plot_table_gene, cols = c(2:ncol(box_plot_table_gene)), names_to = "Tumor_Sample_Barcode", values_to = "Expr_Data")
  box_plot_gene_longer_1<- cbind(box_plot_gene_longer, df_status)
  #subsetting data for t and wilcox test 
  exp_mut <- subset(box_plot_gene_longer_1, box_plot_gene_longer_1$Mutation_Status =="MUT")$Expr_Data
  exp_wt <- subset(box_plot_gene_longer_1, box_plot_gene_longer_1$Mutation_Status=="WT")$Expr_Data
  #p value testing by wilcox and t test
  p<-wilcox.test(exp_mut, exp_wt)$p.val
  p_t<- t.test(exp_mut, exp_wt)$p.val
  if(p<0.05){
    png(paste0("WNT_PC_FAT4_BoxPlot_Bailey", wnt_gene ,".png"))
    print(bp<- ggplot(box_plot_gene_longer_1, aes(x=Mutation_Status, y=log2(Expr_Data))) + 
            geom_boxplot(aes(fill=Mutation_Status)) +
            geom_jitter())
    dev.off()}
  
  #creating csv file with storing the values in iteration one after another in the loop
  y<- rbind.data.frame(p_val, cbind.data.frame(wnt_gene, p))
  write.table(y, file = "genes_wnt_p-value_geneexpr_wt_bailey.csv", append = T, col.names = F, row.names = F, sep = ",")
  z<- rbind.data.frame(p_val, cbind.data.frame(wnt_gene, p_t))
  write.table(z, file = "genes_wnt_p-value_geneexpr_tt_bailey.csv", append = T, col.names = F, row.names = F, sep = ",")
}