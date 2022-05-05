#automate the script
data_patient<- read.delim("D:\\skcm_tcga_pan_can_atlas_2018\\data_clinical_patient.txt", as.is = TRUE)
data_mutation<- read.delim("D:\\skcm_tcga_pan_can_atlas_2018\\data_mutations.txt", as.is = TRUE)
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
os_status<- data_patient[, c(1, 30, 31)]
pfs_status<- data_patient[, c(1, 36, 37)]
mut_gene_uni$Tumor_Sample_Barcode<- gsub(".{3}$", "", mut_gene_uni$Tumor_Sample_Barcode)
os_status$OS_STATUS<- substr(os_status$OS_STATUS, 0, 1)
pfs_status$PFS_STATUS<- substr(pfs_status$PFS_STATUS, 0, 1)
class(os_status$OS_STATUS)= "numeric"
class(os_status$OS_MONTHS)= "numeric"
class(pfs_status$PFS_STATUS)= "numeric"
class(pfs_status$PFS_MONTHS)= "numeric"
mut_gene_uni<- mutate(mut_gene_uni, Variant_Classification=case_when(Variant_Classification=="Silent" ~ "WT"
                                                                     , Variant_Classification=="Intron" ~ "WT"
                                                                     , Variant_Classification=="UTR" ~ "WT"
                                                                     , TRUE ~ "M")
)

colnames(mut_gene_uni)<- gsub("Tumor_Sample_Barcode", "PATIENT_ID", colnames(mut_gene_uni))
os_final<- full_join(os_status, mut_gene_uni, "PATIENT_ID")
os_final<- os_final[, -c(4,6)]
pfs_final<- full_join(pfs_status, mut_gene_uni, "PATIENT_ID")
pfs_final<- pfs_final[, -c(4,6)]
os_final$Variant_Classification[is.na(os_final$Variant_Classification)]<- "WT"
pfs_final$Variant_Classification[is.na(pfs_final$Variant_Classification)]<- "WT"
fit_os<- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Variant_Classification, data=os_final)
p_os<-surv_pvalue(fit_os)$pval
fit_pfs<- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ Variant_Classification, data=pfs_final)
p_pfs<-surv_pvalue(fit_pfs)$pval

#pick only significant values and plot their graphs
#if(p_os<0.05) {
png(file= paste0("OS_MUT_SC_", "FAT4", ".png"))
print(ggsurvplot(
  fit_os, 
  pval = TRUE, 
  risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE)
)
dev.off()#}
#if(p_pfs<0.05) {
png(file= paste0("PFS_MUT_SC_", "FAT4", ".png"))
print(ggsurvplot(
  fit_pfs, 
  pval = TRUE, 
  risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE)
)
dev.off()#}
#}

#creating rna_expr file for .gct and .cls file formats
rna_expr<- read.delim("D:\\skcm_tcga_pan_can_atlas_2018\\data_mrna_seq_v2_rsem.txt", header = TRUE, check.names = FALSE)
names(rna_expr)<- gsub(".{3}$", "", names(rna_expr)) #remove last 3 strings of headers

wt_mut<- os_final[, -c(2,3)]
wt_expr<- filter(wt_mut, Variant_Classification %in% "WT")
mut_expr<- filter(wt_mut, Variant_Classification %in% "M")
wt_vec<- wt_expr$PATIENT_ID
mut_vec<- mut_expr$PATIENT_ID
rna_cols<- names(rna_expr)
wt_vec1<- wt_vec[wt_vec %in% rna_cols]
mut_vec1<- mut_vec[mut_vec %in% rna_cols]
wt_expr_1 <- rna_expr[, wt_vec1]
mut_expr_1 <- rna_expr[, mut_vec1]
a<- cbind(mut_expr_1, wt_expr_1)
expr_cols<- rna_expr[,c(1,2)]
expr_data_sc<- cbind(expr_cols, a)
expr_data_sc$Hugo_Sym[expr_data_sc$Hugo_Sym==""] <- NA

write.csv(expr_data_sc, file = "expr_data_sc.csv")


#HEATMAP GENERATION
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
b<- read.table("D:\\paad_tcga_pan_can_atlas_2018\\mol212490-sup-0003-tables3.csv", sep= ",", header = T, check.names = FALSE)
expr_data_sc<- read.csv("D:\\Expression\\expr_data_sc.csv", header = TRUE, check.names = FALSE)
names(b)<- b[1,]
b<- b[-c(1),]
class(b$`WNT Pathway`)= "numeric"
b_mod<- select(b, 1:2)
b_mod<- filter(b_mod, `WNT Pathway` %in% 1)
gene_vec<- b_mod$`Gene Name`
gene_expr<- expr_data_sc$Hugo_Sym
gene_expr<- gene_expr[gene_expr != ""]
gene_expr1<- gene_expr[gene_expr %in% gene_vec]
gene_expr_wnt <- expr_data_sc[expr_data_sc$Hugo_Sym %in% gene_vec, ]
write.csv(gene_expr_wnt, file = "rna_expr_wnt_sc.csv")
#heatmap
#rm(list=ls())
options(scipen=999)


mat <- read.table("rna_expr_wnt_sc.csv", sep=",", header=T, check.names=F)
mat = mat[,-c(1,2,4)]
rownames(mat) <- mat[,1]
mat[,1] <- NULL 
cols <- colnames(mat)

for(i in 1:nrow(mat)){
  zscore <- (mat[i,] - mean(as.numeric(mat[i,]))) / sd (as.numeric(mat[i,]))
  write.table(zscore, file="zscore_pan_wnt_sc.csv", col.names=F, sep="\t", quote=F, append=T)
}

mat1 <- read.table("D:\\Expression\\zscore_pan_wnt_sc.csv", sep = "\t", header = F)
mat1 <- mat1[,c(2:ncol(mat1))]
colnames(mat1) <- cols
rownames(mat1) <- rownames(mat)
mat1 <- as.matrix(mat1)

df <- data.frame(Condition = c(rep("MUT", 153), rep("WT", 288)))
ha = HeatmapAnnotation(df = df, col = list(Condition = c("MUT" = "#DD292A", "WT" = "#4FBFAD")), show_annotation_name = TRUE)

pdf("Wnt_Heat_Map_FAT4_SC.pdf")
Heatmap(mat1, top_annotation = ha, show_row_names = TRUE, show_column_names=FALSE, heatmap_legend_param = list(title="Z-score"), cluster_columns=T, cluster_rows = T)
dev.off()

#pdf("2.pdf")
#Heatmap(mat1, top_annotation = ha, show_row_names = TRUE, show_column_names=FALSE, heatmap_legend_param = list(title="Z-score"), cluster_columns=FALSE)
#dev.off()

#system("rm x.txt")