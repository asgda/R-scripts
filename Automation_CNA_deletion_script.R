#automate the script
data_patient<- read.delim("C:\\Users\\Avik Sengupta\\Documents\\paad_tcga_pan_can_atlas_2018\\data_clinical_patient.txt", as.is = TRUE)
data_cna<- read.delim("C:\\Users\\Avik Sengupta\\Documents\\paad_tcga_pan_can_atlas_2018\\data_cna.txt", as.is = TRUE)
library(dplyr)
library(survival)
library(survminer)
library(tidyr)
data_cna_mod<- pivot_longer(data_cna, cols = 3:ncol(data_cna), names_to = "Tumor_Sample_Barcode", values_to = "CNA_Type")
#sort(table(data_mutation$Hugo_Symbol), decreasing = TRUE)[1:25] -get the most recurrent genes
my_list <- list("KRAS", "TP53", "MUC16", "SMAD4", "OBSCN", "PLEC", "SYNE1", "RYR1", "RYR3", "FAT4", "COL5A1", "DST") 
#my_list
for (i in 1:length(my_list)) {
  gene<-my_list[i]
  print(gene)
  data_cna_mod$deletion_status<- ifelse(data_cna_mod$Hugo_Symbol==gene, "", "No Deletion")
  cna_gene<- filter(data_cna_mod, Hugo_Symbol %in% gene)
  cna_gene_uni<- cna_gene %>% distinct(Tumor_Sample_Barcode,.keep_all = TRUE)
  cna_gene_uni$deletion_status<- ifelse(cna_gene_uni$CNA_Type==-2, "Deletion", "No Deletion")
  os_status<- data_patient[, c(1, 30, 31)]
  pfs_status<- data_patient[, c(1, 36, 37)]
  cna_gene_uni$Tumor_Sample_Barcode<- gsub(".{3}$", "", cna_gene_uni$Tumor_Sample_Barcode)
  cna_gene_uni$Tumor_Sample_Barcode<- gsub("[.]", "-", cna_gene_uni$Tumor_Sample_Barcode)
  os_status$OS_STATUS<- substr(os_status$OS_STATUS, 0, 1)
  pfs_status$PFS_STATUS<- substr(pfs_status$PFS_STATUS, 0, 1)
  class(os_status$OS_STATUS)= "numeric"
  class(pfs_status$PFS_STATUS)= "numeric"
  colnames(cna_gene_uni)<- gsub("Tumor_Sample_Barcode", "PATIENT_ID", colnames(cna_gene_uni))
  os_final<- full_join(os_status, cna_gene_uni, "PATIENT_ID")
  pfs_final<- full_join(pfs_status, cna_gene_uni, "PATIENT_ID")
  os_final$deletion_status[is.na(os_final$deletion_status)]<- "No Deletion"
  pfs_final$deletion_status[is.na(pfs_final$deletion_status)]<- "No Deletion"
  fit_os<- survfit(Surv(OS_MONTHS, OS_STATUS) ~ deletion_status, data=os_final)
  fit_pfs<- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ deletion_status, data=pfs_final)
  png(file= paste0("OS_DEL_", gene, ".png"))
  print(ggsurvplot(
    fit_os, 
    pval = TRUE, 
    risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE)
    )
  dev.off()
  png(file= paste0("PFS_DEL_", gene, ".png"))
  print(ggsurvplot(
    fit_pfs, 
    pval = TRUE, 
    risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE)
    )
 dev.off()
}