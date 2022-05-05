#automate the script
data_patient<- read.delim("C:\\Users\\Avik Sengupta\\Documents\\paad_tcga_pan_can_atlas_2018\\data_clinical_patient.txt", as.is = TRUE)
data_microbiome<- read.delim("C:\\Users\\Avik Sengupta\\Documents\\paad_tcga_pan_can_atlas_2018\\data_microbiome.txt", as.is = TRUE)
library(dplyr)
library(survival)
library(survminer)
library(tidyr)
data_microbiome<- data_microbiome[, c(-2, -3, -4)]
data_micro_final<- pivot_longer(data_microbiome, cols = 2:ncol(data_microbiome), names_to = "Tumor_Sample_Barcode", values_to = "Reads_Num")
#length(unique(data_micro_final$ENTITY_STABLE_ID))
#sort(table(data_micro_final$ENTITY_STABLE_ID), decreasing = TRUE)[1:25] = get the most recurrent genes
#unique(data_micro_final$ENTITY_STABLE_ID)
#get the most recurrent pathogen
my_list <- as.list(unique(data_micro_final$ENTITY_STABLE_ID))
#my_list
for (i in 1:length(my_list)) {
  micro<-my_list[i]
  print(micro)
  data_micro_final$Infection_Status<- ifelse(data_micro_final$ENTITY_STABLE_ID==micro, "Infected", "Not Infected")
  micro_gene<- filter(data_micro_final, ENTITY_STABLE_ID %in% micro)
  micro_gene_uni<- micro_gene %>% distinct(Tumor_Sample_Barcode,.keep_all = TRUE)
  micro_gene_uni$Lower_Quartile_0.25<-quantile(micro_gene_uni$Reads_Num, prob = c(0.25))
  #micro_gene_uni$Lower_Quartile_0<-quantile(micro_gene_uni$Reads_Num, prob = c(0))
  micro_gene_uni$Upper_Quartile_0.75<-quantile(micro_gene_uni$Reads_Num, prob = c(0.75))
  #micro_gene_uni$Upper_Quartile_1.0<-quantile(micro_gene_uni$Reads_Num, prob = c(1))
  micro_gene_uni<- mutate(micro_gene_uni, Reads_Num_status=case_when(Reads_Num<=Lower_Quartile_0.25 ~ "Low Reads"
                                                                       , Reads_Num>=Upper_Quartile_0.75 ~ "High Reads"
                                                                       , TRUE ~ "NIL")
  )
  os_status<- data_patient[, c(1, 30, 31)]
  pfs_status<- data_patient[, c(1, 36, 37)]
  micro_gene_uni$Tumor_Sample_Barcode<- gsub(".{3}$", "", micro_gene_uni$Tumor_Sample_Barcode)
  micro_gene_uni$Tumor_Sample_Barcode<- gsub("[.]", "-", micro_gene_uni$Tumor_Sample_Barcode)
  os_status$OS_STATUS<- substr(os_status$OS_STATUS, 0, 1)
  pfs_status$PFS_STATUS<- substr(pfs_status$PFS_STATUS, 0, 1)
  class(os_status$OS_STATUS)= "numeric"
  class(pfs_status$PFS_STATUS)= "numeric"
  colnames(micro_gene_uni)<- gsub("Tumor_Sample_Barcode", "PATIENT_ID", colnames(micro_gene_uni))
  os_final<- full_join(os_status, micro_gene_uni, "PATIENT_ID")
  pfs_final<- full_join(pfs_status, micro_gene_uni, "PATIENT_ID")
  os_final$Reads_Num_status[is.na(os_final$Reads_Num_status)]<- "NIL"
  pfs_final$Reads_Num_status[is.na(pfs_final$Reads_Num_status)]<- "NIL"
  os_final<- os_final[!(os_final$Reads_Num_status=="NIL"),]
  pfs_final<- pfs_final[!(pfs_final$Reads_Num_status=="NIL"),]
  fit_os<- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Reads_Num_status, data=os_final)
  p_os<-surv_pvalue(fit_os)$pval
  fit_pfs<- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ Reads_Num_status, data=pfs_final)
  p_pfs<-surv_pvalue(fit_pfs)$pval
  if(p_os<0.05) {
  png(file= paste0("OS_MICRO_", micro, ".png"))
  print(ggsurvplot(
    fit_os, 
    pval = TRUE, 
    risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE)
  )
  dev.off()}
  if(p_pfs<0.05) {
  png(file= paste0("PFS_MICRO_", micro, ".png"))
  print(ggsurvplot(
    fit_pfs, 
    pval = TRUE, 
    risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE)
  )
  dev.off()}
}