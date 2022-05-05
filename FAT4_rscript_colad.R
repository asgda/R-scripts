#automate the script
data_patient<- read.delim("D:\\coadread_tcga_pan_can_atlas_2018\\data_clinical_patient.txt", as.is = TRUE)
data_mutation<- read.delim("D:\\coadread_tcga_pan_can_atlas_2018\\data_mutations.txt", as.is = TRUE)
library(dplyr)
library(survival)
library(survminer)
#sort(table(data_mutation$Hugo_Symbol), decreasing = TRUE)[1:25] -get the most recurrent genes
#my_list <- list("KRAS", "TP53", "MUC16", "SMAD4", "OBSCN", "PLEC", "SYNE1", "RYR1", "RYR3", "FAT4", "COL5A1", "DST") 
#my_list
#for (i in 1:length(my_list)) {
  #gene<-my_list[i]
  #print(gene)
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
  class(pfs_status$PFS_STATUS)= "numeric"
  mut_gene_uni<- mutate(mut_gene_uni, Variant_Classification=case_when(Variant_Classification=="Silent" ~ "WT"
                                                                       , Variant_Classification=="Intron" ~ "WT"
                                                                       , Variant_Classification=="UTR" ~ "WT"
                                                                       , TRUE ~ "M")
  )
  
  colnames(mut_gene_uni)<- gsub("Tumor_Sample_Barcode", "PATIENT_ID", colnames(mut_gene_uni))
  os_final<- full_join(os_status, mut_gene_uni, "PATIENT_ID")
  pfs_final<- full_join(pfs_status, mut_gene_uni, "PATIENT_ID")
  os_final$Variant_Classification[is.na(os_final$Variant_Classification)]<- "WT"
  pfs_final$Variant_Classification[is.na(pfs_final$Variant_Classification)]<- "WT"
  fit_os<- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Variant_Classification, data=os_final)
  p_os<-surv_pvalue(fit_os)$pval
  fit_pfs<- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ Variant_Classification, data=pfs_final)
  p_pfs<-surv_pvalue(fit_pfs)$pval
  #pick only significant values and plot their graphs
  #if(p_os<0.05) {
    png(file= paste0("OS_MUT_CC", "FAT4", ".png"))
    print(ggsurvplot(
      fit_os, 
      pval = TRUE, 
      risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE)
    )
    dev.off()#}
  #if(p_pfs<0.05) {
    png(file= paste0("PFS_MUT_CC", "FAT4", ".png"))
    print(ggsurvplot(
      fit_pfs, 
      pval = TRUE, 
      risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE)
    )
    dev.off()#}
#}