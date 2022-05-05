#automate the script
setwd("D:\\Survival Curves- Pancreatic Cancer")
data_patient<- read.delim("D:\\paad_tcga_pan_can_atlas_2018\\data_clinical_patient.txt", as.is = TRUE, check.names = F)
data_mutation<- read.delim("D:\\paad_tcga_pan_can_atlas_2018\\data_mutations.txt", as.is = TRUE, check.names = F)
library(dplyr)
library(survival)
library(survminer)
#sort(table(data_mutation$Hugo_Symbol), decreasing = TRUE)[1:25] -get the most recurrent genes
my_list1 <- list("KRAS", "TP53", "MUC16", "SMAD4", "OBSCN", "PLEC", "SYNE1", "RYR1", "RYR3", "COL5A1", "DST", "FAT4") 
my_list2 <- list("KRAS", "TP53", "MUC16", "SMAD4", "OBSCN", "PLEC", "SYNE1", "RYR1", "RYR3", "COL5A1", "DST", "FAT4") 
#my_list
for (i in 1:length(my_list1)) {
  for (j in 1:length(my_list2)) {
    gene2<- my_list2[j]
    gene1<- my_list1[i]
    print(paste0(gene1, "_", gene2))
    mut<- select(data_mutation, Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)
    mut1<- mut[, c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode")]
    mut1$mutation_status<- ifelse(mut$Hugo_Symbol==gene1, "M", "WT")
    mut_gene1<- filter(mut1, Hugo_Symbol %in% gene1)
    mut_gene_uni1<- mut_gene1 %>% distinct(Tumor_Sample_Barcode,.keep_all = TRUE)
    mut_gene_uni1$Tumor_Sample_Barcode<- gsub(".{3}$", "", mut_gene_uni1$Tumor_Sample_Barcode)
    mut_gene_uni1<- mutate(mut_gene_uni1, Variant_Classification=case_when(Variant_Classification=="Silent" ~ "WT"
                                                                           , Variant_Classification=="Intron" ~ "WT"
                                                                           , Variant_Classification=="UTR" ~ "WT"
                                                                           , TRUE ~ "M")
    )
    mut_gene_uni_samples1<- mut_gene_uni1$Tumor_Sample_Barcode
    colnames(mut_gene_uni1)<- gsub("Tumor_Sample_Barcode", "PATIENT_ID", colnames(mut_gene_uni1))
    mut_gene_uni_samples1<- mut_gene_uni1$PATIENT_ID
    mut2<- mut[, c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode")]
    mut2$mutation_status<- ifelse(mut$Hugo_Symbol==gene2, "M", "WT")
    mut_gene2<- filter(mut2, Hugo_Symbol %in% gene2)
    mut_gene_uni2<- mut_gene2 %>% distinct(Tumor_Sample_Barcode,.keep_all = TRUE)
    mut_gene_uni2$Tumor_Sample_Barcode<- gsub(".{3}$", "", mut_gene_uni2$Tumor_Sample_Barcode)
    mut_gene_uni2<- mutate(mut_gene_uni2, Variant_Classification=case_when(Variant_Classification=="Silent" ~ "WT"
                                                                           , Variant_Classification=="Intron" ~ "WT"
                                                                           , Variant_Classification=="UTR" ~ "WT"
                                                                           , TRUE ~ "M")
    )
    mut_gene_uni_samples2<- mut_gene_uni2$Tumor_Sample_Barcode
    colnames(mut_gene_uni2)<- gsub("Tumor_Sample_Barcode", "PATIENT_ID", colnames(mut_gene_uni2))
    mut_gene_uni_samples2<- mut_gene_uni2$PATIENT_ID
    #a<- intersect(mut_gene_uni2$PATIENT_ID, mut_gene_uni1$PATIENT_ID)
    mut_gene_uni3<- merge(mut_gene_uni1, mut_gene_uni2, by = "PATIENT_ID")
    mut_gene_uni3$Mutation<- ifelse(mut_gene_uni3$PATIENT_ID %in% mut_gene_uni2$PATIENT_ID && mut_gene_uni3$PATIENT_ID %in% mut_gene_uni1$PATIENT_ID, "M", "WT")
    mut_gene_uni3<- mut_gene_uni3[, -c(2:7)]
    os_status<- data_patient[, c(1, 30, 31)]
    pfs_status<- data_patient[, c(1, 36, 37)]
    os_status$OS_STATUS<- substr(os_status$OS_STATUS, 0, 1)
    pfs_status$PFS_STATUS<- substr(pfs_status$PFS_STATUS, 0, 1)
    class(os_status$OS_STATUS)= "numeric"
    class(pfs_status$PFS_STATUS)= "numeric"
    os_final<- full_join(os_status, mut_gene_uni3, "PATIENT_ID")
    pfs_final<- full_join(pfs_status, mut_gene_uni3, "PATIENT_ID")
    os_final$Mutation[is.na(os_final$Mutation)]<- "WT"
    pfs_final$Mutation[is.na(pfs_final$Mutation)]<- "WT"
    fit_os<- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Mutation, data=os_final)
    p_os<-surv_pvalue(fit_os)$pval
    fit_pfs<- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ Mutation, data=pfs_final)
    p_pfs<-surv_pvalue(fit_pfs)$pval
    if(p_os<0.05){
      png(file= paste0("OS_MUT_PC_", gene1, "_", gene2, ".png"))
      print(ggsurvplot(
        fit_os, 
        pval = TRUE, 
        risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE)
      )
      dev.off()}
    if(p_pfs<0.05){
      png(file= paste0("PFS_MUT_", gene1,"_", gene2, ".png"))
      print(ggsurvplot(
        fit_pfs, 
        pval = TRUE, 
        risk.table = TRUE, risk.table.col = "strata", conf.int = TRUE)
      )
      dev.off()}
  }
}