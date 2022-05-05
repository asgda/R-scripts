library(tidyr)
library(dplyr)
library(ggplot2)
box_plot_table<- read.csv("expr_data_pc.csv", header = T, check.names = F)
box_plot_table<- box_plot_table[,-c(1,3)]
df_status <- data.frame(Status = c(rep("MUT", 6), rep("WT", 171)))
emt_genes <- read.table("emt_gene_list.csv", sep=",", header=T, check.names=F)
genes<- emt_genes$common_genes
bc<-box_plot_table[box_plot_table$Hugo_Sym %in% genes,]
genes<- bc$Hugo_Sym
box_plot_list <- as.list(genes)
p_val<- data.frame()
#creating empty table with column names
for (i in 1:length(box_plot_list)) {
  emt_gene<-box_plot_list[i]
  print(emt_gene)
  box_plot_table_gene<- box_plot_table[box_plot_table$Hugo_Sym %in% emt_gene,]
  box_plot_gene_longer<- pivot_longer(box_plot_table_gene, cols = c(2:ncol(box_plot_table_gene)), names_to = "Tumor_Sample_Barcode", values_to = "Expr_data")
  box_plot_gene_longer_1<- cbind(box_plot_gene_longer, df_status)
  #subsetting data for t and wilcox test
  exp_mut <- subset(box_plot_gene_longer_1, box_plot_gene_longer_1$Status=="MUT")$Expr_data
  exp_wt <- subset(box_plot_gene_longer_1, box_plot_gene_longer_1$Status=="WT")$Expr_data
  #p value testing by wilcox and t test
  p<-wilcox.test(exp_mut, exp_wt)$p.val
  p_t<- t.test(exp_mut, exp_wt)$p.val
  if(p<0.05){
    png(paste0("EMT_PC_FAT4_BoxPlot_", emt_gene ,".png"))
    print(bp<- ggplot(box_plot_gene_longer_1, aes(x=Status, y=log2(Expr_data))) + 
            geom_boxplot(aes(fill=Status)))
    dev.off()}
   
  
################################################## 
  #creating csv file with storing the values in iteration one after another in the loop
  y<- rbind.data.frame(p_val, cbind.data.frame(emt_gene, p))
  write.table(y, file = "genes_emt_p-value_geneexpr_wt.csv", append = T, col.names = F, row.names = F, sep = ",")
  z<- rbind.data.frame(p_val, cbind.data.frame(emt_gene, p_t))
  write.table(z, file = "genes_emt_p-value_geneexpr_tt.csv", append = T, col.names = F, row.names = F, sep = ",")
}