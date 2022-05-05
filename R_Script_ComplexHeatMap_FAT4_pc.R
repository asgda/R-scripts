if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
b<- read.table("D:\\paad_tcga_pan_can_atlas_2018\\mol212490-sup-0003-tables3.csv", sep= ",", header = T, check.names = FALSE)
expr_data_pc<- read.csv("D:\\Expression\\expr_data_pc.csv", header = TRUE, check.names = FALSE)
names(b)<- b[1,]
b<- b[-c(1),]
class(b$`WNT Pathway`)= "numeric"
b_mod<- select(b, 1:2)
b_mod<- filter(b_mod, `WNT Pathway` %in% 1)
gene_vec<- b_mod$`Gene Name`
gene_expr<- expr_data_pc$Hugo_Sym
gene_expr<- gene_expr[gene_expr != ""]
gene_expr1<- gene_expr[gene_expr %in% gene_vec]
gene_expr_wnt <- expr_data_pc[expr_data_pc$Hugo_Sym %in% gene_vec, ]
write.csv(gene_expr_wnt, file = "rna_expr_wnt.csv")
#heatmap
rm(list=ls())
options(scipen=999)


mat <- read.table("rna_expr_wnt.csv", sep=",", header=T, check.names=F)
mat = mat[,-c(1,2,4)]
rownames(mat) <- mat[,1]
mat[,1] <- NULL 
cols <- colnames(mat)
cols<-cols[-c(1)]
for(i in 1:nrow(mat)){
  zscore <- (mat[i,] - mean(as.numeric(mat[i,]))) / sd (as.numeric(mat[i,]))
  write.table(zscore, file="zscore_pan_wnt.csv", col.names=F, sep="\t", quote=F, append=T)
}

mat1 <- read.table("D:\\Expression\\zscore_pan_wnt.csv", sep = "\t", header = F)
mat1 <- mat1[,c(2:ncol(mat1))]
colnames(mat1) <- cols
rownames(mat1) <- rownames(mat)
mat1 <- as.matrix(mat1)

df <- data.frame(Condition = c(rep("MUT", 6), rep("WT", 171)))
ha = HeatmapAnnotation(df = df, col = list(Condition = c("MUT" = "#DD292A", "WT" = "#4FBFAD")), show_annotation_name = TRUE)

pdf("Wnt_Heat_Map_FAT4_PC.pdf")
Heatmap(mat1, top_annotation = ha, show_row_names = TRUE, show_column_names=FALSE, heatmap_legend_param = list(title="Z-score"), cluster_columns=T, cluster_rows = T)
dev.off()

#pdf("2.pdf")
#Heatmap(mat1, top_annotation = ha, show_row_names = TRUE, show_column_names=FALSE, heatmap_legend_param = list(title="Z-score"), cluster_columns=FALSE)
#dev.off()

#system("rm x.txt")