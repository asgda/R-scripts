setwd("D:\\Expression")
abc<- read.csv("C:\\Users\\Avik Sengupta\\Downloads\\prediction table 12_04_2022.csv", header = T, check.names = F)
def<- read.csv("C:\\Users\\Avik Sengupta\\Downloads\\Samples 295.csv", header = T, check.names = F)
colnames(abc)[4] <- "Status"
colnames(def)[4] <- "Status"
library(ggplot2)
library(tidyr)
#install.packages("ggalluvial")
library(ggalluvial)
Sample_No<- c(1:381)
Sample<- c(1:294)
df<- data.frame(Sample_No)
df1<- data.frame(Sample)
abc<- cbind(df, abc)
def<- cbind(df1, def)
abc$Class<-sub("^","A_",abc$Class)
abc$Class_MLP<- sub("^", "P_", abc$Class_MLP)
def$Class<-sub("^","A_",def$Class)
def$Class_MLP<- sub("^", "P_", def$Class_MLP)
#data("majors")

png("alluvial_dlbcl_set1.png")
ggplot(as.data.frame(abc),
       aes(y = Sample_No, axis1 = Class, axis2 = Class_MLP)) +
  geom_alluvium(aes(fill = Status), curve_type = "sigmoid", width = 1/12) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Class", "Class_MLP"), expand = c(.15, .15)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  ggtitle("Matches and Mismatches by Actual and Predicted data", "DLBCL")+
  scale_fill_viridis_d()+
  theme_bw()+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
dev.off()

png("alluvial_dlbcl_set2.png")
ggplot(as.data.frame(def),
       aes(y = Sample, axis1 = Class, axis2 = Class_MLP)) +
  geom_alluvium(aes(fill = Status), curve_type = "sigmoid", width = 1/12) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Class", "Class_MLP"), expand = c(.15, .15)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  ggtitle("Matches and Mismatches by Actual and Predicted data", "DLBCL")+
  scale_fill_viridis_d()+
  theme_bw()+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
dev.off()