Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

###加载所需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(ggrepel)
library(reshape2)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

setwd("D:\\Work\\grant\\")
scRNA_harmony=readRDS("HM28_scRNA_harmony_monocytes.anno.rds")
Idents(scRNA_harmony)="group.a"
Idents(scRNA_harmony)="anno"
mono=subset(scRNA_harmony,idents =c("Monocytes-like-CCL3L1+","Monocytes-like-RFLNB+",
                                    "Intermediate monocytes-like"))

degdf2 <- FindMarkers(mono,ident.1 = "OR",ident.2 = "SOR", logfc.threshold = 0.01,group.by = "group.a")

degdf2$symbol <- rownames(degdf2)
logFC_t=0
P.Value_t = 5e-2
degdf2$change = ifelse(degdf2$p_val_adj < P.Value_t & degdf2$avg_log2FC < 0,"down",
                       ifelse(degdf2$p_val_adj < P.Value_t & degdf2$avg_log2FC > 0,"up","stable"))
p=ggplot(degdf2, aes(avg_log2FC,  -log10(p_val_adj))) +
  geom_point(alpha=0.4, size=2.8, aes(color=change)) +
  ylab("-log10(P_val_adj)")+xlab("avg_log2FC")+
  scale_color_manual(values=c("blue", "grey","#FF0000"))+
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  #geom_vline(xintercept = -log10(logFC_t),lty=4,col="black",lwd=0.8) +
  #geom_vline(xintercept = log10(logFC_t),lty=4,col="black",lwd=0.8) +
  theme_bw()

write.csv(degdf2, "HM28_mono.3clusters_OR_vs_SOR.csv", quote = F)

p+theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))



markers=c(#"COMMD8","C1orf216","DPY19L4","HGH1","FKBP5",
  "SHMT1","NAGA","FAHD2A","CD3EAP","HCLS1","ELP3","CFAP97","TULP3","ZBTB16","C3AR1"
)#OR组Monocytes-like-CCL3L1+细胞亚群中低表达的基因有较好生存期
markers=c("IGHD","CARD19","RPS29","MNDA","HIF1A","MZB1","MT-ND3","IFNGR2","JARID2","MPEG1"
          
)#OR组在"Monocytes-like-CCL3L1+","Monocytes-like-RFLNB+","Intermediate monocytes-like"细胞亚群中表达的基因有较好生存期

data_selected <- degdf2[markers,]
p + geom_label_repel(data=data_selected,
                     aes(label=rownames(data_selected)))+
  geom_point(size = 5, shape =21, data = data_selected,colour="blue",alpha = 1/10)+theme_classic()