setwd("~/Desktop/PhD_Project_related/GENERATION R DATA RESULTS/Data/Data/ONLY THESE RESULTS TO USE WITH CORRECT COEFFICIENT/RESULTS FROM IRENE/All the models Irene ran")


library(dplyr)

## Model 1
pm25_adj=read.table(file="GenR_M2_pm25_glmFit_18787.csv",header = T,stringsAsFactors = F,sep=" ")
pm10_adj=read.table(file="GenR_M2_pm10_glmFit_18787.csv",header = T,stringsAsFactors = F)
no2_adj=read.table(file="GenR_M2_no2_glmFit_18787.csv",header = T,stringsAsFactors = F,sep=" ")

library(org.Hs.eg.db)

#pm25_adj$Gene=rownames(pm25_adj)
#pm10_adj$Gene=rownames(pm10_adj)
#no2_adj$Gene=rownames(no2_adj)

pm25_adj$p.adjusted=p.adjust(pm25_adj$PValue,method = "BH")
no2_adj$p.adjusted=p.adjust(no2_adj$PValue,method = "BH")
pm10_adj$p.adjusted=p.adjust(pm10_adj$PValue,method = "BH")



pm25_adj$Neg_log_p_val=-log10(pm25_adj$p.adjusted)
pm10_adj$Neg_log_p_val=-log10(pm10_adj$p.adjusted)
no2_adj$Neg_log_p_val=-log10(no2_adj$p.adjusted)

pm25_adj$Gene=rownames(pm25_adj)
pm10_adj$Gene=rownames(pm10_adj)
no2_adj$Gene=rownames(no2_adj)


library(stringr)
pm25_adj[c('Gene', 'Dot')] <- str_split_fixed(pm25_adj$Gene, '\\.', 2)
pm10_adj[c('Gene', 'Dot')] <- str_split_fixed(pm10_adj$Gene, '\\.', 2)
no2_adj[c('Gene', 'Dot')] <- str_split_fixed(no2_adj$Gene, '\\.', 2)



pm25_adj$Entrez <- mapIds(org.Hs.eg.db, pm25_adj$Gene,keytype="ENSEMBL", column="ENTREZID")
pm25_adj$Symbol <- mapIds(org.Hs.eg.db, pm25_adj$Entrez,keytype="ENTREZID", column="SYMBOL")

pm10_adj$Entrez <- mapIds(org.Hs.eg.db, pm10_adj$Gene,keytype="ENSEMBL", column="ENTREZID")
pm10_adj$Symbol <- mapIds(org.Hs.eg.db, pm10_adj$Entrez,keytype="ENTREZID", column="SYMBOL")

no2_adj$Entrez <- mapIds(org.Hs.eg.db, no2_adj$Gene,keytype="ENSEMBL", column="ENTREZID")
no2_adj$Symbol <- mapIds(org.Hs.eg.db, no2_adj$Entrez,keytype="ENTREZID", column="SYMBOL")


pm25_adj$Symbol=as.character(pm25_adj$Symbol)
pm10_adj$Symbol=as.character(pm10_adj$Symbol)
no2_adj$Symbol=as.character(no2_adj$Symbol)


## get most sig probe
pm25_up=subset(pm25_adj,pm25_adj$p.adjusted < 0.25 & pm25_adj$logFC > 0)
pm10_up=subset(pm10_adj,pm10_adj$p.adjusted < 0.05 & pm10_adj$logFC > 0)
no2_up=subset(no2_adj,no2_adj$p.adjusted < 0.25 & no2_adj$logFC > 0)



pm25_dn=subset(pm25_adj,pm25_adj$p.adjusted < 0.25 & pm25_adj$logFC < 0)
pm10_dn=subset(pm10_adj,pm10_adj$p.adjusted < 0.05 & pm10_adj$logFC < 0)
no2_dn=subset(no2_adj,no2_adj$p.adjusted < 0.25 & no2_adj$logFC < 0)





intersect(pm25_up$Gene,no2_up$Gene)
intersect(pm25_dn$Gene,no2_dn$Gene)


pm25_adj=pm25_adj[order(pm25_adj$p.adjusted),]
pm10_adj=pm10_adj[order(pm10_adj$p.adjusted),]
no2_adj=no2_adj[order(no2_adj$p.adjusted),]


#pm25_top_20=pm25_adj[1:20,]
#pm10_top_20=pm10_adj[1:2,]
#no2_top_20=no2_adj[1:1,]



pm25_top_20$Symbol=as.character(pm25_top_20$Symbol)
pm10_top_20$Symbol=as.character(pm10_top_20$Symbol)
#no2_top_20$Symbol=as.character(no2_top_20$Symbol)




#pm25_up_1$Symbol[pm25_up_1$Gene=="ENSG00000226121"]="AHCTF1P1"

pm25_to_highlight=subset(pm25_adj,pm25_adj$Symbol %in% c("GFI1B","TSC22D1","GNG11"))


library(ggrepel)
library(ggplot2)

summary(pm25_adj$logFC)
p_pm25=ggplot(pm25_adj, aes(logFC, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=pm25_adj, aes(x=logFC, y=Neg_log_p_val), colour="grey", size=.5)+
  geom_point(data=pm25_up, aes(x=logFC, y=Neg_log_p_val,label=Symbol), colour="red", size=0.5)+
  geom_point(data=pm25_dn, aes(x=logFC, y=Neg_log_p_val,label=Symbol), colour="blue", size=0.5)+
  geom_text_repel(data = pm25_to_highlight, aes(logFC, Neg_log_p_val,label=Symbol,size=5,fontface="italic"),point.padding = 0.25,nudge_y = 0.1)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))

p1_pm25 <- p_pm25 +labs(x= expression (`Beta coefficient for`~PM[2.5]))+ylab("-Log10(FDR)")+xlim(-0.9,0.75)+geom_hline(yintercept = -log10(0.25), linetype="dotted", color = "grey", size=1.5)+ggtitle("")+geom_vline(xintercept = 0, linetype="dotted", color = "grey", size=1.5)+ylim(0,1.3)+theme(plot.title = element_text(hjust = 0.5,size = 10))
p1_pm25


tiff(filename="PM2.5 volcano plot for presentation.tiff",width = 3000,height = 2000,res = 300)
p1_pm25
dev.off()



summary(pm10_adj$logFC)
p_pm10=ggplot(pm10_adj, aes(logFC, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=pm10_adj, aes(x=logFC, y=Neg_log_p_val), colour="grey", size=2)+
  geom_point(data=pm10_up, aes(x=logFC, y=Neg_log_p_val,label=Symbol), colour="red", size=2)+
  geom_point(data=pm10_dn, aes(x=logFC, y=Neg_log_p_val,label=Symbol), colour="blue", size=2)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))

p1_pm10 <- p_pm10 +xlab("Beta coefficient for PM10")+ylab("-Log10(p value)")+xlim(-0.75,0.75)+geom_hline(yintercept = 1.30103, linetype="dotted", color = "grey", size=1.5)+ggtitle("Volcano plot for association of PM10 with gene expression (with FDR < 0.05)")+geom_vline(xintercept = 0, linetype="dotted", color = "grey", size=1.5)+ylim(0,3)+theme(plot.title = element_text(hjust = 0.5,size=10))+geom_hline(yintercept = 2, linetype="dotted", color = "grey", size=1.5)
p1_pm10




summary(no2_adj$logFC)
p_no2=ggplot(no2_adj, aes(logFC, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=no2_adj, aes(x=logFC, y=Neg_log_p_val), colour="grey", size=2)+
  geom_point(data=no2_up, aes(x=logFC, y=Neg_log_p_val,label=Symbol), colour="red", size=2)+
  geom_point(data=no2_dn, aes(x=logFC, y=Neg_log_p_val,label=Symbol), colour="blue", size=2)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))

p1_no2 <- p_no2 +xlab("Beta coefficient for NO2")+ylab("-Log10(p value)")+xlim(-0.1,0.08)+geom_hline(yintercept = 0.60206, linetype="dotted", color = "grey", size=1.5)+ggtitle("Volcano plot for association of NO2 with gene expression expression(gender, age and batch adjusted) with FDR < 0.05")+geom_vline(xintercept = 0, linetype="dotted", color = "grey", size=1.5)+ylim(0,10)+theme(plot.title = element_text(hjust = 0.5,size = 10))
p1_no2


















