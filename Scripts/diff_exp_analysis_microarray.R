setwd("~/Desktop/PhD_Project_related/DATA_FROM_ALSPAC/NEW_DATA_13TH_MAY_2022/RESULTS NOT ADJUSTED FOR PC USING BRYOIS/FINAL RESULTS TO USE WITH HIERARCHIAL CLUSTERING FOR SAMPLE KEEPING AND COMPLETE CASE ANALYSIS/RESULTS AFTER REMOVAL OF CHR X AND Y")

library(illuminaHumanv3.db)
library(dplyr)
library(foreign)
library(limma)
# Log2 - transformed expression signals were then normalized with quantile normalization of the replicates of each individual followed by quantile normalization across all individuals.
## I did the quantile normalisation across all individuals


bryois=read.csv(file="bryois.csv",header = T,stringsAsFactors = F)



rownames(bryois)=bryois$X

#con = illuminaHumanv3.db
#dbGetQuery(con,"select * ffrom type_dict;")

x <- illuminaHumanv3ENSEMBL
gene_symbol <- illuminaHumanv3GENENAME

mapped_probes <- mappedkeys(x)
mapped_probes_symbol <- mappedkeys(gene_symbol)



xx <- as.list(x[mapped_probes])
#xx=as.data.frame(xx)
#xx=t(xx)
#xx=as.data.frame(xx)
#xx$Illumina_id=rownames(xx)
#colnames(xx)[1]=c("Gene_Name")
table_probe_gene=matrix(0,nrow=28919,ncol=2)
colnames(table_probe_gene)=c("Probe","Gene")
for(index in 1:length(xx)){
  probe_gene=xx[index]
  probe_name=names(probe_gene)
  probe_gene=as.character(probe_gene)
  ## get the gene for a probe
  table_probe_gene[index,1]=probe_name
  table_probe_gene[index,2]=probe_gene
  
}
table_probe_gene=as.data.frame(table_probe_gene)


table_probe_gene=read.csv(file="Probes_mapping_to_one_gene_bryois.csv",header = T,stringsAsFactors = F)
hg18=read.table(file="hg18_ensembl_genes.txt",header = T,stringsAsFactors = F,sep="\t")

`%ni%`=Negate(`%in%`)
chr_keep=c(1:22)
chr_keep=as.character(chr_keep)
hg18=subset(hg18,hg18$Chromosome.Name %in% chr_keep)
table(hg18$Chromosome.Name)

table_probe_gene=subset(table_probe_gene,table_probe_gene$Probe!="")
#table_probe_gene=table_probe_gene[grep("ILMN",table_probe_gene$Probe),]
#table_probe_gene_x=table_probe_gene[-grep(",",table_probe_gene$Gene),]

#library(tidyverse)
#table_probe_gene2 <- table_probe_gene %>% separate(Gene, c('new_col_1','new_col_2'), sep=",")

table_probe_gene_keep=subset(table_probe_gene,table_probe_gene$Gene %in% hg18$Ensembl.Gene.ID)


bryois_probes=bryois$X

length(intersect(bryois_probes,table_probe_gene$Probe)) ## 26944

all_known_probes=intersect(bryois_probes,table_probe_gene_keep$Probe)


bryois_known_probe=bryois[all_known_probes,]
bryois_known_probe$X=NULL

#exp_for_mds=t(bryois_known_probe)
## Removing samples as per plotMDS 
mds.obj.all.948=plotMDS(bryois_known_probe,top=1000,gene.selection = "common")
eigen.vector=mds.obj.all.948$eigen.vectors
rownames(eigen.vector)=colnames(bryois_known_probe)
View(mds.obj.all.948$eigen.vectors)
`%ni%`=Negate(`%in%`)
bryois_known_probe=bryois_known_probe[,colnames(bryois_known_probe) %ni% c("X40633A","X35039A")]
## so plotMDS also suggets that "X40633A","X35039A" are to be removed. 946 now
mds.obj.all.946=plotMDS(bryois_known_probe,top=1000,gene.selection = "common")
eigen.vector.946=mds.obj.all.946$eigen.vectors
rownames(eigen.vector.946)=colnames(bryois_known_probe)
var.all.pCs=apply(eigen.vector.946,2,var)

## MUST GET RID OF THE TWO SAMPLES "X40633A","X35039A" AS THEY MESS UP THE FOLLOWING ANALYSIS
boxplot(bryois_known_probe[,c(1:50)])



## QC on raw data in form of MDS plot
bryois_known_probe_5=t(bryois_known_probe)
bryois_known_probe_dist=dist(bryois_known_probe_5)
bryois_known_probe_dist_cmdscale=cmdscale(bryois_known_probe_dist,eig=TRUE,k=20)
mds.log2.plot=bryois_known_probe_dist_cmdscale$points
mds.log2.plot.2.dim <- data.frame(Sample.Name = rownames(bryois_known_probe_dist_cmdscale$points), Component.1 = bryois_known_probe_dist_cmdscale$points[,1], Component.2 = bryois_known_probe_dist_cmdscale$points[,2])
colnames(mds.log2.plot)=paste("Component",1:20,sep="")

`%ni%`=Negate(`%in%`)
mds.log2.plot=as.data.frame(mds.log2.plot)
mds.log2.plot$Sample.Name=rownames(mds.log2.plot)


library(ggplot2)
p <- ggplot(mds.log2.plot, aes(x = Component1, y = Component2)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle("MDS Plot of 946 samples")
p


bryois_known_probe_5=t(bryois_known_probe_5)
#cor(bryois_known_probe_5[c(1:5),],use="p")

## Using IAC in case WGCNA is used
IAC = cor(bryois_known_probe_5,use="p")
hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
# IAC mean = 0.979

cluster1 = hclust(as.dist(1-IAC),method="average")
plot(cluster1,cex=0.7,labels=dimnames(bryois_known_probe_5)[[2]])


## Outlier detection is done and now quantile normalisation
library(preprocessCore)
bryois_known_probe_5_use=t(bryois_known_probe_5)
bryois_known_probe_norm <- as.data.frame(normalize.quantiles(as.matrix(bryois_known_probe_5_use),keep.names = TRUE))


boxplot(bryois_known_probe_norm[,c(1:50)],main="Random subset of 50 samples of all quantile normalised samples")




library(foreign)
pheno=read.spss('Ulven_B3777_18Feb22.sav',to.data.frame = TRUE)
pheno2=read.spss('Ulven_B3777_14July23_9_yr_old.sav',to.data.frame = TRUE)

library(dplyr)
all_pheno_combined=left_join(pheno,pheno2,by=c("cidB3777"))
table(all_pheno_combined$kz021.x)
table(all_pheno_combined$kz021.y)
prop.table(table(pheno$kz021))*100
pheno2$age_9y= as.numeric(levels(pheno2$age_9y))[pheno2$age_9y]
mean(pheno2$age_9y[!is.na(pheno2$age_9y)])
sd(pheno2$age_9y[!is.na(pheno2$age_9y)])
prop.table(table(pheno2$c645a))*100
table(pheno$c645a)
table(pheno2$c645a)
prop.table(table(pheno2$c645a))*100
table(pheno2$areases_quint_9)
prop.table(table(pheno2$areases_quint_9))*100

pheno2$pm25_9= as.numeric(levels(pheno2$pm25_9))[pheno2$pm25_9]
pheno2$no2_9= as.numeric(levels(pheno2$no2_9))[pheno2$no2_9]

mean(pheno2$no2_9[!is.na(pheno2$no2_9)])
sd(pheno2$no2_9[!is.na(pheno2$no2_9)])


ids=read.csv(file='B3777_connecting_Genetic_data_phenotypic_data.csv',header = T,stringsAsFactors = F)

ids$Exp_id=paste("X",ids$ge_ht12_g1,"A",sep="")

bryois_known_probe_5_use=t(bryois_known_probe_5_use)
ids_for_sample_with_gene_exp=subset(ids,ids$Exp_id %in% colnames(bryois_known_probe))

`%ni%`=Negate(`%in%`)
colnames.pheno2_only=colnames(pheno2)[colnames(pheno2) %ni% colnames(pheno)]
pheno2_use=pheno2[,c("cidB3777","pm10_9","no2_9","pm25_9")]


pheno_final=left_join(pheno,pheno2_use,by=c("cidB3777"))
View(as.data.frame(table(pheno$cidB3777)))

pheno_for_samples_with_gene_exp=left_join(ids_for_sample_with_gene_exp,pheno_final,by=c("cidB3777"))
length(unique(pheno_for_samples_with_gene_exp$cidB3777))
summary(!is.na(pheno_for_samples_with_gene_exp$pm25_9))
summary(pheno_for_samples_with_gene_exp$pm25_9)
no_na_values = pheno_for_samples_with_gene_exp$no2_9 [!is.na(pheno_for_samples_with_gene_exp$no2_9)]
sd(no_na_values)

sum(is.na(pheno_for_samples_with_gene_exp$pm25_9))
sum(is.na(pheno_for_samples_with_gene_exp$no2_9))
sum(is.na(pheno_for_samples_with_gene_exp$pm10_9))
sum(is.na(pheno_for_samples_with_gene_exp$pm25_9))*100/nrow(pheno_for_samples_with_gene_exp)
sum(is.na(pheno_for_samples_with_gene_exp$no2_9))*100/nrow(pheno_for_samples_with_gene_exp)
sum(is.na(pheno_for_samples_with_gene_exp$pm10_9))*100/nrow(pheno_for_samples_with_gene_exp)
sum(is.na(pheno_for_samples_with_gene_exp$age_9y))*100/nrow(pheno_for_samples_with_gene_exp)

## first running with gender, age and then a
air_pol_lipid=pheno_for_samples_with_gene_exp[,c("ge_ht12_g1","cidB3777","Exp_id","no2_9","pm25_9","pm10_9","kz021","age_9y","CHOL_F9","trig_f9","HDL_f9","LDL_f9","VLDL_f9","CRP_f9","IL6_f9")]
rownames(air_pol_lipid)=air_pol_lipid$Exp_id
View(as.data.frame(table(air_pol_lipid$cidB3777)))



cor=stats::cor
air_pol=air_pol_lipid[,c("no2_9","pm25_9","pm10_9")]
lipid=air_pol_lipid[,c("CHOL_F9","trig_f9","HDL_f9","LDL_f9","VLDL_f9")]

lipid$CHOL_F9= as.numeric(levels(lipid$CHOL_F9))[lipid$CHOL_F9]
lipid$trig_f9= as.numeric(levels(lipid$trig_f9))[lipid$trig_f9]
lipid$HDL_f9= as.numeric(levels(lipid$HDL_f9))[lipid$HDL_f9]
lipid$LDL_f9= as.numeric(levels(lipid$LDL_f9))[lipid$LDL_f9]
lipid$VLDL_f9= as.numeric(levels(lipid$VLDL_f9))[lipid$VLDL_f9]

air_pol$no2_9= as.numeric(levels(air_pol$no2_9))[air_pol$no2_9]
air_pol$pm25_9= as.numeric(levels(air_pol$pm25_9))[air_pol$pm25_9]

air_pol_lipid_Cor = cor(air_pol,lipid,use="complete.obs",method="pearson")

no2_tc = cor.test(air_pol$no2_9,lipid$CHOL_F9,method="spearman",alternative = "two.sided") ## spearman sig
pm25_tc = cor.test(air_pol$pm25_9,lipid$CHOL_F9,method="spearman",alternative = "two.sided") ## sig -ve for Pearson and spearman
pm10_tc = cor.test(air_pol$pm10_9,lipid$CHOL_F9,method="spearman",alternative = "two.sided")

no2_trig = cor.test(air_pol$no2_9,lipid$trig_f9,method="spearman",alternative = "two.sided") ## 
pm25_trig = cor.test(air_pol$pm25_9,lipid$trig_f9,method="spearman",alternative = "two.sided") ## 
pm10_trig = cor.test(air_pol$pm10_9,lipid$trig_f9,method="spearman",alternative = "two.sided")

no2_HDL = cor.test(air_pol$no2_9,lipid$HDL_f9,method="spearman",alternative = "two.sided") ## 
pm25_HDL = cor.test(air_pol$pm25_9,lipid$HDL_f9,method="spearman",alternative = "two.sided") ## 
pm10_HDL = cor.test(air_pol$pm10_9,lipid$HDL_f9,method="spearman",alternative = "two.sided")

no2_LDL = cor.test(air_pol$no2_9,lipid$LDL_f9,method="spearman",alternative = "two.sided") ## 
pm25_LDL = cor.test(air_pol$pm25_9,lipid$LDL_f9,method="spearman",alternative = "two.sided") ## sig -ve for Pearson and Spearman
pm10_LDL = cor.test(air_pol$pm10_9,lipid$LDL_f9,method="spearman",alternative = "two.sided")



x=as.data.frame(table(pheno_for_samples_with_gene_exp$Exp_id))


`%ni%` = Negate(`%in%`)
#pheno_for_samples_with_gene_exp=subset(pheno_for_samples_with_gene_exp,pheno_for_samples_with_gene_exp$Exp_id %ni% samples_remove)
### Minimally adjusted for gender and smoked during pregnancy or not
## kz021= sex of the subject
## dw042= pre pregnancy BMI
## bestgest = maternal gestational age 
## b663  = tobacco smoked during pregnancy
## c645a = education level
## c800 = race
## areases_quint_preg = Area SES
pheno_for_samples_with_gene_exp$areases_quint_preg=as.factor(pheno_for_samples_with_gene_exp$areases_quint_preg)
pheno_for_samples_with_gene_exp$dw042 = as.numeric(levels(pheno_for_samples_with_gene_exp$dw042))[pheno_for_samples_with_gene_exp$dw042]

pheno_for_samples_with_gene_exp$pm25_9= as.numeric(levels(pheno_for_samples_with_gene_exp$pm25_9))[pheno_for_samples_with_gene_exp$pm25_9]
#pheno_for_samples_with_gene_exp$pm10_9= as.numeric(levels(pheno_for_samples_with_gene_exp$pm10_9))[pheno_for_samples_with_gene_exp$pm10_9]
pheno_for_samples_with_gene_exp$no2_9= as.numeric(levels(pheno_for_samples_with_gene_exp$no2_9))[pheno_for_samples_with_gene_exp$no2_9]
pheno_for_samples_with_gene_exp$age_9y= as.numeric(levels(pheno_for_samples_with_gene_exp$age_9y))[pheno_for_samples_with_gene_exp$age_9y]

areases_quint_9=pheno2[,c("cidB3777","areases_quint_9")]
pheno_for_samples_with_gene_exp=left_join(pheno_for_samples_with_gene_exp,areases_quint_9,by=c("cidB3777"))
## gender_gestational_age_smoking_during_pregnancy_areases
table(pheno_for_samples_with_gene_exp$areases_quint_9)
prop.table(table(pheno_for_samples_with_gene_exp$areases_quint_9))*100
# CSE - 
# Vocational - secondary and higher 
# O level - Ordinary level
# A level - Advanced level
# Degree - highest

pm25=pheno_for_samples_with_gene_exp[,c("Exp_id","pm25_9","age_9y","kz021","areases_quint_9","c645a")] ## 111 missing
pm10=pheno_for_samples_with_gene_exp[,c("Exp_id","pm10_9","age_9y","kz021","areases_quint_9","c645a")] ## 87 missing
no2=pheno_for_samples_with_gene_exp[,c("Exp_id","no2_9","age_9y","kz021","areases_quint_9","c645a")]  ## 111 missing


mean(pheno_for_samples_with_gene_exp$pm25_9[!is.na(pheno_for_samples_with_gene_exp$pm25_9)])
sd(pheno_for_samples_with_gene_exp$pm25_9[!is.na(pheno_for_samples_with_gene_exp$pm25_9)])

mean(pheno_for_samples_with_gene_exp$pm10_9[!is.na(pheno_for_samples_with_gene_exp$pm10_9)])
sd(pheno_for_samples_with_gene_exp$pm10_9[!is.na(pheno_for_samples_with_gene_exp$pm10_9)])

mean(pheno_for_samples_with_gene_exp$no2_9[!is.na(pheno_for_samples_with_gene_exp$no2_9)])
sd(pheno_for_samples_with_gene_exp$no2_9[!is.na(pheno_for_samples_with_gene_exp$no2_9)])

table(pheno_for_samples_with_gene_exp$kz021[!is.na(pheno_for_samples_with_gene_exp$no2_9)])
table(pheno_for_samples_with_gene_exp$kz021[!is.na(pheno_for_samples_with_gene_exp$pm10_9)])
table(pheno_for_samples_with_gene_exp$kz021[!is.na(pheno_for_samples_with_gene_exp$pm25_9)])


## mean and sd for gesttational age for those with air pol data
## mean and sd for BMI for those with air pol data
x=pheno_for_samples_with_gene_exp$dw042[!is.na(pheno_for_samples_with_gene_exp$no2_preg)]
mean(x[!is.na(x)])
sd(x[!is.na(x)])

x=pheno_for_samples_with_gene_exp$dw042[!is.na(pheno_for_samples_with_gene_exp$pm25_preg)]
mean(x[!is.na(x)])
sd(x[!is.na(x)])

x=pheno_for_samples_with_gene_exp$dw042[!is.na(pheno_for_samples_with_gene_exp$pm10_preg)]
mean(x[!is.na(x)])
sd(x[!is.na(x)])


#pheno_for_samples_with_gene_exp$b862=as.numeric(levels(pheno_for_samples_with_gene_exp$b862))[pheno_for_samples_with_gene_exp$b862]
## physical activity classification for those with air pol data
#x=pheno_for_samples_with_gene_exp$b862[!is.na(pheno_for_samples_with_gene_exp$no2_preg)]
#mean(x[!is.na(x)])
#sd(x[!is.na(x)])

#x=pheno_for_samples_with_gene_exp$b862[!is.na(pheno_for_samples_with_gene_exp$pm25_preg)]
#mean(x[!is.na(x)])
#sd(x[!is.na(x)])

#x=pheno_for_samples_with_gene_exp$b862[!is.na(pheno_for_samples_with_gene_exp$pm10_preg)]
#mean(x[!is.na(x)])
#sd(x[!is.na(x)])


## education mother
## CSE  = 
x=pheno_for_samples_with_gene_exp$c645a[!is.na(pheno_for_samples_with_gene_exp$no2_preg)]
table(x)

x=pheno_for_samples_with_gene_exp$c645a[!is.na(pheno_for_samples_with_gene_exp$pm25_preg)]
table(x)


x=pheno_for_samples_with_gene_exp$c645a[!is.na(pheno_for_samples_with_gene_exp$pm10_preg)]
table(x)






no2_tab_not_na_no2=no2[!is.na(no2$no2_9),]
pm10_tab_not_na_pm10=pm10[!is.na(pm10$pm10_9),]
pm25_tab_not_na_pm25=pm25[!is.na(pm25$pm25_9),]

## Mean for NO2 for those with NO2 value not missing
mean(no2_tab_not_na_no2$no2_9[!is.na(no2_tab_not_na_no2$no2_9)])
sd(no2_tab_not_na_no2$no2_9[!is.na(no2_tab_not_na_no2$no2_9)])

median(no2_tab_not_na_no2$no2_9[!is.na(no2_tab_not_na_no2$no2_9)])


## Mean for PM10 for those with PM10 value not missing
mean(pm10_tab_not_na_pm10$pm10_9[!is.na(pm10_tab_not_na_pm10$pm10_9)])
sd(pm10_tab_not_na_pm10$pm10_9[!is.na(pm10_tab_not_na_pm10$pm10_9)])

median(pm10_tab_not_na_pm10$pm10_9[!is.na(pm10_tab_not_na_pm10$pm10_9)])


mean(pm25_tab_not_na_pm25$pm25_9[!is.na(pm25_tab_not_na_pm25$pm25_9)])
sd(pm25_tab_not_na_pm25$pm25_9[!is.na(pm25_tab_not_na_pm25$pm25_9)])

median(pm25_tab_not_na_pm25$pm25_9[!is.na(pm25_tab_not_na_pm25$pm25_9)])





table(pheno_for_samples_with_gene_exp$kz021[!is.na(pheno_for_samples_with_gene_exp$pm25_9)])
table(pheno_for_samples_with_gene_exp$b663[!is.na(pheno_for_samples_with_gene_exp$pm25_9)])



mean(pheno_for_samples_with_gene_exp$dw042[!is.na(pheno_for_samples_with_gene_exp$dw042)])
sd(pheno_for_samples_with_gene_exp$dw042[!is.na(pheno_for_samples_with_gene_exp$dw042)])







sum(is.na(pm25$pm25_9)) ## 109
sum(is.na(pm10$pm10_9)) ## 149
sum(is.na(no2$no2_9)) ## 109


table(pm25$c645a)
table(pm10$c645a)
table(no2$c645a)

table(pheno_for_samples_with_gene_exp$c645a)
table(pheno_for_samples_with_gene_exp$kz021)
prop.table(table(pheno_for_samples_with_gene_exp$kz021))*100
mean(pheno_for_samples_with_gene_exp$age_9y)
sd(pheno_for_samples_with_gene_exp$age_9y)

prop.table(table(pheno_for_samples_with_gene_exp$c645a))*100
### Missing data analysis
library(finalfit)
table_for_missing_var_plot=pheno_for_samples_with_gene_exp
library(data.table)
#table_for_missing_var_plot=setnames(table_for_missing_var_plot,                              # Apply setnames function
 #                                   c("kz021", "bestgest","b663","c645a","areases_quint_preg","pm25_9","pm10_9","no2_9"),
 #                                   c("Gender", "Gestational_Age","Smoking_during_pregnancy","Mother_Education","Area_SES","PM2.5_9_years","PM10_9_years","NO2_9_years"))

explanatory = c("kz021", "bestgest", 
                "b663", "c645a","areases_quint_preg")
dependent = "pm25_9"
table_for_missing_var_plot %>% 
  missing_pairs(dependent, explanatory,title="Measure of covariates among samples with or without PM2.5 exposure data")


dependent = "pm10_9"
table_for_missing_var_plot %>% 
  missing_pairs(dependent, explanatory,title="Measure of covariates among samples with or without PM10 exposure data")

dependent = "no2_9"
table_for_missing_var_plot %>% 
  missing_pairs(dependent, explanatory,title="Measure of covariates among samples with or without NO2 exposure data")


sum(is.na(pheno_for_samples_with_gene_exp$c645a))

pm25_cmpl=pm25[complete.cases(pm25),] ## 827/946, 87%
pm10_cmpl=pm10[complete.cases(pm10),] ## 715/946, 75.2%
no2_cmpl=no2[complete.cases(no2),] ## 827/946, 87%

#pm25_no_na=pm25[!is.na(pm25$pm25_preg),] ## 708, 81%
#pm10_no_na=pm10[!is.na(pm10$pm10_preg),] ## 667, 82.4%
#no2_no_na=no2[!is.na(no2$no2_preg),] ## 708, 81%


#sum(is.na(pm25_no_na$c645a))
#sum(is.na(pm10_no_na$c645a))
#sum(is.na(no2_no_na$c645a))






## Multiple imputation avoided because address connecting to air polutant and that cannot be imputed

table(pm25_cmpl$kz021)
table(pm25_cmpl$b663)
table(pm25_cmpl$c645a)
summary(pm25_cmpl$bestgest)
mean(pm25_cmpl$bestgest)
sd(pm25_cmpl$bestgest)
mean(pm25_cmpl$dw042)
sd(pm25_cmpl$dw042)

mean(pm25_cmpl$dw042)
sd(pm25_cmpl$dw042)

table(no2_cmpl$kz021)
table(no2_cmpl$b663)

table(pm10_cmpl$kz021)
table(pm10_cmpl$b663)


prop.table(table(no2_cmpl$c645a))*100
prop.table(table(pm25_cmpl$c645a))*100

table(no2_cmpl$c645a)
table(pm25_cmpl$c645a)


table(pheno_for_samples_with_gene_exp$c645a)



pm25_cmpl_x=pm25_cmpl

pm25_cmpl$kz021=droplevels(pm25_cmpl$kz021)
pm25_cmpl$areases_quint_9=droplevels(pm25_cmpl$areases_quint_9)
pm25_cmpl$c645a=droplevels(pm25_cmpl$c645a)
#pm25_cmpl$areases_quint_preg=droplevels(pm25_cmpl$areases_quint_preg)
#pm25_cmpl$dw042=as.numeric(levels(pm25_cmpl$dw042))[pm25_cmpl$dw042]

pm10_cmpl$kz021=droplevels(pm10_cmpl$kz021)
pm10_cmpl$areases_quint_9=droplevels(pm10_cmpl$areases_quint_9)
pm10_cmpl$c645a=droplevels(pm10_cmpl$c645a)
#pm10_cmpl$areases_quint_preg=droplevels(pm10_cmpl$areases_quint_preg)
#pm10_cmpl$dw042=as.numeric(levels(pm10_cmpl$dw042))[pm10_cmpl$dw042]


no2_cmpl$kz021=droplevels(no2_cmpl$kz021)
no2_cmpl$areases_quint_9=droplevels(no2_cmpl$areases_quint_9)
no2_cmpl$c645a=droplevels(no2_cmpl$c645a)
#no2_cmpl$areases_quint_preg=droplevels(no2_cmpl$areases_quint_preg)
#no2_cmpl$dw042=as.numeric(levels(no2_cmpl$dw042))[no2_cmpl$dw042]


#no2_cmpl$areases_quint_preg=droplevels(no2_cmpl$areases_quint_preg)


## Air pollutant does not divide gene exp by quantile

pm25_samples=pm25_cmpl$Exp_id ## 827
pm10_samples=pm10_cmpl$Exp_id ## 715
no2_samples=no2_cmpl$Exp_id ## 827


# Use bryois_known_probe_keep to sleect samples for linear model
#bryois_known_probe_keep=t(bryois_known_probe_keep)





# Gordon Smyth recommends not removing outlier probes
#outlier_no2=read.table(file="Outlier_probes_from_NO2_model.txt",header = F,stringsAsFactors = F)
#outlier_pm25=read.table(file="Outlier_probes_from_PM25_model.txt",header = F,stringsAsFactors = F)
#outlier_pm10=read.table(file="Outlier_probes_from_PM10_model.txt",header = F,stringsAsFactors = F)

#outlier_no2=outlier_no2$V1
#outlier_pm25=outlier_pm25$V1
#outlier_pm10=outlier_pm10$V1



#pm25_cmpl$Exp_id=NULL
#pm10_cmpl$Exp_id=NULL
#no2_cmpl$Exp_id=NULL
rownames(pm25_cmpl)=pm25_cmpl$Exp_id
rownames(pm10_cmpl)=pm10_cmpl$Exp_id
rownames(no2_cmpl)=no2_cmpl$Exp_id


samples_pm25=pm25_cmpl$Exp_id

hist(pm25_cmpl$pm25_9)
hist(pm10_cmpl$pm10_9)
hist(no2_cmpl$no2_9)

mean(pm25_cmpl$pm25_9)
sd(pm25_cmpl$pm25_9)

mean(pm10_cmpl$pm10_9)
sd(pm10_cmpl$pm10_9)

mean(no2_cmpl$no2_9)
sd(no2_cmpl$no2_9)


pm25_cmpl$pm25_9_center=scale(pm25_cmpl$pm25_9,center = TRUE)
pm10_cmpl$pm10_9_center=scale(pm10_cmpl$pm10_9,center = TRUE)
no2_cmpl$no2_9_center=scale(no2_cmpl$no2_9,center = TRUE)

pm25_cmpl$age_9y_centered=scale(pm25_cmpl$age_9y,center = TRUE)
pm10_cmpl$age_9y_centered=scale(pm10_cmpl$age_9y,center = TRUE)
no2_cmpl$age_9y_centered=scale(no2_cmpl$age_9y,center = TRUE)

no2_cmpl$areases_quint_9


#colnames(pm25_cmpl)[7]="pm25_9_centered"
#colnames(pm10_cmpl)[7]="pm10_9_centered"
#colnames(no2_cmpl)[7]="no2_9_centered"

colnames(pm25_cmpl)[8]="age_9y_centered"
colnames(pm10_cmpl)[8]="age_9y_centered"
colnames(no2_cmpl)[8]="age_9y_centered"

colnames(pm25_cmpl)[7]="pm25_9_center"
colnames(pm10_cmpl)[7]="pm10_9_center"
colnames(no2_cmpl)[7]="no2_9_center"


pm25_design_mat=model.matrix( ~ pm25_9_center + kz021 + age_9y_centered ,data=pm25_cmpl)
pm10_design_mat=model.matrix( ~ pm10_9_center + kz021 + age_9y_centered  ,data=pm10_cmpl)
no2_design_mat=model.matrix( ~ no2_9_center + kz021 + age_9y_centered ,data=no2_cmpl)

#rownames(pm25_design_mat)=pm25_cmpl$Exp_id
#rownames(pm10_design_mat)=pm10_cmpl$Exp_id
#rownames(no2_design_mat)=no2_cmpl$Exp_id


boxplot(bryois_known_probe_norm[,1:50])

bryois_pm25=bryois_known_probe_norm[,pm25_cmpl$Exp_id]
bryois_pm10=bryois_known_probe_norm[,pm10_cmpl$Exp_id]
bryois_no2=bryois_known_probe_norm[,no2_cmpl$Exp_id]


par(mfrow = c(1,3))
boxplot(bryois_no2[,1:20],main="NO2")
boxplot(bryois_pm10[,1:20],main="PM10")
boxplot(bryois_pm25[,1:20],main="PM2.5")


no2_g_order=colnames(bryois_no2) # 801
pm25_g_order=colnames(bryois_pm25) ## 801
pm10_g_order=colnames(bryois_pm10) ## 752

no2_c_order=rownames(no2_design_mat) ## 801
pm25_c_order=rownames(pm25_design_mat) ##  801
pm10_c_order=rownames(pm10_design_mat) ## 752


library(preprocessCore)

#perform quantile normalization

bryois_no2_norm <- normalize.quantiles(as.matrix(bryois_no2))
rownames(bryois_no2_norm)=rownames(bryois_no2)
colnames(bryois_no2_norm)=colnames(bryois_no2)

bryois_pm10_norm <- normalize.quantiles(as.matrix(bryois_pm10))
rownames(bryois_pm10_norm)=rownames(bryois_pm10)
colnames(bryois_pm10_norm)=colnames(bryois_pm10)

bryois_pm25_norm <- normalize.quantiles(as.matrix(bryois_pm25))
rownames(bryois_pm25_norm)=rownames(bryois_pm25)
colnames(bryois_pm25_norm)=colnames(bryois_pm25)



boxplot(bryois_no2_norm[,1:50],main="NO2")
boxplot(bryois_pm10_norm[,1:50],main="PM10")
boxplot(bryois_pm25_norm[,1:50],main="PM2.5")

no2_g_order=colnames(bryois_no2_norm) # 814
pm25_g_order=colnames(bryois_pm25_norm) ## 814
pm10_g_order=colnames(bryois_pm10_norm) ## 833

table(no2_g_order==no2_c_order) ## all true
table(pm25_g_order==pm25_c_order) ## all true
table(pm10_g_order==pm10_c_order) ## all true


### THIS IS WHERE THE ACTUAL MODEL FITTING STARTS
table(colnames(bryois_no2_norm)==rownames(no2_design_mat))
no2_lm=lmFit(bryois_no2,no2_design_mat) ## 814
no2_lmfit_eb=eBayes(no2_lm[no2_lm$Amean > 7.5,],trend = TRUE)
summary(decideTests(no2_lmfit_eb)) ## nothing

View(no2_lmfit_eb$coefficients)
test_tab=no2_lmfit_eb$coefficients
test_tab=as.data.frame(test_tab)
test_tab$Probe=rownames(test_tab)

no2.table <- topTable(no2_lmfit_eb,n=14054,coef = "ns(no2_9_center, 3)1") ## check results using MA plot

no2.res.matrix <- residuals(no2_lm, bryois_no2)

plotSA(no2_lmfit_eb,main="Plot of residual standard deviation vs Avg exp for fitted model for NO2",col=c("black","red"))
plotSA(no2_lm,main="Plot of residual standard deviation vs average expression for fitted model using NO2",col=c("black","red"))
volcanoplot(no2_lmfit_eb,coef = "no2_9_center",highlight = 20)



table(colnames(bryois_pm25_norm)==rownames(pm25_design_mat))
pm25_lm=lmFit(bryois_pm25,pm25_design_mat) ## 814
pm25_lmfit_eb=eBayes(pm25_lm[pm25_lm$Amean > 7.5,],trend = TRUE)
View(pm25_lmfit_eb$coefficients)
pm25.table <- topTable(pm25_lmfit_eb, n = 14054,coef = "ns(pm25_9_center, 3)2")
pm25.table$adj.P.Val.Bonferroni=p.adjust(pm25.table$P.Value,method = "bonferroni")
summary(decideTests(pm25_lmfit_eb)) ## 17 for PM25

pm25.res.matrix <- residuals(pm25_lm, bryois_pm25)


plotSA(pm25_lmfit_eb,main="Plot of residual standard deviation vs Avg exp for fitted model using PM25")
plotSA(pm25_lm,main="Plot of residual standard deviation vs Avg exp for fitted model using PM25")

volcanoplot(pm25_lmfit_eb,coef = "ns(pm25_9_center, 3)2")

table(colnames(bryois_pm10_norm)==rownames(pm10_design_mat))
pm10_lm=lmFit(bryois_pm10,pm10_design_mat) ##  
pm10_lmfit_eb=eBayes(pm10_lm[pm10_lm$Amean > 7.5,],trend = TRUE) 
View(pm10_lmfit_eb$coefficients)
pm10.table <- topTable(pm10_lmfit_eb,n = 14086,coef = "ns(pm10_9_center, 3)3")
summary(decideTests(pm10_lmfit_eb)) ## some for PM10

pm10.res.matrix <- residuals(pm10_lm, bryois_pm10)


plotSA(pm10_lmfit_eb,main="Plot of residual standard deviation vs Avg exp for fitted model using PM10")
plotSA(pm10_lm,main="Plot of residual standard deviation vs average expression for fitted model using PM10")

volcanoplot(pm10_lmfit_eb,coef=11)


no2.table$Direction=ifelse(no2.table$logFC>0,"Up","Down")
pm25.table$Direction=ifelse(pm25.table$logFC>0,"Up","Down")
pm10.table$Direction=ifelse(pm10.table$logFC>0,"Up","Down")

no2.table$Check_sig=ifelse(no2.table$adj.P.Val < 0.05,"Sig","Not_sig")
pm25.table$Check_sig=ifelse(pm25.table$adj.P.Val < 0.05,"Sig","Not_sig")
pm10.table$Check_sig=ifelse(pm10.table$adj.P.Val < 0.05,"Sig","Not_sig")

table(no2.table$Check_sig)
table(pm25.table$Check_sig)
table(pm10.table$Check_sig)


no2.table$Probe=rownames(no2.table)
pm25.table$Probe=rownames(pm25.table)
pm10.table$Probe=rownames(pm10.table)



library(snpStats)


pm10_pval_arr <- pm10.table %>% 
  arrange(adj.P.Val) %>%
  mutate(quantile=1:nrow(pm10.table)/nrow(pm10.table))

pm10_pval_arr <- pm10_pval_arr %>% 
  mutate(exp_pval = quantile,
         neglog10_rna = -log10(adj.P.Val),
         neglog10_expected = -log10(exp_pval))


pm25_pval_arr <- pm25.table %>% 
  arrange(adj.P.Val) %>%
  mutate(quantile=1:nrow(pm25.table)/nrow(pm25.table))

pm25_pval_arr <- pm25_pval_arr %>% 
  mutate(exp_pval = quantile,
         neglog10_rna = -log10(adj.P.Val),
         neglog10_expected = -log10(exp_pval))


no2_pval_arr <- no2.table %>% 
  arrange(adj.P.Val) %>%
  mutate(quantile=1:nrow(no2.table)/nrow(no2.table))

no2_pval_arr <- no2_pval_arr %>% 
  mutate(exp_pval = quantile,
         neglog10_rna = -log10(adj.P.Val),
         neglog10_expected = -log10(exp_pval))




library(ggplot2)

p <- ggplot(pm10_pval_arr, aes(x = neglog10_expected, y = neglog10_rna)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("-Log10(Expected p value)") + ylab("-Log10(Observed p value)") + ggtitle("QQ Plot for association of PM10 with gene expression")+geom_abline(slope=1, intercept=0)
p



p <- ggplot(pm25_pval_arr, aes(x = neglog10_expected, y = neglog10_rna)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("-Log10(Expected p value)") + ylab("-Log10(Observed p value)") + ggtitle("QQ Plot for association of PM25 with gene expression")+geom_abline(slope=1, intercept=0)
p


p <- ggplot(no2_pval_arr, aes(x = neglog10_expected, y = neglog10_rna)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("-Log10(Expected p value)") + ylab("-Log10(Observed p value)") + ggtitle("QQ Plot for association of NO2 with gene expression")+geom_abline(slope=1, intercept=0)
p






### Overlap venn diagram
no2_probe_gene=left_join(no2.table,table_probe_gene,by=c("Probe"))
pm10_probe_gene=left_join(pm10.table,table_probe_gene,by=c("Probe"))
pm25_probe_gene=left_join(pm25.table,table_probe_gene,by=c("Probe"))


View(no2_probe_gene)
View(pm10_probe_gene)
View(pm25_probe_gene)


no2_probe_gene_highest_F <- no2_probe_gene %>%
  dplyr::arrange(dplyr::desc(abs(F))) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

pm25_probe_gene_highest_F <- pm25_probe_gene %>%
  dplyr::arrange(dplyr::desc(abs(F))) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

pm10_probe_gene_highest_F <- pm10_probe_gene %>%
  dplyr::arrange(dplyr::desc(abs(F))) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

no2_up=subset(no2_probe_gene,no2_probe_gene$logFC > 0 & no2_probe_gene$adj.P.Val < 0.05)
pm25_up=subset(pm25_probe_gene,pm25_probe_gene$logFC > 0 & pm25_probe_gene$adj.P.Val < 0.05)
pm10_up=subset(pm10_probe_gene,pm10_probe_gene$logFC > 0 & pm10_probe_gene$adj.P.Val < 0.05)

no2_dn=subset(no2_probe_gene,no2_probe_gene$logFC < 0 & no2_probe_gene$adj.P.Val < 0.05)
pm25_dn=subset(pm25_probe_gene,pm25_probe_gene$logFC < 0 & pm25_probe_gene$adj.P.Val < 0.05)
pm10_dn=subset(pm10_probe_gene,pm10_probe_gene$logFC < 0 & pm10_probe_gene$adj.P.Val < 0.05)


library(VennDiagram)


pdf('Overlap_genes_up_NO2_PM25_PM10_9yr_old_Exposure_data.pdf')
p1=venn.diagram(list(`NO2`=no2_up$Gene,`PM10`=pm10_up$Gene,`PM25`=pm25_up$Gene),filename=NULL,cex=2,fill=c("yellow","green","orange"),main.cex=c(1),cat.cex=c(1.5,1.5,1.5),cat.dist=c(.025,.025,.025))
grid.draw(p1)
dev.off()

library(eulerr)
up.overlap <- euler(c("NO2" = 8, "PM10" = 32,"PM25"=40,"NO2&PM25"=36,"PM10&PM25"=4,"NO2&PM10"=36,"NO2&PM10&PM25"=38))
plot(up.overlap, counts = TRUE, font=1, cex=4, alpha=0.5,
     fill=c("yellow", "green","orange"),main="Overlap of positively associated genes after NO2,PM10 & PM25 exposure",labels=list(fontsize=15,cex=2),quantities=list(cex=3))

pdf('Overlap_genes_dn_NO2_PM25_PM10_9yr_old_Exposure_data.pdf')
p1=venn.diagram(list(`NO2`=no2_dn$Gene,`PM10`=pm10_dn$Gene,`PM25`=pm25_dn$Gene),filename=NULL,cex=2,fill=c("yellow","green","orange"),main.cex=c(1),cat.cex=c(1.5,1.5,1.5),cat.dist=c(.025,.025,.025))
grid.draw(p1)
dev.off()


dn.overlap <- euler(c("NO2" = 10, "PM10" = 24,"PM25"=64,"NO2&PM25"=54,"PM10&PM25"=5,"NO2&PM10"=24,"NO2&PM10&PM25"=34))
plot(dn.overlap, counts = TRUE, font=1, cex=4, alpha=0.5,
     fill=c("yellow", "green","orange"),main="Overlap of negatively associated genes after NO2,PM10 & PM25 exposure",labels=list(fontsize=15,cex=2),quantities=list(cex=3))



#write.table(colnames(bryois_known_probe_norm),file="833_samples_sued_for_association_after_removing_chr_X_and_Y.txt",col.names = T,row.names = F,sep="\t",quote = F)


#write.table(no2.res.matrix,file="Association_of_NO2_with_gene_expression_adjusted_gender_gestational_age_smoking_during_pregnancy_R1_RESIDUALS.txt",col.names = T,row.names = T,sep="\t",quote = F)
#write.table(pm10.res.matrix,file="Association_of_PM10_with_gene_expression_adjusted_gender_gestational_age_smoking_during_pregnancy_R1_RESIDUALS.txt",col.names = T,row.names = T,sep="\t",quote = F)
#write.table(pm25.res.matrix,file="Association_of_PM25_with_gene_expression_adjusted_gender_gestational_age_smoking_during_pregnancy_R1_RESIDUALS.txt",col.names = T,row.names = T,sep="\t",quote = F)


write.table(no2_samples,file="Samples_used_for_NO2_association_from_837_samples.txt",col.names = F,row.names = F,sep="\t",quote = F)
write.table(pm10_samples,file="Samples_used_for_PM10_association_from_837_samples.txt",col.names = F,row.names = F,sep="\t",quote = F)
write.table(pm25_samples,file="Samples_used_for_PM25_association_from_837_samples.txt",col.names = F,row.names = F,sep="\t",quote = F)



library(org.Hs.eg.db)

pm25_probe_gene$Entrez <- mapIds(org.Hs.eg.db, pm25_probe_gene$Gene,keytype="ENSEMBL", column="ENTREZID")
pm25_probe_gene$Symbol <- mapIds(org.Hs.eg.db, pm25_probe_gene$Entrez,keytype="ENTREZID", column="SYMBOL")

pm10_probe_gene$Entrez <- mapIds(org.Hs.eg.db, pm10_probe_gene$Gene,keytype="ENSEMBL", column="ENTREZID")
pm10_probe_gene$Symbol <- mapIds(org.Hs.eg.db, pm10_probe_gene$Entrez,keytype="ENTREZID", column="SYMBOL")

no2_probe_gene$Entrez <- mapIds(org.Hs.eg.db, no2_probe_gene$Gene,keytype="ENSEMBL", column="ENTREZID")
no2_probe_gene$Symbol <- mapIds(org.Hs.eg.db, no2_probe_gene$Entrez,keytype="ENTREZID", column="SYMBOL")



pm25_probe_gene$Symbol=as.character(pm25_probe_gene$Symbol)
pm10_probe_gene$Symbol=as.character(pm10_probe_gene$Symbol)
no2_probe_gene$Symbol=as.character(no2_probe_gene$Symbol)

pm25_probe_gene$Entrez=as.character(pm25_probe_gene$Entrez)
pm10_probe_gene$Entrez=as.character(pm10_probe_gene$Entrez)
no2_probe_gene$Entrez=as.character(no2_probe_gene$Entrez)



write.table(no2_probe_gene,file="Association_of_NO2_with_gene_expression_adjusted_age9y_gender_maternaledu_areaSES_no_chr_X_Y.txt",col.names = T,row.names = F,sep="\t",quote = F)
write.table(pm25_probe_gene,file="Association_of_PM25_with_gene_expression_adjusted_age9y_gender_maternaledu_areaSES_no_chr_X_Y.txt",col.names = T,row.names = F,sep="\t",quote = F)
write.table(pm10_probe_gene,file="Association_of_PM10_with_gene_expression_adjusted_age9y_gender_maternaledu_areaSES_no_chr_X_Y.txt",col.names = T,row.names = F,sep="\t",quote = F)


#pm10_sig=subset(pm10_probe_gene,pm10_probe_gene$adj.P.Val < 0.05)

#pm10_all_genes=unique(pm10_probe_gene$Entrez)
#pm10_sig_genes=unique(pm10_sig$Entrez)


## KEGG enrichment
#library(clusterProfiler)
#kk <- enrichKEGG(gene         = pm10_sig_genes,
 #                organism     = 'hsa',
 #                pvalueCutoff = 0.05)



#pm10_go <- enrichGO(gene          = pm10_sig_genes,
      #              universe      = pm10_all_genes,
      #              OrgDb         = org.Hs.eg.db,
      #              ont           = "ALL",
       #             pAdjustMethod = "BH",
       #             pvalueCutoff  = 0.01,
       #             qvalueCutoff  = 0.05,
        #            readable      = TRUE,
        #            keyType = "ENTREZID")
##View(pm10_go@result)
#View(kk@result)


#library(ReactomePA)
#x <- enrichPathway(gene=pm10_sig_genes, pvalueCutoff = 0.05, readable=TRUE)
#head(x)
#View(x@result)






































