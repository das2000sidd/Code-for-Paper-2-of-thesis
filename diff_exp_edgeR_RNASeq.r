#RDatasetwd("~/data_das")

setwd("~/Desktop/PhD_Project_related/GENERATION R DATA RESULTS/Data/Data/THESE RESULTS FROM LAST NOVEMBER CHECK FOR LOGICAL CONCLUSIONS")
library(foreign)
library(dplyr)

pheno = read.spss('CHILD-ALLGENERALDATA_12112020.sav',to.data.frame = TRUE)
pheno2 = read.spss('CHILDSERUM9_01082017.sav',to.data.frame = TRUE)
pheno3 = read.spss('Lifecycle_GENR-9523_ap_backextrapolated_preg_post_2021-04-15.sav',to.data.frame = TRUE)
print('Line1')
raw_counts = read.table(file="GenR_counts_03-11-2021.txt",header=T,sep="\t",stringsAsFactors = F)
ses_info = read.spss(file="GENRsocialcontext_25112019.sav",to.data.frame = TRUE)
print('Line2')
raw_counts_no_X_no_Y_no_M = raw_counts[-grep("chrX|chrY|chrM",raw_counts$Chr),] ## removing all genes in Chr X and Y and M
raw_counts_no_X_no_Y_no_M=raw_counts_no_X_no_Y_no_M[,-c(2:7)]

load('GenR_RNAseq_9y_20011_176.RData')

rna.a = GenR_RNAseq_9y_20011_176

attributes(rna.a)

exp=rna.a@assayData$exprs
alignment_stats = rna.a@phenoData@data
pheno_meta = rna.a@phenoData@varMetadata
pheno_rna = rna.a@phenoData@data
length(intersect(rownames(pheno_rna),colnames(exp)))
## rownames in pheno rna are the same as in colnames of exp matrix


pheno_rna_starting_with_400=pheno_rna[grep("400([0-9]+)",rownames(pheno_rna)),]
pheno_rna_not_starting_with_400=pheno_rna[ - grep("400([0-9]+)",rownames(pheno_rna)),]

pheno_rna_starting_with_400$RNA_Id = paste("X",rownames(pheno_rna_starting_with_400),sep="")
pheno_rna_not_starting_with_400$RNA_Id = paste("X00",rownames(pheno_rna_not_starting_with_400),sep="")
print('Line3')
length(intersect(pheno_rna_starting_with_400$RNA_Id,colnames(raw_counts))) ## 96
print('Line4')
length(intersect(pheno_rna_not_starting_with_400$RNA_Id,colnames(raw_counts))) ## 80

pheno_of_interest=rbind(pheno_rna_starting_with_400,pheno_rna_not_starting_with_400)
print("Checking if bug after this")
print('Line5')

rownames(raw_counts_no_X_no_Y_no_M)=raw_counts_no_X_no_Y_no_M$Geneid
print('Line6')


print('Line7')

raw_counts_no_X_no_Y_no_M=raw_counts_no_X_no_Y_no_M[,-c(1)]

print("Checking if bug below command")

print('Line8')

raw_counts_no_X_no_Y_no_M_pheno_of_int=raw_counts_no_X_no_Y_no_M[,pheno_of_interest$RNA_Id]


rna_order=colnames(raw_counts_no_X_no_Y_no_M_pheno_of_int)
pheno_order = pheno_of_interest$RNA_Id

check_order=0
for(index in 1:length(rna_order)){
  if(rna_order[index]==pheno_order[index]){
    check_order=check_order+1
  }
}




## Need th file with batch information

## Need th file with batch information
colnames(pheno3)[grep("no2",colnames(pheno3))]
colnames(pheno2)[grep("no2",colnames(pheno2))]
colnames(pheno)[grep("no2",colnames(pheno))]


## A case insensitive search for no2
grep("NO2",colnames(pheno3),ignore.case = TRUE,value=TRUE)
grep("neigh",colnames(pheno2),ignore.case = TRUE,value=TRUE)




## NO2 - 
## First model covariates are NO2, batch - figure out from MDS, sex, age
idc_ids_use = pheno_of_interest$IDC
exp_ids_use = pheno_of_interest$RNA_Id

gene_ids=rownames(raw_counts_no_X_no_Y_no_M_pheno_of_int)

library(biomaRt)

#ensembl =useMart("ensembl")
#ensembl = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
#attributes_bm = listAttributes(ensembl)

#anno_info=getBM(attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","hgnc_symbol"),mart = ensembl)
#table(anno_info$chromosome_name)
#chrX_Y_genes=subset(anno_info,anno_info$chromosome_name %in% c("X","Y"))
#chrX_Y_genes=chrX_Y_genes[,c("chromosome_name","ensembl_gene_id")]

## Checking for gender swap

counts_use = raw_counts_no_X_no_Y_no_M_pheno_of_int

counts_use$Gene=rownames(counts_use) 

split_gene_ids = data.frame(do.call('rbind',strsplit(as.character(counts_use$Gene),'.',fixed=TRUE)))
colnames(split_gene_ids)=c("Ensembl_ID","Version")
# checking if gene ids unique,"
length(unique(split_gene_ids$Ensembl_ID))

counts_use$split_gene_ids=split_gene_ids$Ensembl_ID

check_for_repeat_ids = as.data.frame(table(counts_use$split_gene_ids))

#counts_X_Y_genes=subset(counts_use,counts_use$split_gene_ids %in% chrX_Y_genes$ensembl_gene_id)

#counts_X_Y_genes_with_chr_info=left_join(counts_X_Y_genes,chrX_Y_genes,by=c("split_gene_ids"="ensembl_gene_id"))

#rownames(counts_X_Y_genes_with_chr_info)=counts_X_Y_genes_with_chr_info$Gene


#counts_X_Y_genes_with_chr_info=counts_X_Y_genes_with_chr_info[,-c(177,178)]

#X_chr_samples=subset(counts_X_Y_genes_with_chr_info,counts_X_Y_genes_with_chr_info$chromosome_name=="X")
#Y_chr_samples=subset(counts_X_Y_genes_with_chr_info,counts_X_Y_genes_with_chr_info$chromosome_name=="Y")

#X_chr_samples=X_chr_samples[,-c(177)]
#Y_chr_samples=Y_chr_samples[,-c(177)]

#X_chr_samples=t(X_chr_samples)
#Y_chr_samples=t(Y_chr_samples)

#sum_X_chr_reads_by_sample=apply(X_chr_samples,1,sum)
#sum_Y_chr_reads_by_sample=apply(Y_chr_samples,1,sum)

#sum_X_chr_reads_by_sample=as.data.frame(sum_X_chr_reads_by_sample)
#sum_Y_chr_reads_by_sample=as.data.frame(sum_Y_chr_reads_by_sample)

#sum_X_Y_chr_reads_by_sample=cbind(sum_X_chr_reads_by_sample,sum_Y_chr_reads_by_sample)

#sum_X_Y_chr_reads_by_sample$Ratio_of_X_by_Y = sum_X_Y_chr_reads_by_sample$sum_X_chr_reads_by_sample / sum_X_Y_chr_reads_by_sample$sum_Y_chr_reads_by_sample


#plot(sum_X_Y_chr_reads_by_sample$sum_Y_chr_reads_by_sample,sum_X_Y_chr_reads_by_sample$Ratio_of_X_by_Y)

#sum_X_Y_chr_reads_by_sample$Sample_Id=rownames(sum_X_Y_chr_reads_by_sample)

#pheno_of_interest_X_Y_assesment = left_join(pheno_of_interest,sum_X_Y_chr_reads_by_sample,by=c("RNA_Id"="Sample_Id"))

#pheno_of_interest_X_Y_assesment$GENDER=as.factor(pheno_of_interest_X_Y_assesment$GENDER)

#boxplot(pheno_of_interest_X_Y_assesment$GENDER,pheno_of_interest_X_Y_assesment$Ratio_of_X_by_Y) ## Based on this, 1 is male and 2 is female


#gender_1=subset(pheno_of_interest_X_Y_assesment,pheno_of_interest_X_Y_assesment$GENDER==1)
#gender_2=subset(pheno_of_interest_X_Y_assesment,pheno_of_interest_X_Y_assesment$GENDER==2)


#summary(gender_1$Ratio_of_X_by_Y)
#summary(gender_2$Ratio_of_X_by_Y)


#pheno_of_interest_X_Y_assesment$GENDER=as.integer(pheno_of_interest_X_Y_assesment$GENDER)

#plot(pheno_of_interest_X_Y_assesment$GENDER,pheno_of_interest_X_Y_assesment$Ratio_of_X_by_Y) ## 


#length(intersect(pheno_of_interest_X_Y_assesment$RNA_Id,colnames(counts_use[,-c(177:178)])))


#rownames(pheno_of_interest_X_Y_assesment)=pheno_of_interest_X_Y_assesment$RNA_Id

pheno_order=pheno_of_interest$RNA_Id

rna_order=colnames(counts_use[,-c(177:178)])


check_order=0
for(index in 1:length(rna_order)){
  if(rna_order[index]==pheno_order[index]){
    check_order=check_order+1
  }
}

age_info=pheno2[,c("IDC","agechild9")]

## First DEseq2 and thn edgeR
cov_mod1=left_join(pheno_of_interest,age_info,by=c("IDC"))
rownames(cov_mod1)=cov_mod1$RNA_Id
#year9_no2_pm10_no2_nox = pheno3[,c("IDC","")]
colnames(pheno3)[grep("year9no2",colnames(pheno3))] ## 
colnames(pheno3)[grep("value_year9no2",colnames(pheno3))] ## 
colnames(pheno)[grep("no2",colnames(pheno))] ## nothing
colnames(pheno2)[grep("no2",colnames(pheno))] ## nothing

sum(is.na(cov_mod1$GENDER))
sum(is.na(cov_mod1$agechild9))
# Confirmed no NA values for covariates in first model




cov_mod1$GENDER=as.factor(cov_mod1$GENDER)

## first analysing with DESeq2

library(edgeR)

pheno3_no2=pheno3[,c(grep("value_year9no2_be",colnames(pheno3)))]
pheno3_no2=as.data.frame(pheno3_no2)
pheno3_no2$IDC=pheno3$IDC
#pheno3_no2=pheno3_no2[,c(11,1:10)]
#year_9_avg_pm=apply(pheno3_no2[,c(2:11)],1,mean)
#pheno3_no2$year_9_avg_pm=year_9_avg_pm
#pheno3_no2=pheno3_no2[,c("IDC","year_9_avg_pm")]

cov_mod2=left_join(cov_mod1,pheno3_no2,by=c("IDC"))

cov_mod2$year_9_avg_pm_centered = scale(cov_mod2$pheno3_no2,center = TRUE)
cov_mod2$age_child_9_centered = scale(cov_mod2$agechild9,center = TRUE)
colnames(cov_mod2)
colnames(cov_mod2)[25] = "year_9_avg_pm_centered"
rownames(cov_mod2)=cov_mod2$RNA_Id

table(colnames(counts_use[,c(1:176)])==rownames(cov_mod2))

cov_mod2$flowcellID=as.factor(cov_mod2$flowcellID)
library(DESeq2)
## testing to see if gender drives any differences
dds = DESeqDataSetFromMatrix(countData = counts_use[,c(1:176)],colData = cov_mod2,design = ~ GENDER + year_9_avg_pm_centered + age_child_9_centered + flowcellID)


## removed two samples
remove_samp = c("X4006673642","X0011103401")

`%ni%`=Negate(`%in%`)
raw_counts_no_outlier = counts_use[,colnames(counts_use) %ni% remove_samp]
cov_mod2_no_outlier = subset(cov_mod2, cov_mod2$RNA_Id %ni% remove_samp)


rna_order=colnames(raw_counts_no_outlier[,c(1:174)])
pheno_order=cov_mod2_no_outlier$RNA_Id


check_order=0
for(index in 1:length(rna_order)){
  if(rna_order[index]==pheno_order[index]){
    check_order=check_order+1
  }
}


dds = DESeqDataSetFromMatrix(countData = raw_counts_no_outlier[,c(1:174)],colData = cov_mod2_no_outlier,design = ~ GENDER + year_9_avg_pm_centered + age_child_9_centered + flowcellID)


MDs.with.edgeR=plotMDS(dds)


MDS.table=MDs.with.edgeR$eigen.vectors

PC1_PC2=MDS.table[,c(1,2)]

## PC plot looks okay after removing those two outliers


#rownames(PC1_PC2)=colnames(counts_use_no_outlier)[1:174]
#PC1_PC2=as.data.frame(PC1_PC2)

#plot(PC1_PC2$V1,PC1_PC2$V2,main="Using plotmDS funciton of egdeR")




#PC1_PC2$Exp_Id=rownames(PC1_PC2)

#Adjust for PC1 and PC2
#colnames(PC1_PC2)[1:2]=c("PC1","PC2")


count_as_matrix = as.matrix(raw_counts_no_outlier[,c(1:174)])
rownames(count_as_matrix)=rownames(raw_counts_no_outlier[,c(1:174)])
colnames(count_as_matrix)=colnames(raw_counts_no_outlier[,c(1:174)])

## differential expression analysis model 1
dds = DGEList(counts = raw_counts_no_outlier[,c(1:174)],remove.zeros = TRUE)
keep = rowSums(cpm(dds) > 1)  >= 2
dds <- dds[keep,keep.lib.sizes=FALSE]
dds <- calcNormFactors(dds)

table(colnames(dds)==rownames(cov_mod2_no_outlier))

design_mat = model.matrix( ~ GENDER + year_9_avg_pm_centered + age_child_9_centered + flowcellID , data = cov_mod2_no_outlier)

dds = estimateGLMCommonDisp(dds, design_mat)

dds = estimateGLMTrendedDisp(dds, design_mat)

dds = estimateGLMTagwiseDisp(dds, design_mat)



fit_edgeR = glmFit (dds, design_mat)
lrt.test = glmLRT(fit_edgeR, coef="year_9_avg_pm_centered") ## previously i had 2 here
lrt.test.table = lrt.test$table

lrt.test.table$p.adjusted=p.adjust(lrt.test.table$PValue,method = "BH")


write.table(lrt.test.table,file="Asso_gene_exp_NO2_adjusted_for_gender_age_flowcellid_as_batch_variable.txt",col.names=T,row.names=T,sep="\t",quote=FALSE)

save.image("Latest_session_for_analysis_EdgeR.RData")

load('Latest_session_for_analysis_EdgeR.RData')





ses9=ses_info[,c("IDC","areases_quint_9","areases_tert_9")]
mother_edu=pheno[,c("IDC","EDUCM")]

cov_mod2_no_outlier_ses=left_join(cov_mod2_no_outlier,ses9,by=c("IDC"))
cov_mod2_no_outlier_ses_educm=left_join(cov_mod2_no_outlier_ses,mother_edu,by=c("IDC"))

cov_mod2_no_outlier_ses_educm=cov_mod2_no_outlier_ses_educm[complete.cases(cov_mod2_no_outlier_ses_educm),]

edu_1=subset(cov_mod2_no_outlier_ses_educm,cov_mod2_no_outlier_ses_educm$EDUCM %in% c("primary","secondary, phase 1","secondary, phase 2"))
edu_2=subset(cov_mod2_no_outlier_ses_educm,cov_mod2_no_outlier_ses_educm$EDUCM %in% c("higher, phase 1","higher, phase 2"))

edu_1$EDUCM="Low"
edu_2$EDUCM="High"

cov_mod2_no_outlier_ses_educm_new=rbind(edu_1,edu_2)
rownames(cov_mod2_no_outlier_ses_educm_new)=cov_mod2_no_outlier_ses_educm_new$RNA_Id
cov_mod2_no_outlier_ses_educm_new$areases_tert_9=as.factor(cov_mod2_no_outlier_ses_educm_new$areases_tert_9)

samples_use=rownames(cov_mod2_no_outlier_ses_educm_new)
raw_counts_no_outlier_2=raw_counts_no_outlier[,samples_use]

table(colnames(raw_counts_no_outlier_2)==rownames(cov_mod2_no_outlier_ses_educm_new))

dds2 = DGEList(counts = raw_counts_no_outlier_2,remove.zeros = TRUE)
keep = rowSums(cpm(dds2) > 1)  >= 2
dds2 <- dds2[keep,keep.lib.sizes=FALSE]
dds2 <- calcNormFactors(dds2)

table(colnames(dds2)==rownames(cov_mod2_no_outlier_ses_educm_new))

design_mat = model.matrix( ~ GENDER + year_9_avg_pm_centered + age_child_9_centered + flowcellID + EDUCM + areases_tert_9 , data = cov_mod2_no_outlier_ses_educm_new)


dds2 = estimateGLMCommonDisp(dds2, design_mat)

dds2 = estimateGLMTrendedDisp(dds2, design_mat)

dds2 = estimateGLMTagwiseDisp(dds2, design_mat)


fit_edgeR = glmFit (dds2, design_mat)
lrt.test = glmLRT(fit_edgeR, coef="year_9_avg_pm_centered") ## previously i had 2 here
lrt.test.table = lrt.test$table

lrt.test.table$p.adjusted=p.adjust(lrt.test.table$PValue,method = "BH")


write.table(lrt.test.table,file="Asso_gene_exp_NO2_adjusted_for_gender_age_flowcellid_as_batch_variable_maternal_edu_as_SES_and_areaSES.txt",col.names=T,row.names=T,sep="\t",quote=FALSE)



celltype_IDC_connector = read.spss('Selection_GENR_450kmeth_release3_9y_20220816.sav',to.data.frame = TRUE) ## connects IDC with meth ID
celltype = read.table(file='GENR_450kmeth_release3_9y_Houseman_default.txt',header=TRUE,sep="\t",row.names=1) ## 
rownames(celltype_IDC_connector) = celltype_IDC_connector$Sample_ID
celltype_with_celltype_idc_coonect = cbind(celltype,celltype_IDC_connector[rownames(celltype),])
celltype_with_celltype_idc_coonect=celltype_with_celltype_idc_coonect[complete.cases(celltype_with_celltype_idc_coonect),]
celltype_with_celltype_idc_coonect=celltype_with_celltype_idc_coonect[,-c(7,8,10,11)]

## after this adjusted for NO2, batch, sex, age, SES, celltype
cov_mod2_no_outlier_with_cell_type=left_join(cov_mod2_no_outlier_ses_educm_new,celltype_with_celltype_idc_coonect,by=c("IDC"="IDC"))
cov_mod2_no_outlier_with_cell_type=cov_mod2_no_outlier_with_cell_type[complete.cases(cov_mod2_no_outlier_with_cell_type),]
#cov_mod2_cell_type_not_missing_ses=left_join(cov_mod2_cell_type_not_missing,ses9,by=c("IDC")) # 128
#cov_mod2_cell_type_not_missing_ses$areases_quint_9 = as.factor(cov_mod2_cell_type_not_missing_ses$areases_quint_9)
#cov_mod2_cell_type_not_missing_ses$areases_tert_9 = as.factor(cov_mod2_cell_type_not_missing_ses$areases_tert_9)

## centering and scaling celltype vars
cov_mod2_no_outlier_with_cell_type$CD8T_centered = scale(cov_mod2_no_outlier_with_cell_type$CD8T,center = TRUE)
cov_mod2_no_outlier_with_cell_type$CD4T_centered = scale(cov_mod2_no_outlier_with_cell_type$CD4T,center = TRUE)
cov_mod2_no_outlier_with_cell_type$NK_centered = scale(cov_mod2_no_outlier_with_cell_type$NK,center = TRUE)
cov_mod2_no_outlier_with_cell_type$Bcell_centered = scale(cov_mod2_no_outlier_with_cell_type$Bcell,center = TRUE)
cov_mod2_no_outlier_with_cell_type$Mono_centered = scale(cov_mod2_no_outlier_with_cell_type$Mono,center = TRUE)
cov_mod2_no_outlier_with_cell_type$Gran_centered = scale(cov_mod2_no_outlier_with_cell_type$Gran,center = TRUE)

rownames(cov_mod2_no_outlier_with_cell_type)=cov_mod2_no_outlier_with_cell_type$RNA_Id

View(cov_mod2_no_outlier_with_cell_type)

#cov_mod2_cell_type_not_missing_ses=cov_mod2_cell_type_not_missing_ses[complete.cases(cov_mod2_cell_type_not_missing_ses),] ## since SES misisng
#dim(cov_mod2_cell_type_not_missing_ses)

cov_mod2_no_outlier_with_cell_type$areases_quint_9 = as.factor(cov_mod2_no_outlier_with_cell_type$areases_quint_9)

design_mat_2 = model.matrix( ~ GENDER + year_9_avg_pm_centered + age_child_9_centered + flowcellID + areases_quint_9 + CD8T_centered + CD4T_centered + NK_centered + Bcell_centered + Mono_centered + Gran_centered + EDUCM , data = cov_mod2_no_outlier_with_cell_type)
rna_seq_samp_use = cov_mod2_no_outlier_with_cell_type$RNA_Id
raw_counts_no_outlier_use=raw_counts_no_outlier[,rna_seq_samp_use]
dim(raw_counts_no_outlier_use)
dim(design_mat_2)
rownames(design_mat_2) = cov_mod2_no_outlier_with_cell_type$RNA_Id

table(rownames(design_mat_2) == colnames(raw_counts_no_outlier_use)) ## All true

## differential expression analysis model 2


dds = DGEList(counts = raw_counts_no_outlier_use,remove.zeros = TRUE)
keep = rowSums(cpm(dds) > 1)  >= 2
dds <- dds[keep,keep.lib.sizes=FALSE]
dds <- calcNormFactors(dds)

dds = estimateGLMCommonDisp(dds, design_mat_2)

dds = estimateGLMTrendedDisp(dds, design_mat_2)

dds = estimateGLMTagwiseDisp(dds, design_mat_2)


fit_edgeR = glmFit (dds, design_mat_2)
lrt.test = glmLRT(fit_edgeR, coef="year_9_avg_pm_centered")
lrt.test.table = lrt.test$table


lrt.test.table$p.adjusted=p.adjust(lrt.test.table$PValue,method = "BH")


write.table(lrt.test.table,file="Asso_gene_exp_NO2_adjusted_for_gender_age_flowcellid_as_batch_variable_maternal_edu_as_SES_and_areaSES_celltype.txt",col.names=T,row.names=T,sep="\t",quote=FALSE)


## differential expression analysis model 2






