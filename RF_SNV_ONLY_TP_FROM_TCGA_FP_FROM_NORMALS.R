#NUevo script par SNV con FP solo los de los empty bams no de TCGA sera en principio mejor

# RANDOM FOREST FOR SNV
#Comments: adaptation of Romina's script and my script for InDels for TCGA data. 
#Esto hacerlo como si tuviera el archivo de TCGA de indels pero para snvs como el de montse. 
library(viridis)
library(data.table)
library(randomForest)
library(UpSetR)
library(caret)
library(randomForestExplainer)
library(reprtree)
library(catspec)
library(devtools)
library(tidyr)
library(dplyr)
library(rpart)
library(rpart.plot)
library(rattle)
library(ggplot2)
library(e1071)
library(openssl)
library(sqldf)
library(plyr)
library(rattle)
library(export)
#devtools::install_github('araastat/reprtree')
#install.packages("openssl")

#snvs_result <- get(load("snvs_empty_bams.RData"))





SNVS_TCGA <- "/media/computer/HDD/MOUNT/SNVS_UCEC_CESC.tab"
snvs_tcga <- fread(SNVS_TCGA, header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))
colnames(snvs_tcga)
#DNA1 <- snvs_result
# fileGoldenT4 <- "/media/computer/HDD/MOUNT/Golden_all"
# goldenT4 <- fread(fileGoldenT4,header=TRUE)


#/gpfs/scratch/bsc96/bsc96452/GENOMICS/peris/benchmarking/data_donors/WES_NV_N/VC
#/gpfs/scratch/bsc96/bsc96452/GENOMICS/peris/benchmarking/data_donors/WES_NV_N/VC/WES_NV_N_1_vs_WES_NV_N_2/final/combine/indels/indels.PASS.ANNOT.GRCh37.p13.RefSeq.tab



moreFP_1 <- "/media/computer/HDD/MOUNT/moreFP/snvs.WES_NV_N_1_vs_WES_NV_N_2.PASS.ANNOT.GRCh37.p13.RefSeq.tab"
moreFP_2 <- "/media/computer/HDD/MOUNT/moreFP/snvs.WES_NV_N_1_vs_WES_NV_N_3.PASS.ANNOT.GRCh37.p13.RefSeq.tab"
moreFP_3 <- "/media/computer/HDD/MOUNT/moreFP/snvs.WES_NV_N_2_vs_WES_NV_N_1.PASS.ANNOT.GRCh37.p13.RefSeq.tab"
moreFP_4 <- "/media/computer/HDD/MOUNT/moreFP/snvs.WES_NV_N_2_vs_WES_NV_N_3.PASS.ANNOT.GRCh37.p13.RefSeq.tab"
moreFP_5 <- "/media/computer/HDD/MOUNT/moreFP/snvs.WES_NV_N_3_vs_WES_NV_N_1.PASS.ANNOT.GRCh37.p13.RefSeq.tab"
moreFP_6 <- "/media/computer/HDD/MOUNT/moreFP/snvs.WES_NV_N_3_vs_WES_NV_N_2.PASS.ANNOT.GRCh37.p13.RefSeq.tab"

MoreFP_1 <- fread(moreFP_1, header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))
MoreFP_2 <- fread(moreFP_2, header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))
MoreFP_3 <- fread(moreFP_3, header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))
MoreFP_4 <- fread(moreFP_4, header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))
MoreFP_5 <- fread(moreFP_5, header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))
MoreFP_6 <- fread(moreFP_6, header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))

colnames(MoreFP_1) <- gsub(paste0("WES_NV_N_1", "."), "N_", colnames(MoreFP_1))
colnames(MoreFP_1) <- gsub(paste0("WES_NV_N_2", "."), "T_", colnames(MoreFP_1))

colnames(MoreFP_2) <- gsub(paste0("WES_NV_N_1", "."), "N_", colnames(MoreFP_2))
colnames(MoreFP_2) <- gsub(paste0("WES_NV_N_3", "."), "T_", colnames(MoreFP_2))

colnames(MoreFP_3) <- gsub(paste0("WES_NV_N_2", "."), "N_", colnames(MoreFP_3))
colnames(MoreFP_3) <- gsub(paste0("WES_NV_N_1", "."), "T_", colnames(MoreFP_3))

colnames(MoreFP_4) <- gsub(paste0("WES_NV_N_2", "."), "N_", colnames(MoreFP_4))
colnames(MoreFP_4) <- gsub(paste0("WES_NV_N_3", "."), "T_", colnames(MoreFP_4))

colnames(MoreFP_5) <- gsub(paste0("WES_NV_N_3", "."), "N_", colnames(MoreFP_5))
colnames(MoreFP_5) <- gsub(paste0("WES_NV_N_1", "."), "T_", colnames(MoreFP_5))

colnames(MoreFP_6) <- gsub(paste0("WES_NV_N_3", "."), "N_", colnames(MoreFP_6))
colnames(MoreFP_6) <- gsub(paste0("WES_NV_N_2", "."), "T_", colnames(MoreFP_6))

MoreFP_1$T_SAMPLE <- "N_1_vs_N_2"
MoreFP_2$T_SAMPLE <- "N_1_vs_N_3"
MoreFP_3$T_SAMPLE <- "N_2_vs_N_1"
MoreFP_4$T_SAMPLE <- "N_2_vs_N_3"
MoreFP_5$T_SAMPLE <- "N_3_vs_N_1"
MoreFP_6$T_SAMPLE <- "N_3_vs_N_2"


ncol(MoreFP_1)
ncol(snvs_tcga)
colnames(MoreFP_6)
colnames(snvs_tcga)

table(MoreFP_6$ALT)
#Aquí se observa que la mayoria de los nt mutados son bien C o G y mutan a A o T mayoritariamente
#Esto porque sera? pq la estructura quimica? por el entorno? la CG tienen mas fuerza que AT cómo puede ser? 



all<- rbind(snvs_tcga, MoreFP_1, MoreFP_2, MoreFP_3,MoreFP_4,MoreFP_5,MoreFP_6)
table(snvs_tcga$T_GT)
table(MoreFP_1$T_GT)
table(snvs_tcga$N_GT)
table(MoreFP_1$N_GT)


#Golden:



fileGoldenT4 <- "/media/computer/HDD/MOUNT/Golden_all"
goldenT4 <- fread(fileGoldenT4,header=TRUE)
table(goldenT4$Variant_Type)
GOLDEN_TCGA_SNVS <- goldenT4[goldenT4$Variant_Type=="SNP",]

GOLDEN_TCGA_SNVS2 <- goldenT4[goldenT4$Variant_Type=="SNP" & goldenT4$mutval_status_targeted =="validated_powered" ,]
GOLDEN_TCGA_SNVS3 <- goldenT4[goldenT4$Variant_Type=="SNP" & goldenT4$mutval_status_targeted =="unvalidated_powered" ,]

table(GOLDEN_TCGA_SNVS$Reference_Allele)
table(GOLDEN_TCGA_SNVS$mutval_status_targeted)
# > table(GOLDEN_TCGA_SNVS$Reference_Allele)
# A     C     G     T 
# 3534 18065 18303  3648 

# > table(GOLDEN_TCGA_SNVS$Tumor_Seq_Allele1)
# A     C     G     T 
# 3534 18065 18303  3648 

# > table(GOLDEN_TCGA_SNVS$Tumor_Seq_Allele2)
# A     C     G     T 
# 15524  6275  6188 15563 

#Aquí también se observa que la mayoria de los nt mutados son bien C o G y mutan a A o T mayoritariamente. Los FP eran menos (769)
#Esto porque sera? pq la estructura quimica? por el entorno? la CG tienen mas fuerza que AT cómo puede ser? 


GOLDEN_TCGA_SNVS$mut <- paste0(GOLDEN_TCGA_SNVS$Chromosome,"_",GOLDEN_TCGA_SNVS$Start_Position,"_",GOLDEN_TCGA_SNVS$Reference_Allele,"_",GOLDEN_TCGA_SNVS$Tumor_Seq_Allele2,"_",GOLDEN_TCGA_SNVS$Tumor_Sample_Barcode)




DNA0 <- as.data.frame(all)

#Here, I think that you havae to merge colums with _ CHROM and POS in order to create a new column to compare later
#Due to the fact that several sample are being used, it is also required to merge using T_SAMPLE
DNA0$mut <- paste0(DNA0$CHROM,"_",DNA0$POSITION,"_",DNA0$REF,"_",DNA0$ALT,"_",DNA0$T_SAMPLE)

DNA1 <- DNA0

#Sería mejor hacerlo con el file, ya directamente....
head(DNA1)
a <- colnames(DNA1)
list(colnames(DNA1))

mode(a)

#Aquí hace falta seleccionar las columnas adecuadas de snvs no de indels:
#Que hay de la normalizacion de DP? no see me tocara probar las diferentes posibilidades primero sin dnormaliar ya que no confio en eso

selectcols_snv <- c("CHROM", "POSITION", "REF", "ALT", "T_NPROGS", "T_PROGS", "N_DP_avg", "N_DP_std", "N_DP_REF_avg", 
                    "N_DP_ALT_avg", "N_VAF_avg", "N_DP_REF_cav", "N_DP_ALT_cav", "N_DP_REF_g4m2", "N_DP_ALT_g4m2", 
                    "N_MBQ_g4m2", "N_MMQ_g4m2", "N_DP_REF_str", "N_DP_ALT_str", "N_DP_REF_mus", "N_DP_ALT_mus", 
                    "N_BQ_REF_mus", "N_BQ_ALT_mus", "N_DP_REF_lan", "N_DP_ALT_lan", "T_NPROGS", "T_PROGS", "T_DP_avg", 
                    "T_DP_std", "T_DP_REF_avg", "T_DP_ALT_avg", "T_VAF_avg", "T_DP_REF_cav", "T_DP_ALT_cav", "T_VAF_cav", 
                    "T_ASMD_cav", "T_CLPM_cav", "T_DP_REF_g4m2", "T_DP_ALT_g4m2", "T_VAF_g4m2", "T_MBQ_g4m2", "T_MMQ_g4m2", 
                    "T_DP_REF_str", "T_DP_ALT_str", "T_VAF_str", "T_MQ_str", "T_MQ0_str", "T_SomaticEVS_str", "T_DP_REF_mus", 
                    "T_DP_ALT_mus", "T_VAF_mus", "T_BQ_REF_mus", "T_BQ_ALT_mus", "T_DP_REF_lan", "T_DP_ALT_lan", "T_VAF_lan", 
                    "T_FETS_lan", "T_SB_lan", "T_SAMPLE", "mut")


#Se eliminan las que tienen el mismo valor en todas las rows muchas de las N_ y las de GT ya que o eran 0/0 en N o 0/1 en T
#Las de los T_FILTER tb fuera ya que estan totalmente correlacionadas con el resto en el setndio de que solo tendran PASS cuando el resto tengan info

selectCols <- c("mut", "CHROM","POSITION","REF","ALT","N_GT", "T_NPROGS","T_PROGS",
                "N_DP_avg","N_DP_std","N_DP_REF_avg",           "N_DP_ALT_avg",           "N_VAF_avg",  "N_DP_REF_g4m2",          "N_DP_ALT_g4m2",          "N_VAF_g4m2",
                "N_MBQ_g4m2", "N_MMQ_g4m2", "N_IN_PON_g4m2",          "N_FILTER_g4m2",          "N_DP_REF_str",           "N_DP_ALT_str",          
                "N_VAF_str",  "N_MQ_str",   "N_MQ0_str",  "N_SomaticEVS_str",       "N_RU_str",   "N_RC_str",  
                "N_IC_str",   "N_IHP_str",  "N_FILTER_str",           "N_DP_REF_lan",           "N_DP_ALT_lan",           "N_VAF_lan", 
                "N_FETS_lan", "N_SB_lan",   "N_FILTER_lan",          "N_DP_REF_pin",           "N_DP_ALT_pin",           "N_VAF_pin",  "N_QUAL_pin", "N_S1_pin",              
                "N_S2_pin",         "N_FILTER_pin",     "N_DP_REF_pla",     "N_DP_ALT_pla",     "N_VAF_pla",        "N_GQ_pla",        
                "N_GOF_pla",        "N_MMLQ_pla",       "N_HP_pla",         "N_MQ_pla",         "N_QD_pla",         "N_FILTER_pla",    
                "N_DP_REF_sva",     "N_DP_ALT_sva",     "N_VAF_sva",        "N_LO_sva",         "N_LR_sva",         "N_MAPQ_sva",      
                "N_LOD_sva",        "N_REPSEQ_sva",     "N_FILTER_sva",     "T_GT",             "T_NPROGS",         "T_PROGS",         
                "T_DP_avg",         "T_DP_std",         "T_DP_REF_avg",     "T_DP_ALT_avg",     "T_VAF_avg",        "T_DP_REF_g4m2",   
                "T_DP_ALT_g4m2",    "T_VAF_g4m2",       "T_MBQ_g4m2",       "T_MMQ_g4m2",       "T_IN_PON_g4m2",    "T_FILTER_g4m2",   
                "T_DP_REF_str",     "T_DP_ALT_str",     "T_VAF_str",        "T_MQ_str",         "T_MQ0_str",        "T_SomaticEVS_str",
                "T_RU_str",         "T_RC_str",         "T_IC_str",         "T_IHP_str",        "T_FILTER_str",     "T_DP_REF_lan",    
                "T_DP_ALT_lan",     "T_VAF_lan",        "T_FETS_lan",       "T_SB_lan",         "T_FILTER_lan",     "T_DP_REF_pin",     "T_DP_ALT_pin",     "T_VAF_pin",       
                "T_QUAL_pin",       "T_S1_pin",         "T_S2_pin",         "T_FILTER_pin",     "T_DP_REF_pla",     "T_DP_ALT_pla",    
                "T_VAF_pla",        "T_GQ_pla",         "T_GOF_pla",        "T_MMLQ_pla",       "T_HP_pla",         "T_MQ_pla",        
                "T_QD_pla",         "T_FILTER_pla",     "T_DP_REF_sva",     "T_DP_ALT_sva",     "T_VAF_sva",        "T_LO_sva",        
                "T_LR_sva",         "T_MAPQ_sva",       "T_LOD_sva",        "T_REPSEQ_sva",     "T_FILTER_sva",     "T_SAMPLE")



muts34 <- DNA1[,names(DNA1) %in% selectcols_snv]
table(muts34$T_PROGS)
summary(muts34)

muts34$T_PROGS <- as.factor(muts34$T_PROGS)

summary(muts34)
#Con esto ya estan los datos preparados a falta de: DP normalization, eliminar correlacionadas, etc
#Para el training habra que dejar solo las columnas del variant calling nada de posiciones etc...



muts2 <- muts34
ncol(muts2) #58
summary(muts2)

validated <- GOLDEN_TCGA_SNVS[GOLDEN_TCGA_SNVS$mutval_status_targeted == "validated_powered",] # 42781
unvalidated <- GOLDEN_TCGA_SNVS[GOLDEN_TCGA_SNVS$mutval_status_targeted == "unvalidated_powered",] # 769

table(muts2$T_PROGS)


progs_1_more <- c("1","2","3","4","5")
progs_2_more <- c("2","3","4","5")
progs_3_more <- c("3","4","5")
progs_4_more <- c("4","5")
progs_5_more <- c("5")
progs_gatk4m2_more <- c("gatk4m2")
progs_lancet_more  <- c("lancet")
progs_caveman_more <- c("caveman")
progs_muse_more <- c("muse")
progs_strelka_more <- c("strelka")



muts_1 <- muts2[muts2$T_NPROGS %in% progs_1_more,]
muts_2 <- muts2[muts2$T_NPROGS %in% progs_2_more,]
muts_3 <- muts2[muts2$T_NPROGS %in% progs_3_more,]
muts_4 <- muts2[muts2$T_NPROGS %in% progs_4_more,]
muts_5 <- muts2[muts2$T_NPROGS %in% progs_5_more,]

muts_gatk4m2 <- muts2[muts2$T_PROGS %like% progs_gatk4m2_more,]
muts_lancet  <- muts2[muts2$T_PROGS %like% progs_lancet_more,]
muts_caveman <- muts2[muts2$T_PROGS %like% progs_caveman_more,]
muts_muse    <- muts2[muts2$T_PROGS %like% progs_muse_more,]
muts_strelka <- muts2[muts2$T_PROGS %like% progs_strelka_more,]





TP_1a <- validated[validated$mut %in% muts_1$mut,] # 42108 #Orden no importa aquí
TP_1 <- muts_1[muts_1$mut %in% validated$mut,] #esto da lo mismo 42108
FN_1 <- anti_join(validated, TP_1a) # 673 #Orden si importa aquí
# FP_TCGA_1 <- muts_1[muts_1$mut %in% unvalidated$mut,] # 411 
FP_empty_bams_1 <- muts_1[muts_1$T_SAMPLE %like% "_vs_",]
FP_1 <- FP_empty_bams_1
nrow(FP_1) #7798
TN_1 <- unvalidated[!(unvalidated$mut %in% muts_1$mut),] #358
Recall_1 <- nrow(TP_1)/nrow(validated) *100
#If TN from the normal-normal datasets are not included, specificity values will heavily depend on FP rates
Specificity_1 <-nrow(TN_1) / (nrow(TN_1) + nrow(FP_1)) *100
Precision_1 <-nrow(TP_1)/(nrow(TP_1)+nrow(FP_1)) *100
#The number of FP is 6-fold lower than the number of TP!!
# > nrow(FP_1)
# [1] 7798
# > nrow(TP_1)
# [1] 42108
Accuracy_1 <-100*(nrow(TN_1)+nrow(TP_1))/(nrow(TN_1)+nrow(FN_1)+nrow(TP_1)+nrow(FP_1))
F1_Score_1 <- 2 * Precision_1 * Recall_1 / (Precision_1 +Recall_1)
MCC_1 <- (as.numeric(nrow(TP_1))*as.numeric(nrow(TN_1)) -as.numeric(nrow(FP_1))*as.numeric(nrow(FN_1))) / sqrt((as.numeric(nrow(TP_1))+as.numeric(nrow(FP_1)))*(as.numeric(nrow(TP_1))+as.numeric(nrow(FN_1)))*(as.numeric(nrow(TN_1))+as.numeric(nrow(FP_1)))*(as.numeric(nrow(TN_1))+as.numeric(nrow(FN_1))))


TP_2a <- validated[validated$mut %in% muts_2$mut,] 
TP_2 <- muts_2[muts_2$mut %in% validated$mut,]
FN_2 <- anti_join(validated, TP_2a) 
# FP_TCGA_2 <- muts_2[muts_2$mut %in% unvalidated$mut,] 
FP_empty_bams_2 <- muts_2[muts_2$T_SAMPLE %like% "_vs_",]
FP_2 <- FP_empty_bams_2
TN_2 <- unvalidated[!(unvalidated$mut %in% muts_2$mut),] 
Recall_2 <- nrow(TP_2)/nrow(validated) *100
Specificity_2 <-nrow(TN_2) / (nrow(TN_2) + nrow(FP_2)) *100
Precision_2 <-nrow(TP_2)/(nrow(TP_2)+nrow(FP_2)) *100
Accuracy_2 <-100*(nrow(TN_2)+nrow(TP_2))/(nrow(TN_2)+nrow(FN_2)+nrow(TP_2)+nrow(FP_2))
F1_Score_2 <- 2 * Precision_2 * Recall_2 / (Precision_2 +Recall_2)
MCC_2 <- (as.numeric(nrow(TP_2))*as.numeric(nrow(TN_2)) -as.numeric(nrow(FP_2))*as.numeric(nrow(FN_2))) / sqrt((as.numeric(nrow(TP_2))+as.numeric(nrow(FP_2)))*(as.numeric(nrow(TP_2))+as.numeric(nrow(FN_2)))*(as.numeric(nrow(TN_2))+as.numeric(nrow(FP_2)))*(as.numeric(nrow(TN_2))+as.numeric(nrow(FN_2))))

TP_3a <- validated[validated$mut %in% muts_3$mut,] 
TP_3 <- muts_3[muts_3$mut %in% validated$mut,]
FN_3 <- anti_join(validated, TP_3a) 
# FP_TCGA_3 <- muts_3[muts_3$mut %in% unvalidated$mut,] 
FP_empty_bams_3 <- muts_3[muts_3$T_SAMPLE %like% "_vs_",]
FP_3 <- FP_empty_bams_3
TN_3 <- unvalidated[!(unvalidated$mut %in% muts_3$mut),] 
Recall_3 <- nrow(TP_3)/nrow(validated) *100
Specificity_3 <-nrow(TN_3) / (nrow(TN_3) + nrow(FP_3)) *100
Precision_3 <-nrow(TP_3)/(nrow(TP_3)+nrow(FP_3)) *100
Accuracy_3 <-100*(nrow(TN_3)+nrow(TP_3))/(nrow(TN_3)+nrow(FN_3)+nrow(TP_3)+nrow(FP_3))
F1_Score_3 <- 2 * Precision_3 * Recall_3 / (Precision_3 +Recall_3)
MCC_3 <- (as.numeric(nrow(TP_3))*as.numeric(nrow(TN_3)) -as.numeric(nrow(FP_3))*as.numeric(nrow(FN_3))) / sqrt((as.numeric(nrow(TP_3))+as.numeric(nrow(FP_3)))*(as.numeric(nrow(TP_3))+as.numeric(nrow(FN_3)))*(as.numeric(nrow(TN_3))+as.numeric(nrow(FP_3)))*(as.numeric(nrow(TN_3))+as.numeric(nrow(FN_3))))

TP_4a <- validated[validated$mut %in% muts_4$mut,] 
TP_4 <- muts_4[muts_4$mut %in% validated$mut,]
FN_4 <- anti_join(validated, TP_4a) 
# FP_TCGA_4 <- muts_4[muts_4$mut %in% unvalidated$mut,] 
FP_empty_bams_4 <- muts_4[muts_4$T_SAMPLE %like% "_vs_",]
FP_4 <- FP_empty_bams_4
TN_4 <- unvalidated[!(unvalidated$mut %in% muts_4$mut),] 
Recall_4 <- nrow(TP_4)/nrow(validated) *100
Specificity_4 <-nrow(TN_4) / (nrow(TN_4) + nrow(FP_4)) *100
Precision_4 <-nrow(TP_4)/(nrow(TP_4)+nrow(FP_4)) *100
Accuracy_4 <-100*(nrow(TN_4)+nrow(TP_4))/(nrow(TN_4)+nrow(FN_4)+nrow(TP_4)+nrow(FP_4))
F1_Score_4 <- 2 * Precision_4 * Recall_4 / (Precision_4 +Recall_4)
MCC_4 <- (as.numeric(nrow(TP_4))*as.numeric(nrow(TN_4)) -as.numeric(nrow(FP_4))*as.numeric(nrow(FN_4))) / sqrt((as.numeric(nrow(TP_4))+as.numeric(nrow(FP_4)))*(as.numeric(nrow(TP_4))+as.numeric(nrow(FN_4)))*(as.numeric(nrow(TN_4))+as.numeric(nrow(FP_4)))*(as.numeric(nrow(TN_4))+as.numeric(nrow(FN_4))))

TP_5a <- validated[validated$mut %in% muts_5$mut,] 
TP_5 <- muts_5[muts_5$mut %in% validated$mut,]
FN_5 <- anti_join(validated, TP_5a) 
# FP_TCGA_5 <- muts_5[muts_5$mut %in% unvalidated$mut,] 
FP_empty_bams_5 <- muts_5[muts_5$T_SAMPLE %like% "_vs_",]
FP_5 <- FP_empty_bams_5
TN_5 <- unvalidated[!(unvalidated$mut %in% muts_5$mut),] 
Recall_5 <- nrow(TP_5)/nrow(validated) *100
Specificity_5 <-nrow(TN_5) / (nrow(TN_5) + nrow(FP_5)) *100
Precision_5 <-nrow(TP_5)/(nrow(TP_5)+nrow(FP_5)) *100
Accuracy_5 <-100*(nrow(TN_5)+nrow(TP_5))/(nrow(TN_5)+nrow(FN_5)+nrow(TP_5)+nrow(FP_5))
F1_Score_5 <- 2 * Precision_5 * Recall_5 / (Precision_5 +Recall_5)
MCC_5 <- (as.numeric(nrow(TP_5))*as.numeric(nrow(TN_5)) -as.numeric(nrow(FP_5))*as.numeric(nrow(FN_5))) / sqrt((as.numeric(nrow(TP_5))+as.numeric(nrow(FP_5)))*(as.numeric(nrow(TP_5))+as.numeric(nrow(FN_5)))*(as.numeric(nrow(TN_5))+as.numeric(nrow(FP_5)))*(as.numeric(nrow(TN_5))+as.numeric(nrow(FN_5))))


#Falta la de los programas


TP_gatk4m2a <- validated[validated$mut %in% muts_gatk4m2$mut,] 
TP_gatk4m2 <- muts_gatk4m2[muts_gatk4m2$mut %in% validated$mut,] 
FN_gatk4m2 <- anti_join(validated, TP_gatk4m2a) 
# FP_TCGA_gatk4m2 <- muts_gatk4m2[muts_gatk4m2$mut %in% unvalidated$mut,] 
FP_empty_bams_gatk4m2 <- muts_gatk4m2[muts_gatk4m2$T_SAMPLE %like% "_vs_",]
FP_gatk4m2 <- FP_empty_bams_gatk4m2
TN_gatk4m2 <- unvalidated[!(unvalidated$mut %in% muts_gatk4m2$mut),] 
Recall_gatk4m2 <- nrow(TP_gatk4m2)/nrow(validated) *100
Specificity_gatk4m2 <-nrow(TN_gatk4m2) / (nrow(TN_gatk4m2) + nrow(FP_gatk4m2)) *100
Precision_gatk4m2 <-nrow(TP_gatk4m2)/(nrow(TP_gatk4m2)+nrow(FP_gatk4m2)) *100
Accuracy_gatk4m2 <-100*(nrow(TN_gatk4m2)+nrow(TP_gatk4m2))/(nrow(TN_gatk4m2)+nrow(FN_gatk4m2)+nrow(TP_gatk4m2)+nrow(FP_gatk4m2))
F1_Score_gatk4m2 <- 2 * Precision_gatk4m2 * Recall_gatk4m2 / (Precision_gatk4m2 +Recall_gatk4m2)
MCC_gatk4m2 <- (as.numeric(nrow(TP_gatk4m2))*as.numeric(nrow(TN_gatk4m2)) -as.numeric(nrow(FP_gatk4m2))*as.numeric(nrow(FN_gatk4m2))) / sqrt((as.numeric(nrow(TP_gatk4m2))+as.numeric(nrow(FP_gatk4m2)))*(as.numeric(nrow(TP_gatk4m2))+as.numeric(nrow(FN_gatk4m2)))*(as.numeric(nrow(TN_gatk4m2))+as.numeric(nrow(FP_gatk4m2)))*(as.numeric(nrow(TN_gatk4m2))+as.numeric(nrow(FN_gatk4m2))))


TP_lanceta <- validated[validated$mut %in% muts_lancet$mut,] 
TP_lancet <- muts_lancet[muts_lancet$mut %in% validated$mut,] 
FN_lancet <- anti_join(validated, TP_lanceta) 
# FP_TCGA_lancet <- muts_lancet[muts_lancet$mut %in% unvalidated$mut,] 
FP_empty_bams_lancet <- muts_lancet[muts_lancet$T_SAMPLE %like% "_vs_",]
FP_lancet <- FP_empty_bams_lancet
TN_lancet <- unvalidated[!(unvalidated$mut %in% muts_lancet$mut),] 
Recall_lancet <- nrow(TP_lancet)/nrow(validated) *100
Specificity_lancet <-nrow(TN_lancet) / (nrow(TN_lancet) + nrow(FP_lancet)) *100
Precision_lancet <-nrow(TP_lancet)/(nrow(TP_lancet)+nrow(FP_lancet)) *100
Accuracy_lancet <-100*(nrow(TN_lancet)+nrow(TP_lancet))/(nrow(TN_lancet)+nrow(FN_lancet)+nrow(TP_lancet)+nrow(FP_lancet))
F1_Score_lancet <- 2 * Precision_lancet * Recall_lancet / (Precision_lancet +Recall_lancet)
MCC_lancet <- (as.numeric(nrow(TP_lancet))*as.numeric(nrow(TN_lancet)) -as.numeric(nrow(FP_lancet))*as.numeric(nrow(FN_lancet))) / sqrt((as.numeric(nrow(TP_lancet))+as.numeric(nrow(FP_lancet)))*(as.numeric(nrow(TP_lancet))+as.numeric(nrow(FN_lancet)))*(as.numeric(nrow(TN_lancet))+as.numeric(nrow(FP_lancet)))*(as.numeric(nrow(TN_lancet))+as.numeric(nrow(FN_lancet))))


TP_cavemana <- validated[validated$mut %in% muts_caveman$mut,] 
TP_caveman <- muts_caveman[muts_caveman$mut %in% validated$mut,] 
FN_caveman <- anti_join(validated, TP_cavemana) 
# FP_TCGA_caveman <- muts_caveman[muts_caveman$mut %in% unvalidated$mut,] 
FP_empty_bams_caveman <- muts_caveman[muts_caveman$T_SAMPLE %like% "_vs_",]
FP_caveman <- FP_empty_bams_caveman
TN_caveman <- unvalidated[!(unvalidated$mut %in% muts_caveman$mut),] 
Recall_caveman <- nrow(TP_caveman)/nrow(validated) *100
Specificity_caveman <-nrow(TN_caveman) / (nrow(TN_caveman) + nrow(FP_caveman)) *100
Precision_caveman <-nrow(TP_caveman)/(nrow(TP_caveman)+nrow(FP_caveman)) *100
Accuracy_caveman <-100*(nrow(TN_caveman)+nrow(TP_caveman))/(nrow(TN_caveman)+nrow(FN_caveman)+nrow(TP_caveman)+nrow(FP_caveman))
F1_Score_caveman <- 2 * Precision_caveman * Recall_caveman / (Precision_caveman +Recall_caveman)
MCC_caveman <- (as.numeric(nrow(TP_caveman))*as.numeric(nrow(TN_caveman)) -as.numeric(nrow(FP_caveman))*as.numeric(nrow(FN_caveman))) / sqrt((as.numeric(nrow(TP_caveman))+as.numeric(nrow(FP_caveman)))*(as.numeric(nrow(TP_caveman))+as.numeric(nrow(FN_caveman)))*(as.numeric(nrow(TN_caveman))+as.numeric(nrow(FP_caveman)))*(as.numeric(nrow(TN_caveman))+as.numeric(nrow(FN_caveman))))


TP_musea <- validated[validated$mut %in% muts_muse$mut,] 
TP_muse <- muts_muse[muts_muse$mut %in% validated$mut,] 
FN_muse <- anti_join(validated, TP_musea) 
# FP_TCGA_muse <- muts_muse[muts_muse$mut %in% unvalidated$mut,] 
FP_empty_bams_muse <- muts_muse[muts_muse$T_SAMPLE %like% "_vs_",]
FP_muse <- FP_empty_bams_muse
TN_muse <- unvalidated[!(unvalidated$mut %in% muts_muse$mut),] 
Recall_muse <- nrow(TP_muse)/nrow(validated) *100
Specificity_muse <-nrow(TN_muse) / (nrow(TN_muse) + nrow(FP_muse)) *100
Precision_muse <-nrow(TP_muse)/(nrow(TP_muse)+nrow(FP_muse)) *100
Accuracy_muse <-100*(nrow(TN_muse)+nrow(TP_muse))/(nrow(TN_muse)+nrow(FN_muse)+nrow(TP_muse)+nrow(FP_muse))
F1_Score_muse <- 2 * Precision_muse * Recall_muse / (Precision_muse +Recall_muse)
MCC_muse <- (as.numeric(nrow(TP_muse))*as.numeric(nrow(TN_muse)) -as.numeric(nrow(FP_muse))*as.numeric(nrow(FN_muse))) / sqrt((as.numeric(nrow(TP_muse))+as.numeric(nrow(FP_muse)))*(as.numeric(nrow(TP_muse))+as.numeric(nrow(FN_muse)))*(as.numeric(nrow(TN_muse))+as.numeric(nrow(FP_muse)))*(as.numeric(nrow(TN_muse))+as.numeric(nrow(FN_muse))))



TP_strelkaa <- validated[validated$mut %in% muts_strelka$mut,] 
TP_strelka <- muts_strelka[muts_strelka$mut %in% validated$mut,] 
FN_strelka <- anti_join(validated, TP_strelkaa) 
# FP_TCGA_strelka <- muts_strelka[muts_strelka$mut %in% unvalidated$mut,] 
FP_empty_bams_strelka <- muts_strelka[muts_strelka$T_SAMPLE %like% "_vs_",]
FP_strelka <- FP_empty_bams_strelka
TN_strelka <- unvalidated[!(unvalidated$mut %in% muts_strelka$mut),] 
Recall_strelka <- nrow(TP_strelka)/nrow(validated) *100
Specificity_strelka <-nrow(TN_strelka) / (nrow(TN_strelka) + nrow(FP_strelka)) *100
Precision_strelka <-nrow(TP_strelka)/(nrow(TP_strelka)+nrow(FP_strelka)) *100
Accuracy_strelka <-100*(nrow(TN_strelka)+nrow(TP_strelka))/(nrow(TN_strelka)+nrow(FN_strelka)+nrow(TP_strelka)+nrow(FP_strelka))
F1_Score_strelka <- 2 * Precision_strelka * Recall_strelka / (Precision_strelka +Recall_strelka)
MCC_strelka <- (as.numeric(nrow(TP_strelka))*as.numeric(nrow(TN_strelka)) -as.numeric(nrow(FP_strelka))*as.numeric(nrow(FN_strelka))) / sqrt((as.numeric(nrow(TP_strelka))+as.numeric(nrow(FP_strelka)))*(as.numeric(nrow(TP_strelka))+as.numeric(nrow(FN_strelka)))*(as.numeric(nrow(TN_strelka))+as.numeric(nrow(FP_strelka)))*(as.numeric(nrow(TN_strelka))+as.numeric(nrow(FN_strelka))))






row_names <- c("1 or more","2 or more","3 or more","4 or more","5","gatk4m2","lancet","caveman","muse","strelka")

#OJITOOOOOO PONER EL ORDEN CORRECTO EN LAS MÉTRICAS A CONTINUACIÓN DE CADA PROGRAMA POR SEPARADO!!!

Precision_values <- c(Precision_1,Precision_2,Precision_3,Precision_4,Precision_5,Precision_gatk4m2,Precision_lancet,Precision_caveman,Precision_muse,Precision_strelka)
MCC_values <- c(MCC_1,MCC_2,MCC_3,MCC_4,MCC_5,MCC_gatk4m2,MCC_lancet,MCC_caveman,MCC_muse,MCC_strelka)
F1_Score_values <- c(F1_Score_1,F1_Score_2,F1_Score_3,F1_Score_4,F1_Score_5,F1_Score_gatk4m2,F1_Score_lancet,F1_Score_caveman,F1_Score_muse,F1_Score_strelka)
Accuracy_values <- c(Accuracy_1,Accuracy_2,Accuracy_3,Accuracy_4,Accuracy_5,Accuracy_gatk4m2,Accuracy_lancet,Accuracy_caveman,Accuracy_muse,Accuracy_strelka)
Specificity_values <- c(Specificity_1,Specificity_2,Specificity_3,Specificity_4,Specificity_5,Specificity_gatk4m2,Specificity_lancet,Specificity_caveman,Specificity_muse,Specificity_strelka)
Recall_values <- c(Recall_1,Recall_2,Recall_3,Recall_4,Recall_5,Recall_gatk4m2,Recall_lancet,Recall_caveman,Recall_muse,Recall_strelka)

TN_values <- c(nrow(TN_1),nrow(TN_2),nrow(TN_3),nrow(TN_4),nrow(TN_5),nrow(TN_gatk4m2),nrow(TN_lancet),nrow(TN_caveman),nrow(TN_muse),nrow(TN_strelka))
FP_values <- c(nrow(FP_1),nrow(FP_2),nrow(FP_3),nrow(FP_4),nrow(FP_5),nrow(FP_gatk4m2),nrow(FP_lancet),nrow(FP_caveman),nrow(FP_muse),nrow(FP_strelka))
TP_values <- c(nrow(TP_1),nrow(TP_2),nrow(TP_3),nrow(TP_4),nrow(TP_5),nrow(TP_gatk4m2),nrow(TP_lancet),nrow(TP_caveman),nrow(TP_muse),nrow(TP_strelka))
FN_values <- c(nrow(FN_1),nrow(FN_2),nrow(FN_3),nrow(FN_4),nrow(FN_5),nrow(FN_gatk4m2),nrow(FN_lancet),nrow(FN_caveman),nrow(FN_muse),nrow(FN_strelka))


METRICS <- data.table(row_names,Recall_values,Precision_values,Specificity_values,FP_values,TP_values,FN_values, TN_values,Accuracy_values,MCC_values,F1_Score_values)
METRICS

#Tal y como están ahora los datos sale esto, pero deberían equilibrarse el número de FP que a pesar de ser alto es bajo en comparación a los TP. 
#ya sabes que independientemente de la proporción de clases en el test, en el training debe ser del 50%
#    row_names Recall_values Precision_values Specificity_values FP_values TP_values FN_values TN_values Accuracy_values MCC_values F1_Score_values
# 1: 1 or more      98.42687         85.07526           4.622337      7387     42108       673       358        84.04782 0.07769863        91.26533
# 2: 2 or more      96.97529         99.04269          54.740406       401     41487      1294       485        96.11835 0.36883484        97.99809
# 3: 3 or more      93.94124         99.82365          88.871473        71     40189      2592       567        93.86674 0.38363754        96.79315
# 4: 4 or more      86.74412         99.97306          98.456790        10     37110      5671       638        86.91888 0.29313855        92.88995
# 5:         5      68.30602        100.00000         100.000000         0     29222     13559       702        68.81770 0.18336780        81.16884
# 6:   gatk4m2      88.41542         99.36689          70.893720       241     37825      4956       587        88.08274 0.24300834        93.57181
# 7:    lancet      79.30156         98.95867          61.571582       357     33926      8855       572        78.92473 0.14333267        88.04630
# 8:   caveman      89.70337         95.72940          24.747253      1712     38376      4405       563        86.42356 0.10101941        92.61847
# 9:      muse      91.25780         99.20466          66.198704       313     39041      3740       613        90.72689 0.27630032        95.06544
# 10:  strelka      95.71539         88.64355           7.330860      5246     40948      1833       415        85.38665 0.04652071        92.04383

set.seed(500)
vector_sample <- sample(1:nrow(validated), nrow(validated)*0.2,replace = F)
validated2 <- validated[vector_sample,]
typeof(validated2)


TP_1a <- validated2[validated2$mut %in% muts_1$mut,] # 42108 #Orden no importa aquí
TP_1 <- muts_1[muts_1$mut %in% validated2$mut,] 
FN_1 <- anti_join(validated2, TP_1a) # 673 #Orden si importa aquí
# FP_TCGA_1 <- muts_1[muts_1$mut %in% unvalidated$mut,] # 411 
FP_empty_bams_1 <- muts_1[muts_1$T_SAMPLE %like% "_vs_",]
FP_1 <- FP_empty_bams_1
nrow(FP_1) #7798
TN_1 <- unvalidated[!(unvalidated$mut %in% muts_1$mut),] #358
Recall_1 <- nrow(TP_1)/nrow(validated2) *100
#If TN from the normal-normal datasets are not included, specificity values will heavily depend on FP rates
Specificity_1 <-nrow(TN_1) / (nrow(TN_1) + nrow(FP_1)) *100
Precision_1 <-nrow(TP_1)/(nrow(TP_1)+nrow(FP_1)) *100
#The number of FP is 6-fold lower than the number of TP!!
# > nrow(FP_1)
# [1] 7798
# > nrow(TP_1)
# [1] 42108
Accuracy_1 <-100*(nrow(TN_1)+nrow(TP_1))/(nrow(TN_1)+nrow(FN_1)+nrow(TP_1)+nrow(FP_1))
F1_Score_1 <- 2 * Precision_1 * Recall_1 / (Precision_1 +Recall_1)
MCC_1 <- (as.numeric(nrow(TP_1))*as.numeric(nrow(TN_1)) -as.numeric(nrow(FP_1))*as.numeric(nrow(FN_1))) / sqrt((as.numeric(nrow(TP_1))+as.numeric(nrow(FP_1)))*(as.numeric(nrow(TP_1))+as.numeric(nrow(FN_1)))*(as.numeric(nrow(TN_1))+as.numeric(nrow(FP_1)))*(as.numeric(nrow(TN_1))+as.numeric(nrow(FN_1))))


TP_2a <- validated2[validated2$mut %in% muts_2$mut,] 
TP_2 <- muts_2[muts_2$mut %in% validated2$mut,] 
FN_2 <- anti_join(validated2, TP_2a) 
# FP_TCGA_2 <- muts_2[muts_2$mut %in% unvalidated$mut,] 
FP_empty_bams_2 <- muts_2[muts_2$T_SAMPLE %like% "_vs_",]
FP_2 <- FP_empty_bams_2
TN_2 <- unvalidated[!(unvalidated$mut %in% muts_2$mut),] 
Recall_2 <- nrow(TP_2)/nrow(validated2) *100
Specificity_2 <-nrow(TN_2) / (nrow(TN_2) + nrow(FP_2)) *100
Precision_2 <-nrow(TP_2)/(nrow(TP_2)+nrow(FP_2)) *100
Accuracy_2 <-100*(nrow(TN_2)+nrow(TP_2))/(nrow(TN_2)+nrow(FN_2)+nrow(TP_2)+nrow(FP_2))
F1_Score_2 <- 2 * Precision_2 * Recall_2 / (Precision_2 +Recall_2)
MCC_2 <- (as.numeric(nrow(TP_2))*as.numeric(nrow(TN_2)) -as.numeric(nrow(FP_2))*as.numeric(nrow(FN_2))) / sqrt((as.numeric(nrow(TP_2))+as.numeric(nrow(FP_2)))*(as.numeric(nrow(TP_2))+as.numeric(nrow(FN_2)))*(as.numeric(nrow(TN_2))+as.numeric(nrow(FP_2)))*(as.numeric(nrow(TN_2))+as.numeric(nrow(FN_2))))

TP_3a <- validated2[validated2$mut %in% muts_3$mut,] 
TP_3 <- muts_3[muts_3$mut %in% validated2$mut,] 
FN_3 <- anti_join(validated2, TP_3a) 
# FP_TCGA_3 <- muts_3[muts_3$mut %in% unvalidated$mut,] 
FP_empty_bams_3 <- muts_3[muts_3$T_SAMPLE %like% "_vs_",]
FP_3 <- FP_empty_bams_3
TN_3 <- unvalidated[!(unvalidated$mut %in% muts_3$mut),] 
Recall_3 <- nrow(TP_3)/nrow(validated2) *100
Specificity_3 <-nrow(TN_3) / (nrow(TN_3) + nrow(FP_3)) *100
Precision_3 <-nrow(TP_3)/(nrow(TP_3)+nrow(FP_3)) *100
Accuracy_3 <-100*(nrow(TN_3)+nrow(TP_3))/(nrow(TN_3)+nrow(FN_3)+nrow(TP_3)+nrow(FP_3))
F1_Score_3 <- 2 * Precision_3 * Recall_3 / (Precision_3 +Recall_3)
MCC_3 <- (as.numeric(nrow(TP_3))*as.numeric(nrow(TN_3)) -as.numeric(nrow(FP_3))*as.numeric(nrow(FN_3))) / sqrt((as.numeric(nrow(TP_3))+as.numeric(nrow(FP_3)))*(as.numeric(nrow(TP_3))+as.numeric(nrow(FN_3)))*(as.numeric(nrow(TN_3))+as.numeric(nrow(FP_3)))*(as.numeric(nrow(TN_3))+as.numeric(nrow(FN_3))))

TP_4a <- validated2[validated2$mut %in% muts_4$mut,] 
TP_4 <- muts_4[muts_4$mut %in% validated2$mut,] 
FN_4 <- anti_join(validated2, TP_4a) 
# FP_TCGA_4 <- muts_4[muts_4$mut %in% unvalidated$mut,] 
FP_empty_bams_4 <- muts_4[muts_4$T_SAMPLE %like% "_vs_",]
FP_4 <- FP_empty_bams_4
TN_4 <- unvalidated[!(unvalidated$mut %in% muts_4$mut),] 
Recall_4 <- nrow(TP_4)/nrow(validated2) *100
Specificity_4 <-nrow(TN_4) / (nrow(TN_4) + nrow(FP_4)) *100
Precision_4 <-nrow(TP_4)/(nrow(TP_4)+nrow(FP_4)) *100
Accuracy_4 <-100*(nrow(TN_4)+nrow(TP_4))/(nrow(TN_4)+nrow(FN_4)+nrow(TP_4)+nrow(FP_4))
F1_Score_4 <- 2 * Precision_4 * Recall_4 / (Precision_4 +Recall_4)
MCC_4 <- (as.numeric(nrow(TP_4))*as.numeric(nrow(TN_4)) -as.numeric(nrow(FP_4))*as.numeric(nrow(FN_4))) / sqrt((as.numeric(nrow(TP_4))+as.numeric(nrow(FP_4)))*(as.numeric(nrow(TP_4))+as.numeric(nrow(FN_4)))*(as.numeric(nrow(TN_4))+as.numeric(nrow(FP_4)))*(as.numeric(nrow(TN_4))+as.numeric(nrow(FN_4))))

TP_5a <- validated2[validated2$mut %in% muts_5$mut,] 
TP_5 <- muts_5[muts_5$mut %in% validated2$mut,] 
FN_5 <- anti_join(validated2, TP_5a) 
# FP_TCGA_5 <- muts_5[muts_5$mut %in% unvalidated$mut,] 
FP_empty_bams_5 <- muts_5[muts_5$T_SAMPLE %like% "_vs_",]
FP_5 <-FP_empty_bams_5
TN_5 <- unvalidated[!(unvalidated$mut %in% muts_5$mut),] 
Recall_5 <- nrow(TP_5)/nrow(validated2) *100
Specificity_5 <-nrow(TN_5) / (nrow(TN_5) + nrow(FP_5)) *100
Precision_5 <-nrow(TP_5)/(nrow(TP_5)+nrow(FP_5)) *100
Accuracy_5 <-100*(nrow(TN_5)+nrow(TP_5))/(nrow(TN_5)+nrow(FN_5)+nrow(TP_5)+nrow(FP_5))
F1_Score_5 <- 2 * Precision_5 * Recall_5 / (Precision_5 +Recall_5)
MCC_5 <- (as.numeric(nrow(TP_5))*as.numeric(nrow(TN_5)) -as.numeric(nrow(FP_5))*as.numeric(nrow(FN_5))) / sqrt((as.numeric(nrow(TP_5))+as.numeric(nrow(FP_5)))*(as.numeric(nrow(TP_5))+as.numeric(nrow(FN_5)))*(as.numeric(nrow(TN_5))+as.numeric(nrow(FP_5)))*(as.numeric(nrow(TN_5))+as.numeric(nrow(FN_5))))


#Falta la de los programas


TP_gatk4m2a <- validated2[validated2$mut %in% muts_gatk4m2$mut,] 
TP_gatk4m2 <- muts_gatk4m2[muts_gatk4m2$mut %in% validated2$mut,]
FN_gatk4m2 <- anti_join(validated2, TP_gatk4m2a) 
# FP_TCGA_gatk4m2 <- muts_gatk4m2[muts_gatk4m2$mut %in% unvalidated$mut,] 
FP_empty_bams_gatk4m2 <- muts_gatk4m2[muts_gatk4m2$T_SAMPLE %like% "_vs_",]
FP_gatk4m2 <- FP_empty_bams_gatk4m2
TN_gatk4m2 <- unvalidated[!(unvalidated$mut %in% muts_gatk4m2$mut),] 
Recall_gatk4m2 <- nrow(TP_gatk4m2)/nrow(validated2) *100
Specificity_gatk4m2 <-nrow(TN_gatk4m2) / (nrow(TN_gatk4m2) + nrow(FP_gatk4m2)) *100
Precision_gatk4m2 <-nrow(TP_gatk4m2)/(nrow(TP_gatk4m2)+nrow(FP_gatk4m2)) *100
Accuracy_gatk4m2 <-100*(nrow(TN_gatk4m2)+nrow(TP_gatk4m2))/(nrow(TN_gatk4m2)+nrow(FN_gatk4m2)+nrow(TP_gatk4m2)+nrow(FP_gatk4m2))
F1_Score_gatk4m2 <- 2 * Precision_gatk4m2 * Recall_gatk4m2 / (Precision_gatk4m2 +Recall_gatk4m2)
MCC_gatk4m2 <- (as.numeric(nrow(TP_gatk4m2))*as.numeric(nrow(TN_gatk4m2)) -as.numeric(nrow(FP_gatk4m2))*as.numeric(nrow(FN_gatk4m2))) / sqrt((as.numeric(nrow(TP_gatk4m2))+as.numeric(nrow(FP_gatk4m2)))*(as.numeric(nrow(TP_gatk4m2))+as.numeric(nrow(FN_gatk4m2)))*(as.numeric(nrow(TN_gatk4m2))+as.numeric(nrow(FP_gatk4m2)))*(as.numeric(nrow(TN_gatk4m2))+as.numeric(nrow(FN_gatk4m2))))


TP_lanceta <- validated2[validated2$mut %in% muts_lancet$mut,] 
TP_lancet <- muts_lancet[muts_lancet$mut %in% validated2$mut,]
FN_lancet <- anti_join(validated2, TP_lanceta) 
# FP_TCGA_lancet <- muts_lancet[muts_lancet$mut %in% unvalidated$mut,] 
FP_empty_bams_lancet <- muts_lancet[muts_lancet$T_SAMPLE %like% "_vs_",]
FP_lancet <- FP_empty_bams_lancet
TN_lancet <- unvalidated[!(unvalidated$mut %in% muts_lancet$mut),] 
Recall_lancet <- nrow(TP_lancet)/nrow(validated2) *100
Specificity_lancet <-nrow(TN_lancet) / (nrow(TN_lancet) + nrow(FP_lancet)) *100
Precision_lancet <-nrow(TP_lancet)/(nrow(TP_lancet)+nrow(FP_lancet)) *100
Accuracy_lancet <-100*(nrow(TN_lancet)+nrow(TP_lancet))/(nrow(TN_lancet)+nrow(FN_lancet)+nrow(TP_lancet)+nrow(FP_lancet))
F1_Score_lancet <- 2 * Precision_lancet * Recall_lancet / (Precision_lancet +Recall_lancet)
MCC_lancet <- (as.numeric(nrow(TP_lancet))*as.numeric(nrow(TN_lancet)) -as.numeric(nrow(FP_lancet))*as.numeric(nrow(FN_lancet))) / sqrt((as.numeric(nrow(TP_lancet))+as.numeric(nrow(FP_lancet)))*(as.numeric(nrow(TP_lancet))+as.numeric(nrow(FN_lancet)))*(as.numeric(nrow(TN_lancet))+as.numeric(nrow(FP_lancet)))*(as.numeric(nrow(TN_lancet))+as.numeric(nrow(FN_lancet))))


TP_cavemana <- validated2[validated2$mut %in% muts_caveman$mut,] 
TP_caveman <- muts_caveman[muts_caveman$mut %in% validated2$mut,]
FN_caveman <- anti_join(validated2, TP_cavemana) 
# FP_TCGA_caveman <- muts_caveman[muts_caveman$mut %in% unvalidated$mut,] 
FP_empty_bams_caveman <- muts_caveman[muts_caveman$T_SAMPLE %like% "_vs_",]
FP_caveman <- FP_empty_bams_caveman
TN_caveman <- unvalidated[!(unvalidated$mut %in% muts_caveman$mut),] 
Recall_caveman <- nrow(TP_caveman)/nrow(validated2) *100
Specificity_caveman <-nrow(TN_caveman) / (nrow(TN_caveman) + nrow(FP_caveman)) *100
Precision_caveman <-nrow(TP_caveman)/(nrow(TP_caveman)+nrow(FP_caveman)) *100
Accuracy_caveman <-100*(nrow(TN_caveman)+nrow(TP_caveman))/(nrow(TN_caveman)+nrow(FN_caveman)+nrow(TP_caveman)+nrow(FP_caveman))
F1_Score_caveman <- 2 * Precision_caveman * Recall_caveman / (Precision_caveman +Recall_caveman)
MCC_caveman <- (as.numeric(nrow(TP_caveman))*as.numeric(nrow(TN_caveman)) -as.numeric(nrow(FP_caveman))*as.numeric(nrow(FN_caveman))) / sqrt((as.numeric(nrow(TP_caveman))+as.numeric(nrow(FP_caveman)))*(as.numeric(nrow(TP_caveman))+as.numeric(nrow(FN_caveman)))*(as.numeric(nrow(TN_caveman))+as.numeric(nrow(FP_caveman)))*(as.numeric(nrow(TN_caveman))+as.numeric(nrow(FN_caveman))))


TP_musea <- validated2[validated2$mut %in% muts_muse$mut,] 
TP_muse <- muts_muse[muts_muse$mut %in% validated2$mut,]
FN_muse <- anti_join(validated2, TP_musea) 
# FP_TCGA_muse <- muts_muse[muts_muse$mut %in% unvalidated$mut,] 
FP_empty_bams_muse <- muts_muse[muts_muse$T_SAMPLE %like% "_vs_",]
FP_muse <- FP_empty_bams_muse
TN_muse <- unvalidated[!(unvalidated$mut %in% muts_muse$mut),] 
Recall_muse <- nrow(TP_muse)/nrow(validated2) *100
Specificity_muse <-nrow(TN_muse) / (nrow(TN_muse) + nrow(FP_muse)) *100
Precision_muse <-nrow(TP_muse)/(nrow(TP_muse)+nrow(FP_muse)) *100
Accuracy_muse <-100*(nrow(TN_muse)+nrow(TP_muse))/(nrow(TN_muse)+nrow(FN_muse)+nrow(TP_muse)+nrow(FP_muse))
F1_Score_muse <- 2 * Precision_muse * Recall_muse / (Precision_muse +Recall_muse)
MCC_muse <- (as.numeric(nrow(TP_muse))*as.numeric(nrow(TN_muse)) -as.numeric(nrow(FP_muse))*as.numeric(nrow(FN_muse))) / sqrt((as.numeric(nrow(TP_muse))+as.numeric(nrow(FP_muse)))*(as.numeric(nrow(TP_muse))+as.numeric(nrow(FN_muse)))*(as.numeric(nrow(TN_muse))+as.numeric(nrow(FP_muse)))*(as.numeric(nrow(TN_muse))+as.numeric(nrow(FN_muse))))



TP_strelkaa <- validated2[validated2$mut %in% muts_strelka$mut,] 
TP_strelka <- muts_strelka[muts_strelka$mut %in% validated2$mut,]
FN_strelka <- anti_join(validated2, TP_strelkaa) 
# FP_TCGA_strelka <- muts_strelka[muts_strelka$mut %in% unvalidated$mut,] 
FP_empty_bams_strelka <- muts_strelka[muts_strelka$T_SAMPLE %like% "_vs_",]
FP_strelka <- FP_empty_bams_strelka
TN_strelka <- unvalidated[!(unvalidated$mut %in% muts_strelka$mut),] 
Recall_strelka <- nrow(TP_strelka)/nrow(validated2) *100
Specificity_strelka <-nrow(TN_strelka) / (nrow(TN_strelka) + nrow(FP_strelka)) *100
Precision_strelka <-nrow(TP_strelka)/(nrow(TP_strelka)+nrow(FP_strelka)) *100
Accuracy_strelka <-100*(nrow(TN_strelka)+nrow(TP_strelka))/(nrow(TN_strelka)+nrow(FN_strelka)+nrow(TP_strelka)+nrow(FP_strelka))
F1_Score_strelka <- 2 * Precision_strelka * Recall_strelka / (Precision_strelka +Recall_strelka)
MCC_strelka <- (as.numeric(nrow(TP_strelka))*as.numeric(nrow(TN_strelka)) -as.numeric(nrow(FP_strelka))*as.numeric(nrow(FN_strelka))) / sqrt((as.numeric(nrow(TP_strelka))+as.numeric(nrow(FP_strelka)))*(as.numeric(nrow(TP_strelka))+as.numeric(nrow(FN_strelka)))*(as.numeric(nrow(TN_strelka))+as.numeric(nrow(FP_strelka)))*(as.numeric(nrow(TN_strelka))+as.numeric(nrow(FN_strelka))))






row_names <- c("1 or more","2 or more","3 or more","4 or more","5","gatk4m2","lancet","caveman","muse","strelka")

#OJITOOOOOO PONER EL ORDEN CORRECTO EN LAS MÉTRICAS A CONTINUACIÓN DE CADA PROGRAMA POR SEPARADO!!!

Precision_values <- c(Precision_1,Precision_2,Precision_3,Precision_4,Precision_5,Precision_gatk4m2,Precision_lancet,Precision_caveman,Precision_muse,Precision_strelka)
MCC_values <- c(MCC_1,MCC_2,MCC_3,MCC_4,MCC_5,MCC_gatk4m2,MCC_lancet,MCC_caveman,MCC_muse,MCC_strelka)
F1_Score_values <- c(F1_Score_1,F1_Score_2,F1_Score_3,F1_Score_4,F1_Score_5,F1_Score_gatk4m2,F1_Score_lancet,F1_Score_caveman,F1_Score_muse,F1_Score_strelka)
Accuracy_values <- c(Accuracy_1,Accuracy_2,Accuracy_3,Accuracy_4,Accuracy_5,Accuracy_gatk4m2,Accuracy_lancet,Accuracy_caveman,Accuracy_muse,Accuracy_strelka)
Specificity_values <- c(Specificity_1,Specificity_2,Specificity_3,Specificity_4,Specificity_5,Specificity_gatk4m2,Specificity_lancet,Specificity_caveman,Specificity_muse,Specificity_strelka)
Recall_values <- c(Recall_1,Recall_2,Recall_3,Recall_4,Recall_5,Recall_gatk4m2,Recall_lancet,Recall_caveman,Recall_muse,Recall_strelka)

TN_values <- c(nrow(TN_1),nrow(TN_2),nrow(TN_3),nrow(TN_4),nrow(TN_5),nrow(TN_gatk4m2),nrow(TN_lancet),nrow(TN_caveman),nrow(TN_muse),nrow(TN_strelka))
FP_values <- c(nrow(FP_1),nrow(FP_2),nrow(FP_3),nrow(FP_4),nrow(FP_5),nrow(FP_gatk4m2),nrow(FP_lancet),nrow(FP_caveman),nrow(FP_muse),nrow(FP_strelka))
TP_values <- c(nrow(TP_1),nrow(TP_2),nrow(TP_3),nrow(TP_4),nrow(TP_5),nrow(TP_gatk4m2),nrow(TP_lancet),nrow(TP_caveman),nrow(TP_muse),nrow(TP_strelka))
FN_values <- c(nrow(FN_1),nrow(FN_2),nrow(FN_3),nrow(FN_4),nrow(FN_5),nrow(FN_gatk4m2),nrow(FN_lancet),nrow(FN_caveman),nrow(FN_muse),nrow(FN_strelka))


METRICS <- data.table(row_names,Recall_values,Precision_values,Specificity_values,FP_values,TP_values,FN_values, TN_values,Accuracy_values,MCC_values,F1_Score_values)
METRICS
#    row_names Recall_values Precision_values Specificity_values FP_values TP_values FN_values TN_values Accuracy_values MCC_values F1_Score_values
# 1: 1 or more      98.60916         51.96797           4.389407      7798      8437       119       358        52.62685 0.09001250        68.06502
# 2: 2 or more      96.90276         92.36854          41.452991       685      8291       265       485        90.23237 0.46772086        94.58134
# 3: 3 or more      94.07433         96.71954          67.500000       273      8049       507       567        91.69860 0.55215214        95.37860
# 4: 4 or more      87.36559         98.14863          81.899872       141      7475      1081       638        86.90948 0.49421902        92.44373
# 5:         5      68.84058         98.87527          91.287386        67      5890      2666       702        70.69169 0.34433029        81.16861
# 6:   gatk4m2      88.75643         94.72371          58.118812       423      7594       962       587        85.52164 0.39102807        91.64303
# 7:    lancet      79.88546         92.50237          50.799290       554      6835      1721       572        76.50279 0.23138346        85.73220
# 8:   caveman      90.07714         80.07273          22.692463      1918      7707       849       563        74.92978 0.15959100        84.78082
# 9:      muse      91.08228         94.32341          56.654344       469      7793       763       613        87.21726 0.43077344        92.67452
# 10:  strelka      95.99112         59.45848           6.899418      5600      8213       343       415        59.21351 0.06408448        73.43198


#OJITO AVERIGUAR PORQUE CAMBIAN LOS TP!!!!!!! Será por la seleccíón aleatoria??????? fijar una seed!!!!!!!!!!!! para garantizar reproducibilidad!!
#YA ESTA FIJADA Y SIEMRE DA ESTO A PESAR DE QUE REPITA EL PASO DE SAMPLE, en principio la hipótesis es que el cambio de los números en los TP se debía a haber repetido el paso de sample si haber fijado una misma seed
#    row_names Recall_values Precision_values Specificity_values FP_values TP_values FN_values TN_values Accuracy_values MCC_values F1_Score_values
# 1: 1 or more      98.45722         53.27936           4.622337      7387      8424       132       358        53.87400 0.09006521        69.14269
# 2: 2 or more      97.14820         95.39768          54.740406       401      8312       244       485        93.16882 0.56686107        96.26498
# 3: 3 or more      94.03927         99.12529          88.871473        71      8046       510       567        93.68066 0.65516534        96.51532
# 4: 4 or more      86.16176         99.86454          98.456790        10      7372      1184       638        87.02738 0.54328442        92.50847
# 5:         5      68.04582        100.00000         100.000000         0      5822      2734       702        70.46878 0.37285733        80.98484
# 6:   gatk4m2      88.35905         96.91065          70.893720       241      7560       996       587        86.81799 0.44879043        92.43749
# 7:    lancet      78.55306         94.95620          61.571582       357      6721      1835       572        76.88983 0.27406976        85.97928
# 8:   caveman      89.65638         81.75424          24.747253      1712      7671       885       563        76.02253 0.17240229        85.52316
# 9:      muse      91.70173         96.16375          66.198704       313      7846       710       613        89.21114 0.49604869        93.87975
# 10:  strelka      95.58205         60.92074           7.330860      5246      8178       378       415        60.44172 0.06213426        74.41310
# 
