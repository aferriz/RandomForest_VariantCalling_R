#Nadaaa no hay golden datset de aml31 para indels es solo snvs, por tanto me toca primero adaptar el script para snv!!!!!!!!!!
#es deir me toca primero entrenar al modelo para que funcione con snvs, y ver que tal va y luego probarlo con AML31



#Calcular las métricas para cada programa con este dataset, y luego pasarle el RF y ver que tal son las métricas. 
#Supongo que habra que hacer trasnformacion de datos? Lo mismo no pq en este caso el golden es el que es y los FP son todo lo que no esten ahi. 


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



AML31_snvs_tab <- "/media/computer/HDD/MOUNT/AML31/snvs.PASS.ANNOT.GRCh37.p13.RefSeq.tab"
AML31_snvs_vc <- fread(AML31_snvs_tab, header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))
nSample <- "SRR2251868_DS"
tSample <- "SRR2415196_DS"
colnames(AML31_snvs_vc) <- gsub(paste0(tSample, "."), "T_", colnames(AML31_snvs_vc))
colnames(AML31_snvs_vc) <- gsub(paste0(nSample, "."), "N_", colnames(AML31_snvs_vc))

AML31_snvs_vc$T_SAMPLE <- "AML31"
AML31_snvs_vc$mut <- paste0(AML31_snvs_vc$CHROM,"_",AML31_snvs_vc$POSITION,"_",AML31_snvs_vc$REF,"_",AML31_snvs_vc$ALT,"_",AML31_snvs_vc$T_SAMPLE)


#No se si es necesario que tengan el mismo numero de columnas pero por si acaso lo voy a hacer. Simplemente es añadir a todos el nomrbe del fichero y al golden tb y ya. 



filePlatinum_AML31 <- "/media/computer/HDD/MOUNT/AML31/snvs-platinum.venn" #Empezar con este que es el mejor
# fileGolden_AML31 <- "/media/computer/HDD/MOUNT/AML31/snvs-golden.venn"
# golden_AML31 <- fread(fileGolden_AML31,header=FALSE)
platinum_AML31 <- fread(filePlatinum_AML31,header=FALSE)
colnames(platinum_AML31) <- "GOLD"
platinum_AML31$GOLD <- paste0(platinum_AML31$GOLD,"_","AML31")
head(platinum_AML31)


#Aquí las mismas que en el último script
#OJITO QUE AQUI SAMPLE NO HAY, POR ESO LO HE METIDO ANTES TANTO EN EL TAB COMO EN EL GOLDEN COMO _AML31 A TODOS IGUAL PARA NO TENER QUE MODIFICAR NADA MÁS
selectcols_snv <- c("CHROM", "POSITION", "REF", "ALT", "T_NPROGS", "T_PROGS", "N_DP_avg", "N_DP_std", "N_DP_REF_avg", 
                    "N_DP_ALT_avg", "N_VAF_avg", "N_DP_REF_cav", "N_DP_ALT_cav", "N_DP_REF_g4m2", "N_DP_ALT_g4m2", 
                    "N_MBQ_g4m2", "N_MMQ_g4m2", "N_DP_REF_str", "N_DP_ALT_str", "N_DP_REF_mus", "N_DP_ALT_mus", 
                    "N_BQ_REF_mus", "N_BQ_ALT_mus", "N_DP_REF_lan", "N_DP_ALT_lan", "T_NPROGS", "T_PROGS", "T_DP_avg", 
                    "T_DP_std", "T_DP_REF_avg", "T_DP_ALT_avg", "T_VAF_avg", "T_DP_REF_cav", "T_DP_ALT_cav", "T_VAF_cav", 
                    "T_ASMD_cav", "T_CLPM_cav", "T_DP_REF_g4m2", "T_DP_ALT_g4m2", "T_VAF_g4m2", "T_MBQ_g4m2", "T_MMQ_g4m2", 
                    "T_DP_REF_str", "T_DP_ALT_str", "T_VAF_str", "T_MQ_str", "T_MQ0_str", "T_SomaticEVS_str", "T_DP_REF_mus", 
                    "T_DP_ALT_mus", "T_VAF_mus", "T_BQ_REF_mus", "T_BQ_ALT_mus", "T_DP_REF_lan", "T_DP_ALT_lan", "T_VAF_lan", 
                    "T_FETS_lan", "T_SB_lan", "T_SAMPLE", "mut")


DNA1 <- as.data.frame(AML31_snvs_vc)
muts2 <- DNA1[,names(DNA1) %in% selectcols_snv]
muts2$T_PROGS <- as.factor(muts2$T_PROGS)
muts2$GOLD <- ifelse(muts2$mut %in% platinum_AML31$GOLD, TRUE, FALSE)
muts2[is.na(muts2)] <- -1
#OJITO QUE AQUI LA ESTRATEGIA CAMBIA Y YA ES TODO LO QUE NO ESTE EN EL GOLDEN ES FP

summary(muts2)

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


TP_1a <- platinum_AML31[platinum_AML31$GOLD %in% muts_1$mut,] 
TP_1 <- muts_1[muts_1$GOLD==TRUE,] 
FN_1 <- anti_join(platinum_AML31, TP_1a) 
FP_1 <- muts_1[muts_1$GOLD==FALSE,]
# TN_1 <- unvalidated[!(unvalidated$mut %in% muts_1$mut),]
TN_1 <- (table(factor(c(1:1000)))) #ESTO NO SE QUE VALOR METERLE PQ DEBERIA SER TODO EL GENOMA
Recall_1 <- nrow(TP_1)/nrow(platinum_AML31) *100
#If TN from the normal-normal datasets are not included, specificity values will heavily depend on FP rates
Specificity_1 <-nrow(TN_1) / (nrow(TN_1) + nrow(FP_1)) *100
Precision_1 <-nrow(TP_1)/(nrow(TP_1)+nrow(FP_1)) *100
Accuracy_1 <-100*(nrow(TN_1)+nrow(TP_1))/(nrow(TN_1)+nrow(FN_1)+nrow(TP_1)+nrow(FP_1))
F1_Score_1 <- 2 * Precision_1 * Recall_1 / (Precision_1 +Recall_1)
MCC_1 <- (as.numeric(nrow(TP_1))*as.numeric(nrow(TN_1)) -as.numeric(nrow(FP_1))*as.numeric(nrow(FN_1))) / sqrt((as.numeric(nrow(TP_1))+as.numeric(nrow(FP_1)))*(as.numeric(nrow(TP_1))+as.numeric(nrow(FN_1)))*(as.numeric(nrow(TN_1))+as.numeric(nrow(FP_1)))*(as.numeric(nrow(TN_1))+as.numeric(nrow(FN_1))))

TP_2a <- platinum_AML31[platinum_AML31$GOLD %in% muts_2$mut,] 
TP_2 <- muts_2[muts_2$GOLD==TRUE,] 
FN_2 <- anti_join(platinum_AML31, TP_2a) 
FP_2 <- muts_2[muts_2$GOLD==FALSE,]
# TN_2 <- unvalidated[!(unvalidated$mut %in% muts_2$mut),]
TN_2 <- (table(factor(c(1:1000)))) #ESTO NO SE QUE VALOR METERLE PQ DEBERIA SER TODO EL GENOMA
Recall_2 <- nrow(TP_2)/nrow(platinum_AML31) *100
#If TN from the normal-normal datasets are not included, specificity values will heavily depend on FP rates
Specificity_2 <-nrow(TN_2) / (nrow(TN_2) + nrow(FP_2)) *100
Precision_2 <-nrow(TP_2)/(nrow(TP_2)+nrow(FP_2)) *100
Accuracy_2 <-100*(nrow(TN_2)+nrow(TP_2))/(nrow(TN_2)+nrow(FN_2)+nrow(TP_2)+nrow(FP_2))
F1_Score_2 <- 2 * Precision_2 * Recall_2 / (Precision_2 +Recall_2)
MCC_2 <- (as.numeric(nrow(TP_2))*as.numeric(nrow(TN_2)) -as.numeric(nrow(FP_2))*as.numeric(nrow(FN_2))) / sqrt((as.numeric(nrow(TP_2))+as.numeric(nrow(FP_2)))*(as.numeric(nrow(TP_2))+as.numeric(nrow(FN_2)))*(as.numeric(nrow(TN_2))+as.numeric(nrow(FP_2)))*(as.numeric(nrow(TN_2))+as.numeric(nrow(FN_2))))

TP_3a <- platinum_AML31[platinum_AML31$GOLD %in% muts_3$mut,] 
TP_3 <- muts_3[muts_3$GOLD==TRUE,] 
FN_3 <- anti_join(platinum_AML31, TP_3a) 
FP_3 <- muts_3[muts_3$GOLD==FALSE,]
# TN_3 <- unvalidated[!(unvalidated$mut %in% muts_3$mut),]
TN_3 <- (table(factor(c(1:1000)))) #ESTO NO SE QUE VALOR METERLE PQ DEBERIA SER TODO EL GENOMA
Recall_3 <- nrow(TP_3)/nrow(platinum_AML31) *100
#If TN from the normal-normal datasets are not included, specificity values will heavily depend on FP rates
Specificity_3 <-nrow(TN_3) / (nrow(TN_3) + nrow(FP_3)) *100
Precision_3 <-nrow(TP_3)/(nrow(TP_3)+nrow(FP_3)) *100
Accuracy_3 <-100*(nrow(TN_3)+nrow(TP_3))/(nrow(TN_3)+nrow(FN_3)+nrow(TP_3)+nrow(FP_3))
F1_Score_3 <- 2 * Precision_3 * Recall_3 / (Precision_3 +Recall_3)
MCC_3 <- (as.numeric(nrow(TP_3))*as.numeric(nrow(TN_3)) -as.numeric(nrow(FP_3))*as.numeric(nrow(FN_3))) / sqrt((as.numeric(nrow(TP_3))+as.numeric(nrow(FP_3)))*(as.numeric(nrow(TP_3))+as.numeric(nrow(FN_3)))*(as.numeric(nrow(TN_3))+as.numeric(nrow(FP_3)))*(as.numeric(nrow(TN_3))+as.numeric(nrow(FN_3))))

TP_4a <- platinum_AML31[platinum_AML31$GOLD %in% muts_4$mut,] 
TP_4 <- muts_4[muts_4$GOLD==TRUE,] 
FN_4 <- anti_join(platinum_AML31, TP_4a) 
FP_4 <- muts_4[muts_4$GOLD==FALSE,]
# TN_4 <- unvalidated[!(unvalidated$mut %in% muts_4$mut),]
TN_4 <- (table(factor(c(1:1000)))) #ESTO NO SE QUE VALOR METERLE PQ DEBERIA SER TODO EL GENOMA
Recall_4 <- nrow(TP_4)/nrow(platinum_AML31) *100
#If TN from the normal-normal datasets are not included, specificity values will heavily depend on FP rates
Specificity_4 <-nrow(TN_4) / (nrow(TN_4) + nrow(FP_4)) *100
Precision_4 <-nrow(TP_4)/(nrow(TP_4)+nrow(FP_4)) *100
Accuracy_4 <-100*(nrow(TN_4)+nrow(TP_4))/(nrow(TN_4)+nrow(FN_4)+nrow(TP_4)+nrow(FP_4))
F1_Score_4 <- 2 * Precision_4 * Recall_4 / (Precision_4 +Recall_4)
MCC_4 <- (as.numeric(nrow(TP_4))*as.numeric(nrow(TN_4)) -as.numeric(nrow(FP_4))*as.numeric(nrow(FN_4))) / sqrt((as.numeric(nrow(TP_4))+as.numeric(nrow(FP_4)))*(as.numeric(nrow(TP_4))+as.numeric(nrow(FN_4)))*(as.numeric(nrow(TN_4))+as.numeric(nrow(FP_4)))*(as.numeric(nrow(TN_4))+as.numeric(nrow(FN_4))))

TP_5a <- platinum_AML31[platinum_AML31$GOLD %in% muts_5$mut,] 
TP_5 <- muts_5[muts_5$GOLD==TRUE,] 
FN_5 <- anti_join(platinum_AML31, TP_5a) 
FP_5 <- muts_5[muts_5$GOLD==FALSE,]
# TN_5 <- unvalidated[!(unvalidated$mut %in% muts_5$mut),]
TN_5 <- (table(factor(c(1:1000)))) #ESTO NO SE QUE VALOR METERLE PQ DEBERIA SER TODO EL GENOMA
Recall_5 <- nrow(TP_5)/nrow(platinum_AML31) *100
#If TN from the normal-normal datasets are not included, specificity values will heavily depend on FP rates
Specificity_5 <-nrow(TN_5) / (nrow(TN_5) + nrow(FP_5)) *100
Precision_5 <-nrow(TP_5)/(nrow(TP_5)+nrow(FP_5)) *100
Accuracy_5 <-100*(nrow(TN_5)+nrow(TP_5))/(nrow(TN_5)+nrow(FN_5)+nrow(TP_5)+nrow(FP_5))
F1_Score_5 <- 2 * Precision_5 * Recall_5 / (Precision_5 +Recall_5)
MCC_5 <- (as.numeric(nrow(TP_5))*as.numeric(nrow(TN_5)) -as.numeric(nrow(FP_5))*as.numeric(nrow(FN_5))) / sqrt((as.numeric(nrow(TP_5))+as.numeric(nrow(FP_5)))*(as.numeric(nrow(TP_5))+as.numeric(nrow(FN_5)))*(as.numeric(nrow(TN_5))+as.numeric(nrow(FP_5)))*(as.numeric(nrow(TN_5))+as.numeric(nrow(FN_5))))

TP_gatk4m2a <- platinum_AML31[platinum_AML31$GOLD %in% muts_gatk4m2$mut,] 
TP_gatk4m2 <- muts_gatk4m2[muts_gatk4m2$GOLD==TRUE,] 
FN_gatk4m2 <- anti_join(platinum_AML31, TP_gatk4m2a) 
FP_gatk4m2 <- muts_gatk4m2[muts_gatk4m2$GOLD==FALSE,]
# TN_gatk4m2 <- unvalidated[!(unvalidated$mut %in% muts_gatk4m2$mut),]
TN_gatk4m2 <- (table(factor(c(1:1000)))) #ESTO NO SE QUE VALOR METERLE PQ DEBERIA SER TODO EL GENOMA
Recall_gatk4m2 <- nrow(TP_gatk4m2)/nrow(platinum_AML31) *100
#If TN from the normal-normal datasets are not included, specificity values will heavily depend on FP rates
Specificity_gatk4m2 <-nrow(TN_gatk4m2) / (nrow(TN_gatk4m2) + nrow(FP_gatk4m2)) *100
Precision_gatk4m2 <-nrow(TP_gatk4m2)/(nrow(TP_gatk4m2)+nrow(FP_gatk4m2)) *100
Accuracy_gatk4m2 <-100*(nrow(TN_gatk4m2)+nrow(TP_gatk4m2))/(nrow(TN_gatk4m2)+nrow(FN_gatk4m2)+nrow(TP_gatk4m2)+nrow(FP_gatk4m2))
F1_Score_gatk4m2 <- 2 * Precision_gatk4m2 * Recall_gatk4m2 / (Precision_gatk4m2 +Recall_gatk4m2)
MCC_gatk4m2 <- (as.numeric(nrow(TP_gatk4m2))*as.numeric(nrow(TN_gatk4m2)) -as.numeric(nrow(FP_gatk4m2))*as.numeric(nrow(FN_gatk4m2))) / sqrt((as.numeric(nrow(TP_gatk4m2))+as.numeric(nrow(FP_gatk4m2)))*(as.numeric(nrow(TP_gatk4m2))+as.numeric(nrow(FN_gatk4m2)))*(as.numeric(nrow(TN_gatk4m2))+as.numeric(nrow(FP_gatk4m2)))*(as.numeric(nrow(TN_gatk4m2))+as.numeric(nrow(FN_gatk4m2))))

TP_lanceta <- platinum_AML31[platinum_AML31$GOLD %in% muts_lancet$mut,] 
TP_lancet <- muts_lancet[muts_lancet$GOLD==TRUE,] 
FN_lancet <- anti_join(platinum_AML31, TP_lanceta) 
FP_lancet <- muts_lancet[muts_lancet$GOLD==FALSE,]
# TN_lancet <- unvalidated[!(unvalidated$mut %in% muts_lancet$mut),]
TN_lancet <- (table(factor(c(1:1000)))) #ESTO NO SE QUE VALOR METERLE PQ DEBERIA SER TODO EL GENOMA
Recall_lancet <- nrow(TP_lancet)/nrow(platinum_AML31) *100
#If TN from the normal-normal datasets are not included, specificity values will heavily depend on FP rates
Specificity_lancet <-nrow(TN_lancet) / (nrow(TN_lancet) + nrow(FP_lancet)) *100
Precision_lancet <-nrow(TP_lancet)/(nrow(TP_lancet)+nrow(FP_lancet)) *100
Accuracy_lancet <-100*(nrow(TN_lancet)+nrow(TP_lancet))/(nrow(TN_lancet)+nrow(FN_lancet)+nrow(TP_lancet)+nrow(FP_lancet))
F1_Score_lancet <- 2 * Precision_lancet * Recall_lancet / (Precision_lancet +Recall_lancet)
MCC_lancet <- (as.numeric(nrow(TP_lancet))*as.numeric(nrow(TN_lancet)) -as.numeric(nrow(FP_lancet))*as.numeric(nrow(FN_lancet))) / sqrt((as.numeric(nrow(TP_lancet))+as.numeric(nrow(FP_lancet)))*(as.numeric(nrow(TP_lancet))+as.numeric(nrow(FN_lancet)))*(as.numeric(nrow(TN_lancet))+as.numeric(nrow(FP_lancet)))*(as.numeric(nrow(TN_lancet))+as.numeric(nrow(FN_lancet))))

TP_cavemana <- platinum_AML31[platinum_AML31$GOLD %in% muts_caveman$mut,] 
TP_caveman <- muts_caveman[muts_caveman$GOLD==TRUE,] 
FN_caveman <- anti_join(platinum_AML31, TP_cavemana) 
FP_caveman <- muts_caveman[muts_caveman$GOLD==FALSE,]
# TN_caveman <- unvalidated[!(unvalidated$mut %in% muts_caveman$mut),]
TN_caveman <- (table(factor(c(1:1000)))) #ESTO NO SE QUE VALOR METERLE PQ DEBERIA SER TODO EL GENOMA
Recall_caveman <- nrow(TP_caveman)/nrow(platinum_AML31) *100
#If TN from the normal-normal datasets are not included, specificity values will heavily depend on FP rates
Specificity_caveman <-nrow(TN_caveman) / (nrow(TN_caveman) + nrow(FP_caveman)) *100
Precision_caveman <-nrow(TP_caveman)/(nrow(TP_caveman)+nrow(FP_caveman)) *100
Accuracy_caveman <-100*(nrow(TN_caveman)+nrow(TP_caveman))/(nrow(TN_caveman)+nrow(FN_caveman)+nrow(TP_caveman)+nrow(FP_caveman))
F1_Score_caveman <- 2 * Precision_caveman * Recall_caveman / (Precision_caveman +Recall_caveman)
MCC_caveman <- (as.numeric(nrow(TP_caveman))*as.numeric(nrow(TN_caveman)) -as.numeric(nrow(FP_caveman))*as.numeric(nrow(FN_caveman))) / sqrt((as.numeric(nrow(TP_caveman))+as.numeric(nrow(FP_caveman)))*(as.numeric(nrow(TP_caveman))+as.numeric(nrow(FN_caveman)))*(as.numeric(nrow(TN_caveman))+as.numeric(nrow(FP_caveman)))*(as.numeric(nrow(TN_caveman))+as.numeric(nrow(FN_caveman))))

TP_musea <- platinum_AML31[platinum_AML31$GOLD %in% muts_muse$mut,] 
TP_muse <- muts_muse[muts_muse$GOLD==TRUE,] 
FN_muse <- anti_join(platinum_AML31, TP_musea) 
FP_muse <- muts_muse[muts_muse$GOLD==FALSE,]
# TN_muse <- unvalidated[!(unvalidated$mut %in% muts_muse$mut),]
TN_muse <- (table(factor(c(1:1000)))) #ESTO NO SE QUE VALOR METERLE PQ DEBERIA SER TODO EL GENOMA
Recall_muse <- nrow(TP_muse)/nrow(platinum_AML31) *100
#If TN from the normal-normal datasets are not included, specificity values will heavily depend on FP rates
Specificity_muse <-nrow(TN_muse) / (nrow(TN_muse) + nrow(FP_muse)) *100
Precision_muse <-nrow(TP_muse)/(nrow(TP_muse)+nrow(FP_muse)) *100
Accuracy_muse <-100*(nrow(TN_muse)+nrow(TP_muse))/(nrow(TN_muse)+nrow(FN_muse)+nrow(TP_muse)+nrow(FP_muse))
F1_Score_muse <- 2 * Precision_muse * Recall_muse / (Precision_muse +Recall_muse)
MCC_muse <- (as.numeric(nrow(TP_muse))*as.numeric(nrow(TN_muse)) -as.numeric(nrow(FP_muse))*as.numeric(nrow(FN_muse))) / sqrt((as.numeric(nrow(TP_muse))+as.numeric(nrow(FP_muse)))*(as.numeric(nrow(TP_muse))+as.numeric(nrow(FN_muse)))*(as.numeric(nrow(TN_muse))+as.numeric(nrow(FP_muse)))*(as.numeric(nrow(TN_muse))+as.numeric(nrow(FN_muse))))

TP_strelkaa <- platinum_AML31[platinum_AML31$GOLD %in% muts_strelka$mut,] 
TP_strelka <- muts_strelka[muts_strelka$GOLD==TRUE,] 
FN_strelka <- anti_join(platinum_AML31, TP_strelkaa) 
FP_strelka <- muts_strelka[muts_strelka$GOLD==FALSE,]
# TN_strelka <- unvalidated[!(unvalidated$mut %in% muts_strelka$mut),]
TN_strelka <- (table(factor(c(1:1000)))) #ESTO NO SE QUE VALOR METERLE PQ DEBERIA SER TODO EL GENOMA
Recall_strelka <- nrow(TP_strelka)/nrow(platinum_AML31) *100
#If TN from the normal-normal datasets are not included, specificity values will heavily depend on FP rates
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

#    row_names Recall_values Precision_values Specificity_values FP_values TP_values FN_values TN_values Accuracy_values  MCC_values F1_Score_values
# 1: 1 or more      4.839911         11.11111           65.78947       520        65      1278      1000        37.19874 -0.36351142        6.742739
# 2: 2 or more      4.616530         52.10084           94.60738        57        62      1281      1000        44.25000 -0.01774809        8.481532
# 3: 3 or more      4.393150         81.94444           98.71668        13        59      1284      1000        44.94907  0.08944694        8.339223
# 4: 4 or more      3.797468         96.22642           99.80040         2        51      1292      1000        44.81876  0.11974897        7.306590
# 5:         5      2.606106         94.59459           99.80040         2        35      1308      1000        44.13646  0.09553002        5.072464
# 6:   gatk4m2      3.871929         37.14286           91.91176        88        52      1291      1000        43.27437 -0.08999287        7.012812
# 7:    lancet      3.350707         30.00000           90.49774       105        45      1298      1000        42.68791 -0.12763901        6.028131
# 8:   caveman      4.169769         35.00000           90.57971       104        56      1287      1000        43.15488 -0.10568950        7.451763
# 9:      muse      4.244229         47.50000           94.07338        63        57      1286      1000        43.93184 -0.03837987        7.792208
# 10:  strelka      4.616530         20.94595           81.03728       234        62      1281      1000        41.21071 -0.22476245        7.565589

selectcols_snv_only_VC_vars <- c("T_NPROGS", "T_PROGS", "N_DP_avg", "N_DP_std", "N_DP_REF_avg", 
                                 "N_DP_ALT_avg", "N_VAF_avg", "N_DP_REF_cav", "N_DP_ALT_cav", "N_DP_REF_g4m2", "N_DP_ALT_g4m2", 
                                 "N_MBQ_g4m2", "N_MMQ_g4m2", "N_DP_REF_str", "N_DP_ALT_str", "N_DP_REF_mus", "N_DP_ALT_mus", 
                                 "N_BQ_REF_mus", "N_BQ_ALT_mus", "N_DP_REF_lan", "N_DP_ALT_lan", "T_NPROGS", "T_PROGS", "T_DP_avg", 
                                 "T_DP_std", "T_DP_REF_avg", "T_DP_ALT_avg", "T_VAF_avg", "T_DP_REF_cav", "T_DP_ALT_cav", "T_VAF_cav", 
                                 "T_ASMD_cav", "T_CLPM_cav", "T_DP_REF_g4m2", "T_DP_ALT_g4m2", "T_VAF_g4m2", "T_MBQ_g4m2", "T_MMQ_g4m2", 
                                 "T_DP_REF_str", "T_DP_ALT_str", "T_VAF_str", "T_MQ_str", "T_MQ0_str", "T_SomaticEVS_str", "T_DP_REF_mus", 
                                 "T_DP_ALT_mus", "T_VAF_mus", "T_BQ_REF_mus", "T_BQ_ALT_mus", "T_DP_REF_lan", "T_DP_ALT_lan", "T_VAF_lan", 
                                 "T_FETS_lan", "T_SB_lan")


muts2_only_VC_vars <- as.data.frame(muts2[,names(muts2) %in% selectcols_snv_only_VC_vars])
muts2_only_gold <- as.factor(muts2$GOLD)

summary(muts2_only_VC_vars)


# rf_snvs <- get(load("rf_snvs.RData"))
rf_snvs <- get(load("/home/computer/rf_snvs_noFPfromTCGA.RData"))

levels(muts2_only_VC_vars$T_PROGS)
# levels_training <- get(load("levels_training.RData"))
levels_training <- get(load("/home/computer/levels_training_noFPfromTCGA.RData"))

levels(muts2_only_VC_vars$T_PROGS) <- levels_training #Esto lo ha resuelto, si fuese alreves, sería necesario pegar estas rows en los datos de training y volverlas a quitar antes de hacer training para que conserve el numero de levels o hacer esta asignacion pero en el sentido opuesto



prediction_rf_snvs_AML31 <- predict(rf_snvs, muts2_only_VC_vars, type = "response")
tab <- table(prediction_rf_snvs_AML31, muts2_only_gold)
tab[1] <- as.numeric(tab[1])
tab[2] <- as.numeric(tab[2])
tab[3] <- as.numeric(tab[3])
tab[4] <- as.numeric(tab[4])






cf <- confusionMatrix(tab, positive = "TRUE")
cf

mcc <- (tab[4]*tab[1]-tab[2]*tab[3])  / (sqrt((tab[4]+tab[2])*(tab[4]+tab[3])*(tab[1]+tab[2])*(tab[1]+tab[3])))

f1_score <- 2 * tab[4] / (2 * tab[4]+ tab[2]+tab[3])



draw_confusion_matrix <- function(cm) {
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('Random Forest Confusion Matrix Results for SNVs', cex.main = 2)
  
  rect(150, 430, 240, 370, col = '#3F97D0')
  text(195, 440, 'No variant', cex = 1.5)
  rect(250, 430, 340, 370, col = '#F7AD50')
  text(295, 440, 'Variant', cex = 1.5)
  text(125, 370, 'Predicted', cex = 1.3, srt = 90, font = 2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col = '#F7AD50')
  rect(250, 305, 340, 365, col = '#3F97D0')
  text(140, 400, 'No variant', cex = 1.5, srt = 90)
  text(140, 335, 'Variant', cex = 1.5, srt = 90)
  
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex = 1.6, font = 2, col = 'white')
  text(195, 385, "TN", cex = 1.6, font = 2, col = 'white')
  text(195, 335, res[2], cex = 1.6, font = 2, col = 'white')
  text(195, 320, "FP", cex = 1.6, font = 2, col = 'white')
  text(295, 400, res[3], cex = 1.6, font = 2, col = 'white')
  text(295, 385, "FN", cex = 1.6, font = 2, col = 'white')
  text(295, 335, res[4], cex = 1.6, font = 2, col = 'white')
  text(295, 320, "TP", cex = 1.6, font = 2, col = 'white')
  
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt = 'n', yaxt = 'n')
  text(10, 85, "Absolute Recall", cex=1.4, font=2)
  text(10, 68, (round(as.numeric(cm$byClass[6])*0.01*Recall_1, 3)), cex=1.4)
  #text(10, 68, (round(as.numeric(cm$byClass[6]), 3) * 0.920), cex=1.4)
  text(30, 85, names(cm$byClass[2]), cex=1.4, font=2)
  text(30, 68, round(as.numeric(cm$byClass[2]), 3), cex=1.4)
  text(50, 85, names(cm$byClass[5]), cex=1.4, font=2)
  text(50, 68, round(as.numeric(cm$byClass[5]), 3), cex=1.4)
  text(70, 85, names(cm$byClass[6]), cex=1.4, font=2)
  text(70, 68, round(as.numeric(cm$byClass[6]), 3), cex=1.4)
  text(90, 85, names(cm$byClass[7]), cex=1.4, font=2)
  text(90, 68, round(as.numeric(cm$byClass[7]), 3), cex=1.4)
  # text(10, 35, "F1-Score", cex=1.4, font=2)
  # text(10, 18, round(f1_score,3), cex=1.4)
  text(30, 35, names(cm$overall[1]), cex=1.4, font=2)
  text(30, 18, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(50, 35, "MCC", cex=1.4, font=2)
  text(50, 18, round(mcc,3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.4, font=2)
  text(70, 18, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  

draw_confusion_matrix(cf)







votes=predict(rf_snvs, muts2_only_VC_vars, type = "vote")

votes <- as.data.frame(votes)
votes$GOLD <- muts2_only_gold
colnames(votes) <- c("NO","YES","GOLD")
rownames(votes) <- c()



#FIGURA VOTOS
ggplot(votes) +
  geom_point(aes(x=NO,y=YES, color=GOLD),alpha=0.4) +ylim(0,1)+xlim(0,1) +
  ggtitle("Proportion of trees that voted the variant as:")+ xlab("Not a variant") + ylab("Variant") + theme_bw() + 
  theme(legend.text=element_text(size=12),axis.text=element_text(size=15),axis.title.y=element_text(size=15),axis.title.x=element_text(size=15),title=element_text(size=15)) +
  scale_color_viridis(option="viridis", begin = 0.3,name = "Benchmarking status",discrete=TRUE, direction = 1)

#CON JITTER
ggplot(votes, aes(x=NO,y=YES, color=GOLD, alpha=0.2)) +
  geom_jitter(width=0.03, height =0.03) +ylim(-0.1,1.05)+xlim(-0.1,1.05) +
  ggtitle("Proportion of trees that voted the variant as:")+ xlab("Not a variant") + ylab("Variant") + theme_bw() + 
  theme(legend.text=element_text(size=12),axis.text=element_text(size=15),axis.title.y=element_text(size=15),axis.title.x=element_text(size=15),title=element_text(size=15)) +
  scale_color_viridis(option="viridis", begin = 0.3,name = "Benchmarking status",discrete=TRUE, direction = 1)



