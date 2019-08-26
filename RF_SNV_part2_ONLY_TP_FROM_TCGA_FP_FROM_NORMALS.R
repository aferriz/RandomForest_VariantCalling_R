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


#Añadir información del estatus de la vairante TP o FP a dataframe con info de las variables
#Ver a ver si lo puedo hacer de otra manera

#Añadir info status, y separar en TRUE y FALSE. Luego de esos coger el 70% de cada uno y crear el tdb y el 30% restante con lo del - y crear el vDB y ya esta

#Comienzo por el que esta igualado


#LO PRIMERO DE TODO AHORA ES COGER EL DE TP_1 Y FP_1 Y SELECCIONAR EL TRAINING Y VALIDATION CON NUMERO QUE PUEDA SER PERMUTADO
#UNA VEZ SEPARADOS SERA DOBLE PR TP Y PARA FP  FUCK TOCA REPETIR ALGUN PASO SI O SI, ASIQUE NO IMPORTA CUAL BAY PRIMERO
#CREAR COLUMAS NUEVAS PARA EL DE TP_1 FP_1 DE GOLD TODA TRUE Y TODA FALSE Y LUEGO HACER RBIND

TP_model <- TP_1
FP_model <- FP_1 

colnames(TP_model)


TP_model[is.na(TP_model)] <- -1
FP_model[is.na(FP_model)] <- -1

TP_model$GOLD <- TRUE
FP_model$GOLD <- FALSE


summary(TP_model)
summary(FP_model)

set.seed(500)
vector_sample_TP <- sample(1:nrow(TP_model), nrow(TP_model)*0.7,replace = F)
training_TP <- TP_model[vector_sample_TP,]
validation_TP <-TP_model[-vector_sample_TP,]

vector_sample_FP <- sample(1:nrow(FP_model), nrow(FP_model)*0.7,replace = F)
training_FP <- FP_model[vector_sample_FP,]
validation_FP <-FP_model[-vector_sample_FP,]

training_snvs <- rbind(training_TP,training_FP)
validation_snvs <- rbind(validation_TP,validation_FP)

selectcols_snv_only_VC_vars <- c("T_NPROGS", "T_PROGS", "N_DP_avg", "N_DP_std", "N_DP_REF_avg", 
                                 "N_DP_ALT_avg", "N_VAF_avg", "N_DP_REF_cav", "N_DP_ALT_cav", "N_DP_REF_g4m2", "N_DP_ALT_g4m2", 
                                 "N_MBQ_g4m2", "N_MMQ_g4m2", "N_DP_REF_str", "N_DP_ALT_str", "N_DP_REF_mus", "N_DP_ALT_mus", 
                                 "N_BQ_REF_mus", "N_BQ_ALT_mus", "N_DP_REF_lan", "N_DP_ALT_lan", "T_NPROGS", "T_PROGS", "T_DP_avg", 
                                 "T_DP_std", "T_DP_REF_avg", "T_DP_ALT_avg", "T_VAF_avg", "T_DP_REF_cav", "T_DP_ALT_cav", "T_VAF_cav", 
                                 "T_ASMD_cav", "T_CLPM_cav", "T_DP_REF_g4m2", "T_DP_ALT_g4m2", "T_VAF_g4m2", "T_MBQ_g4m2", "T_MMQ_g4m2", 
                                 "T_DP_REF_str", "T_DP_ALT_str", "T_VAF_str", "T_MQ_str", "T_MQ0_str", "T_SomaticEVS_str", "T_DP_REF_mus", 
                                 "T_DP_ALT_mus", "T_VAF_mus", "T_BQ_REF_mus", "T_BQ_ALT_mus", "T_DP_REF_lan", "T_DP_ALT_lan", "T_VAF_lan", 
                                 "T_FETS_lan", "T_SB_lan")

training_snvs_only_VC_vars <- as.data.frame(training_snvs[,names(training_snvs) %in% selectcols_snv_only_VC_vars])
training_snvs_only_gold <- as.factor(as.character(training_snvs$GOLD))
# training_snvs_only_VC_vars$T_PROGS <- as.factor(as.character(training_snvs_only_VC_vars$T_PROGS))


summary(training_snvs_only_gold)
summary(training_snvs_only_VC_vars)


nrow(training_snvs_only_VC_vars)
nrow(training_snvs_only_gold)



validation_snvs_only_VC_vars <- as.data.frame(validation_snvs[,names(validation_snvs) %in% selectcols_snv_only_VC_vars])
validation_snvs_only_gold <- as.factor(validation_snvs$GOLD)


#La response debe ser un vector no un dataframe ni lista. 
rf_snvs <-randomForest(training_snvs_only_VC_vars, training_snvs_only_gold, ntree=100, mtry=30 ,maxnodes = 200 )

table(training_snvs_only_gold)
mode(training_snvs_only_VC_vars)
mode(training_snvs_only_gold)

levels_training <- levels(training_snvs_only_VC_vars$T_PROGS)

save(rf_snvs, file = "/home/computer/rf_snvs_noFPfromTCGA.RData")
save(levels_training, file = "/home/computer/levels_training_noFPfromTCGA.RData") #para evitar error por diferencia de niveles



prediction_rf_snvs <- predict(rf_snvs, validation_snvs_only_VC_vars, type = "response")
tab <- table(prediction_rf_snvs, validation_snvs_only_gold)
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

# > cf
# Confusion Matrix and Statistics
# 
# validation_snvs_only_gold
# prediction_rf_snvs FALSE TRUE
# FALSE  2276   46
# TRUE     64 2482
# 
# Accuracy : 0.9774          
# 95% CI : (0.9728, 0.9814)
# No Information Rate : 0.5193          
# P-Value [Acc > NIR] : <2e-16          
# 
# Kappa : 0.9547          
# 
# Mcnemar's Test P-Value : 0.105           
# 
# Sensitivity : 0.9818          
# Specificity : 0.9726          
# Pos Pred Value : 0.9749          
# Neg Pred Value : 0.9802          
# Prevalence : 0.5193          
# Detection Rate : 0.5099          
# Detection Prevalence : 0.5230          
# Balanced Accuracy : 0.9772          
# 
# 'Positive' Class : TRUE     

#Ahora que funciona bien para la validation del propio tipo de datos TCGA y FP de normal bams
#Es cuando hay que probarlo sobre AML31 repitiendo todo pero respecto al RF no hacer training, 
#solo la parte de predict y ver como sale si sale mejor que los programas por separado


votes=predict(rf_snvs, validation_snvs_only_VC_vars, type = "vote")

votes <- as.data.frame(votes)
votes$GOLD <- validation_snvs_only_gold
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



imp <- varImpPlot(rf_snvs)
# this part just creates the data.frame for the plot part
library(dplyr)
imp <- as.data.frame(imp)
imp$varnames <- rownames(imp) # row names to column
rownames(imp) <- NULL  
imp$var_categ <- rep(1:2, 15) # random var category

imp$var_categ <- ifelse(imp$varnames %like% "_pin" | imp$varnames %like% "_lan" | imp$varnames %like% "_str"
                        | imp$varnames %like% "_g4m2" | imp$varnames %like% "_pla" | imp$varnames %like% "_sva","PROG_SPECIFIC","COMMON")



ggplot(imp, aes(x=reorder(varnames, MeanDecreaseGini), y=MeanDecreaseGini, color=as.factor(var_categ))) + 
  geom_point()   +
  geom_segment(aes(x=varnames,xend=varnames,y=0,yend=MeanDecreaseGini)) +
  ylab("MeanDecreaseGini") +
  xlab("Variable Name") +
  coord_flip()+ theme_bw() +scale_color_viridis(option="viridis", begin = 0.1, end= 1,name = "Variable Group",discrete=TRUE, direction = -1)














