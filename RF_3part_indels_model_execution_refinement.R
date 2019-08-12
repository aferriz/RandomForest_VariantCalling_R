# RANDOM FOREST EXECUTION AND REFINEMENT


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

save(tDB_only_VC_variables, file = "tDB_only_VC_variables.RData")
#Pasar al otro script
save(tDB_only_VC_variables_and_golden,file="tDB_only_VC_variables_and_golden.RData")


tDB_only_VC_variables <- get(load("nuevo_tDB.RData"))

set.seed(500)
# rf <- randomForest(x[,seq(5,ncol(x)-2,1)], factor(y), ntree=500)
# prediction <- predict(rf, x_val[,seq(5,ncol(x_val)-1,1)], type = "response")
# table(prediction,truth)

#k = valor de cutoff EN el gráfico de votos
# k=0.25
# rf <- randomForest(tDB_only_VC_variables, tDB_golden_column, ntree=100,cutoff=c(k,1-k), mtry=30 ,maxnodes = 10000 )


nmin <-  sum(tDB_golden_column == "FALSE")
rf <- randomForest(tDB_only_VC_variables, tDB_golden_column, strata=tDB_golden_column, sampsize=rep(nmin,2),ntree=100, mtry=30 ,maxnodes = 300 )
summary(tDB_only_VC_variables)
colnames(tDB_only_VC_variables)

prediction_training <- predict(rf, tDB_only_VC_variables, type = "response")
table(prediction_training, tDB_golden_column)
save(rf, file = "RF_10julio.RData")




############################
############################

rf <- randomForest(tDB_only_VC_variables, tDB_golden_column, ntree=100, mtry=30 ,maxnodes = 200 )
summary(tDB_only_VC_variables)
colnames(tDB_only_VC_variables)

prediction_training <- predict(rf, tDB_only_VC_variables, type = "response")
table(prediction_training, tDB_golden_column)
save(rf, file = "RF_10julio.RData")


sapply(tDB_only_VC_variables, class)
print(rf)
prediction <- predict(rf, vDB_only_VC_variables, type = "response")
table(prediction, vDB_golden_colum)

#              truth
# prediction FALSE TRUE
# FALSE          0    0
# TRUE          44 4795



#               tDB_golden_column
# prediction_training FALSE TRUE
# FALSE                  60    0
# TRUE                    6 7191

# object.size(rf)
# View(factor(y))
# algo <- x[,seq(5,ncol(x)-1,1)]
# summary(algo)
#factor(y) es exactamente igual a tDB_GOLD_colum
#Por tanto la diferencia en el resultado solo viene de x[,seq(5,ncol(x)-1,1)]
#Lo mismo podria eliminarlo en el paso de tDB esas rows y todo iria bien. ESTO ES LO QUE HE HECHO

#localIMP Should casewise importance measure be computed? (Setting this to TRUE will override importance.

# rf2 <- randomForest(GOLD2 ~ N_DP_avg + N_DP_std + T_NPROGS + T_PROGS + T_DP_avg +  T_DP_std + T_VAF_avg +  
                      # T_MBQ_g4m2 + T_MMQ_g4m2 + T_MQ_str + T_MQ0_str + T_SomaticEVS_str + T_RC_str + T_IC_str + 
                      # T_IHP_str + T_FETS_lan + T_SB_lan + T_QUAL_pin + T_S1_pin + T_S2_pin + T_GQ_pla + T_GOF_pla + 
                      # T_MMLQ_pla + T_HP_pla + T_MQ_pla + T_QD_pla + T_LO_sva + T_LR_sva + T_MAPQ_sva + 
                      # T_LOD_sva, tDB, ntree=500, importance =TRUE)


################# este paquete no esta disponible para estas versiones de R
#Cargar este paquete. Me toca instalar una version de R nueva, y que pasa si los paquetes no estan actualizados para funcionar con 3.6?
# Despues de actualizar a 3.6.0. sigue sin funcionar, porqueria de paquete no lo uso, usamos otro,
#install.packages("reprtree")


#Despues de tener tanto el x, como el tDB con lo de as.factor en el T_PROGS ahora da error de NA en esto:
# reprtree:::plot.getTree(rf2,1)
# reprtree:::plot.getTree(rf2,k=300)
# reprtree:::plot.getTree(rf2,500)
# Warning messages:
#   1: In intToBits(x) : NAs introduced by coercion to integer range
#Esto pasa pq tiene que caluclar un montón de números y el número producido es muy grande, más de lo que r CONTEMPLA
#Para solucionarlo no se si se podría cargando el paquete MPFR
# https://cran.r-project.org/web/packages/Rmpfr/index.html
# dentro de la función de plot.getTree
#https://github.com/araastat/reprtree/blob/master/R/plot.getTree.R
##############
#OJITO CON ESTO, QUE ES SOLO LA GRÁFICA, ASIQUE EN PRINCIPIO NO DEBE AFECTAR A LOS CALCULOS PERO ESO...


# library(rpart)
# library(rpart.plot)
# 
# # https://stats.stackexchange.com/questions/41443/how-to-actually-plot-a-sample-tree-from-randomforestgettree
# tree <- getTree(rf, k=1, labelVar=TRUE)
# tree
# realtree <- reprtree:::as.tree(tree, rf)
# #realtree
# 
# print(rf)
# print(rf2)
# 
# library(rattle)

# !!!
# FYI: reprtree:::plot.getTree() seems not to wor k if the random forest is generated using the alternative way 
# e.g. randomForest(iris[,-5], iris[,5])

rf$confusion
head(rf$votes)
head(rf$predicted)
head(rf$importance)

varImpPlot(rf,main="Important variables")
###
# Explain random forest
# https://rawgit.com/geneticsMiNIng/BlackBoxOpener/master/randomForestExplainer/inst/doc/randomForestExplainer.html
###
# importance_frame <- measure_importance(rf)
# plot_multi_way_importance(importance_frame, size_measure = "no_of_nodes", min_no_of_trees = 30)


# "Warning: your forest does not contain information on local importance so 'accuracy_decrease' 
# measure cannot be extracted. To add it regrow the forest with the option localImp = TRUE and 
# run this function again."

#YA ESTA EL factor(y) lo que hace es convertir en factor a la y, la y es la columna del tDB$GOLD. No es nada del validation ni test

# set.seed(500)
# rf <- randomForest(x[,seq(5,ncol(x)-2,1)], factor(y), ntree=500, localImp=TRUE)
# object.size(rf)

importance_frame <- measure_importance(rf)

plot_multi_way_importance(importance_frame, size_measure = "no_of_nodes", min_no_of_trees = 30)
plot_multi_way_importance(importance_frame, x_measure = "accuracy_decrease", y_measure = "gini_decrease", size_measure = "p_value")
plot_importance_ggpairs(importance_frame)

vars <- important_variables(importance_frame, k = 20, measures = c("mean_min_depth", "no_of_trees"))

#This step takes a lot of time, run in terminal with R
#hasta aquí se puede ejecutar en Rstudio, para ejecutar este paso hay que borrar los objetos más pesados desde la vision en Grid no en List
interactions_frame <- min_depth_interactions(rf, vars)
# Warning messages:
# 1: Factor `split var` contains implicit NA, consider using `forcats::fct_explicit_na` 
# 2: Factor `split var` contains implicit NA, consider using `forcats::fct_explicit_na` 

head(interactions_frame[order(interactions_frame$occurrences, decreasing = TRUE), ])
#Resultado hecho en R terminal!!!
# variable root_variable mean_min_depth occurrences         interaction
# 540 T_VAF_avg     T_VAF_avg       2.501305         383 T_VAF_avg:T_VAF_avg
# 20   N_DP_avg     T_VAF_avg       3.726930         370  T_VAF_avg:N_DP_avg
# 60   T_DP_avg     T_VAF_avg       3.207029         369  T_VAF_avg:T_DP_avg
# 40   N_DP_std     T_VAF_avg       4.370052         348  T_VAF_avg:N_DP_std
# 43   T_DP_avg      T_DP_avg       3.621760         326   T_DP_avg:T_DP_avg
# 180 T_IHP_str     T_VAF_avg       4.175055         326 T_VAF_avg:T_IHP_str
# uncond_mean_min_depth
# 540              2.536000
# 20               4.250000
# 60               3.184000
# 40               4.126944
# 43               3.184000
# 180              3.315888


#Resultado hecho en Rstudio
# variable root_variable mean_min_depth occurrences         interaction uncond_mean_min_depth
# 567 T_VAF_avg     T_VAF_avg       2.417989         378 T_VAF_avg:T_VAF_avg                 2.570
# 63   T_DP_avg     T_VAF_avg       2.621270         362  T_VAF_avg:T_DP_avg                 2.858
# 21   N_DP_avg     T_VAF_avg       3.503466         359  T_VAF_avg:N_DP_avg                 3.814
# 189 T_IHP_str     T_VAF_avg       3.469021         357 T_VAF_avg:T_IHP_str                 3.376
# 45   T_DP_avg      T_DP_avg       3.077190         345   T_DP_avg:T_DP_avg                 2.858
# 42   N_DP_std     T_VAF_avg       4.090291         341  T_VAF_avg:N_DP_std                 3.656

#Resultado hecho en Rstudio tras los nuevos cambios
# variable root_variable mean_min_depth occurrences         interaction uncond_mean_min_depth
# 13   N_DP_avg       T_PROGS       2.694517         383    T_PROGS:N_DP_avg                 4.168
# 53   T_DP_avg       T_PROGS       2.305698         376    T_PROGS:T_DP_avg                 2.854
# 593 T_VAF_avg       T_PROGS       2.481501         372   T_PROGS:T_VAF_avg                 2.948
# 600 T_VAF_avg     T_VAF_avg       3.189796         344 T_VAF_avg:T_VAF_avg                 2.948
# 60   T_DP_avg     T_VAF_avg       3.082360         335  T_VAF_avg:T_DP_avg                 2.854
# 20   N_DP_avg     T_VAF_avg       3.854454         331  T_VAF_avg:N_DP_avg                 4.168
# 


plot_min_depth_interactions(interactions_frame) # Esta es la que he generado con R en terminal

###
# https://stackoverflow.com/questions/52200095/how-to-customize-the-importance-plot-generated-by-package-randomforest
###

#rf <- randomForest(x[,seq(5,ncol(x)-1,1)], factor(y), ntree=1000, keep.forest=FALSE,
#                          importance=TRUE)

imp <- varImpPlot(rf)
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

###
# https://cran.rstudio.com/web/packages/randomForestExplainer/vignettes/randomForestExplainer.html
###
# (vars <- important_variables(forest, k = 5, measures = c("mean_min_depth", "no_of_trees"))) # gives the same result as below but takes longer
# library(randomForestExplainer)
# (vars <- important_variables(importance_frame, k = 5, measures = c("mean_min_depth", "no_of_trees")))
# interactions_frame <- min_depth_interactions(forest, vars)

# IncNodePurity ???
# ggplot(imp, aes(x=reorder(varnames, IncNodePurity), y=IncNodePurity, color=as.factor(var_categ))) + 
#   geom_point() +
#   geom_segment(aes(x=varnames,xend=varnames,y=0,yend=IncNodePurity)) +
#   scale_color_discrete(name="Variable Group") +
#   ylab("IncNodePurity") +
#   xlab("Variable Name") +
#   coord_flip()

###


## Create validation series:
# x_val <- vDB[,1:(ncol(vDB)-1)]
# summary(x_val)
# # x[,8] <- as.factor(x[,8])
# 
# #ESTAS SON LAS VARIANTES QUE ESTAN EN EL SUPUESTO GOLDEN DE TARGETED, INCLUYEN P Y N. Como esto es el validation solo hay el 40%. De un total de 111 pues aqui hay 
# truth <- vDB[, ncol(vDB)] #Esto relamente lo que hace es coger la ultima columa que es la de $GOLD y es la que tiene la info de si la variante es TRUE or FALSE
# rrf <- x[,seq(5,ncol(x)-1,1)]
# #OJOOOOOOOOOOOOOO AQUI LO DE TRUTH REALMENTE QUIERE DECIR TEST 40%
# 
# 
# View(truth)

#mirar esta solucion para lo de: Error in predict.randomForest(rf, x[, seq(5, ncol(x) - 1, 1)], type = "response") : New factor levels not present in the training data
#https://stats.stackexchange.com/questions/235764/new-factors-levels-not-present-in-training-data



#Guardar el objeto RF??
# save(rf, file = "RandomForest_60.RData")
# print(rf)
# 
# 
# NUEVO_rf <- get(load("/home/computer/RandomForest_60.RData"))
# print(NUEVO_rf)
# 
# NUEVO_prediction <- predict(NUEVO_rf, x[,seq(5,ncol(x)-1,1)], type = "response")
# table(NUEVO_prediction, truth)
# NUEVO_tab <- table(NUEVO_prediction, truth)
# NUEVO_tab[3] <- NUEVO_tab[3] + nrow(vDB_FN) 
# NUEVO_tab


# prediction <- predict(rf, x_val[,seq(5,ncol(x_val)-1,1)], type = "response")
# prediction_training <- predict(rf, x[,seq(5,ncol(x)-1,1)], type = "response")
# 
# length(x_val)
# training <- tDB[, ncol(tDB)]
# 
# table(prediction_training,training)
# 
# table(prediction, truth)
# 
# length(prediction)
# length(truth)

#confusionMatrix(prediction, truth)
tab <- table(prediction, vDB_golden_colum)
tab[3] <- tab[3] + nrow(vDB_FN) #Esto es totalmente artificial y debo mejorarlo y entenderlo
tab
cf <- confusionMatrix(tab, positive = "TRUE")
cf

draw_confusion_matrix <- function(cm) {
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('Confusion Matrix', cex.main = 2)
  
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
  text(10, 68, (round(as.numeric(cm$byClass[6]), 3) * 0.920), cex=1.4)
  text(30, 85, names(cm$byClass[2]), cex=1.4, font=2)
  text(30, 68, round(as.numeric(cm$byClass[2]), 3), cex=1.4)
  text(50, 85, names(cm$byClass[5]), cex=1.4, font=2)
  text(50, 68, round(as.numeric(cm$byClass[5]), 3), cex=1.4)
  text(70, 85, names(cm$byClass[6]), cex=1.4, font=2)
  text(70, 68, round(as.numeric(cm$byClass[6]), 3), cex=1.4)
  text(90, 85, names(cm$byClass[7]), cex=1.4, font=2)
  text(90, 68, round(as.numeric(cm$byClass[7]), 3), cex=1.4)
  text(30, 35, names(cm$overall[1]), cex=1.4, font=2)
  text(30, 18, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.4, font=2)
  text(70, 18, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  
draw_confusion_matrix(cf)




#GRÁFICAS DE LAS VARIABLES
#Se puede hacer con for, con el tDB_only_VC_variables_and_golden
#Y que haces con los -1?
a__1 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(N_DP_avg),fill='#E69F00')+ theme(legend.position = "none")+ggtitle("TP training")


a__2 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(N_DP_std),fill='#E69F00')+
  labs(fill='Variant Detected')  + theme(legend.position = "none")+ggtitle("TP training")


a__3 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_DP_avg),binwidth=5)+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__4 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_DP_std),binwidth=5)+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")


a__5 <-ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_PROGS),stat="count")+ 
  labs(fill='Variant Detected')+theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1)) + theme(legend.position = "none")+ggtitle("TP training")


a__6 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_NPROGS),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__7 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_VAF_avg),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__8 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_MBQ_g4m2),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__9 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_MMQ_g4m2),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__10 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_MQ_str),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__11 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_MQ0_str),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__12 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_SomaticEVS_str),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__13 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_RC_str),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__14 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_IC_str),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__15 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_IHP_str),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__16 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_FETS_lan),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__17 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_SB_lan),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__18 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_QUAL_pin),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__19 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_S1_pin),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__20 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_S2_pin),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__21 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_GQ_pla),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__22 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_GOF_pla),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__23 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_MMLQ_pla),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__24 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_HP_pla),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__25 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_MQ_pla),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__26 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_QD_pla),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__27 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_LO_sva),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__28 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_LR_sva),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__29 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_MAPQ_sva),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")

a__30 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==TRUE,])+
  geom_histogram(aes(T_LOD_sva),fill='#E69F00')+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("TP training")


####################33 #ESTE COPIO Y PEGO 4 VECES CAMBIANDO LO NECESARIO, b__, COLOR, false Y data.
b__1 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(N_DP_avg))+
  labs(fill='Variant Detected')  + theme(legend.position = "none")+ggtitle("FP training")


b__2 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(N_DP_std))+
  labs(fill='Variant Detected')  + theme(legend.position = "none")+ggtitle("FP training")


b__3 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_DP_avg),binwidth=5)+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__4 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_DP_std),binwidth=5)+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")


b__5 <-ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_PROGS),stat="count")+ 
  labs(fill='Variant Detected')+theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1)) + theme(legend.position = "none")+ggtitle("FP training")


b__6 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_NPROGS))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__7 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_VAF_avg))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__8 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_MBQ_g4m2))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__9 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_MMQ_g4m2))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__10 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_MQ_str))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__11 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_MQ0_str))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__12 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_SomaticEVS_str))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__13 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_RC_str))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__14 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_IC_str))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__15 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_IHP_str))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__16 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_FETS_lan))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__17 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_SB_lan))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__18 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_QUAL_pin))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__19 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_S1_pin))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__20 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_S2_pin))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__21 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_GQ_pla))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__22 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_GOF_pla))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__23 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_MMLQ_pla))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__24 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_HP_pla))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__25 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_MQ_pla))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__26 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_QD_pla))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__27 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_LO_sva))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__28 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_LR_sva))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__29 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_MAPQ_sva))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")

b__30 <- ggplot(tDB_only_VC_variables_and_golden[tDB_only_VC_variables_and_golden$GOLD2==FALSE,])+
  geom_histogram(aes(T_LOD_sva))+
  labs(fill='Variant Detected') + theme(legend.position = "none")+ggtitle("FP training")



#####################

library(gridExtra)
grid.arrange(a__1, a__2, b__1,b__2,nrow = 2, ncol=2)













colnames(tDB_only_VC_variables_and_golden)

# # OJO HI HA COMBINACIONS PROGS QUE NO ESTAVEN AL MB99 !! lancet,muse,strelka
# ## Creating sample data
# values_development=tDB$T_PROGS ## Values used when building the random forest model
# # values_production=vDB_AML31$T_PROGS ## New values to used when using the model
# ## Deleting cases which were not present when developing
# values_production=sapply(as.character(values_production), function(x) if(x %in% values_development) x else NA)
# ## Creating the factor variable, (with the correct NA value level)
# values_production=factor(values_production)
# ## Checking
# values_production # =>  a     b     c  <NA> 
# unique(values_production)
# 
# x2$T_PROGS[!x2$T_PROGS %in% tDB$T_PROGS]
# ### ARREGLO per si hi ha combiniacions que no estaven al training (NO ES EL CAS !!!)
# x2$T_PROGS <- sapply(as.character(values_production), function(x) if(x %in% values_development) x else NA)
# x2$T_PROGS <- as.factor(x2$T_PROGS)
# nrow(x2)
# x2 <- x2[!is.na(x2$T_PROGS),]
# nrow(x2)
# ### ARREGLO PEL QUE FALTEN
# levels(tDB$T_PROGS)
# levels(x2$T_PROGS)
# levels(x2$T_PROGS) <- levels(tDB$T_PROGS)
# levels(x2$T_PROGS)
# ###
# 
# prediction2 <- predict(rf, x2[,seq(5,ncol(x)-1,1)], type = "response")
# table(prediction2, truth2)
# 
# tab2 <- table(prediction2, truth2)
# #tab[3] <- tab[3] + nrow(vDB_FN)
# cf <- confusionMatrix(tab2, positive = "TRUE")
# cf
# 
# draw_confusion_matrix(cf)
# 
# nrow(x2)
# 
# nrow(vDB_FN)
# 
# head(truth)
# head(prediction)
# length(prediction)

###########
# get FN: FN from vDB & FN not in list
###########

# # FN from vDB
# head(vDB)
# head(prediction)
# nrow(vDB)
# length(prediction)
# prediction_df <- as.data.frame(prediction)
# nrow(prediction_df)
# head(prediction_df)
# vDB2 <- vDB
# vDB2$prediction <- prediction_df$prediction
# head(vDB2)
# FN_from_vDB2 <- vDB2[vDB2$prediction=="FALSE" & vDB2$GOLD,]
# nrow(FN_from_vDB2)
# TP_from_vDB2 <- vDB2[vDB2$prediction=="TRUE" & vDB2$GOLD,]
# nrow(TP_from_vDB2)
# 
# ggplot(FN_from_vDB2)+
#   geom_bar(aes(T_PROGS,fill=prediction))+
#   theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))+
#   labs(fill='Variant Detected') +ggtitle("Combinations of programs: FALSE NEGATIVES")

# ggplot(vDB2)+
#   geom_bar(aes(T_PROGS,fill=GOLD))+
#   theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))+
#   labs(fill='Variant in Golden') +ggtitle("Combinations of programs: TRUE VARIANTS") 
# 
# ggplot(vDB2)+
#   geom_bar(aes(T_PROGS,fill=prediction))+
#   theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))+
#   labs(fill='Variant Detected') +ggtitle("Combinations of programs: PREDICTION") 
# 
# FN_from_vDB2 <- vDB2[vDB2$prediction=="FALSE" & vDB2$GOLD,]
# nrow(FN_from_vDB2)
# TP_from_vDB2 <- vDB2[vDB2$prediction=="TRUE" & vDB2$GOLD,]
# nrow(TP_from_vDB2)
# 
# test <- vDB2
# test$AA <- ifelse(vDB2$prediction=="FALSE" & vDB2$GOLD,"FN","UNK")
# table(test$AA)
# test$AA <- ifelse(vDB2$prediction=="TRUE" & vDB2$GOLD,"TP",test$AA)
# test$AA[test$AA=="UNK"] <- "TN"
# table(test$AA)
# 
# ggplot(test)+
#   geom_bar(aes(T_PROGS,fill=AA))+
#   theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))+
#   labs(fill='Variant Detected') +ggtitle("Combinations of programs: PREDICTION")
# 
# ggplot(test)+
#   geom_bar(aes(T_PROGS,fill=AA))+
#   theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))+
#   labs(fill='Variant Detected') +ggtitle("Combinations of programs: PREDICTION") +
#   scale_y_log10() + ylab("log10(count)")
# 
# ggplot(test)+
#   geom_bar(aes(T_DP_std,fill=AA))+
#   theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))+
#   labs(fill='Variant Detected') +ggtitle("Combinations of programs: PREDICTION") +
#   scale_y_log10() + ylab("log10(count)") +
#   xlim(-1,10)
# 
# ggplot(test)+
#   geom_density(aes(T_DP_std,fill=GOLD), alpha=0.5) +
#   xlim(-1,10)


# prob=predict(rf,x,type="prob")
# 
# plot(prob)
# nrow(prob)
# head(prob)

votes=predict(rf,vDB_only_VC_variables,type="vote")
votes <- as.data.frame(votes)
votes$GOLD <- vDB$GOLD
colnames(votes) <- c("NO","YES","GOLD")
rownames(votes) <- c()

votes2=predict(rf,tDB_only_VC_variables,type="vote")
votes2 <- as.data.frame(votes2)
votes2$GOLD <- tDB$GOLD
colnames(votes2) <- c("NO","YES","GOLD")
rownames(votes2) <- c()
# nrow(votes)
# head(votes)
# plot(votes)



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

#CON JITTER
ggplot(votes2, aes(x=NO,y=YES, color=GOLD, alpha=0.2)) +
  geom_jitter(width=0.03, height =0.03) +ylim(-0.1,1.05)+xlim(-0.1,1.05) +
  ggtitle("Proportion of trees that voted the variant as:")+ xlab("Not a variant") + ylab("Variant") + theme_bw() + 
  theme(legend.text=element_text(size=12),axis.text=element_text(size=15),axis.title.y=element_text(size=15),axis.title.x=element_text(size=15),title=element_text(size=15)) +
  scale_color_viridis(option="viridis", begin = 0.3,name = "Benchmarking status",discrete=TRUE, direction = 1)



# votes=predict(rf,x2,type="vote")
# votes <- as.data.frame(votes)
# nrow(votes)
# head(votes)
# plot(votes)

# votes$GOLD <- vDB_AML31$PLAT
# colnames(votes) <- c("NO","YES","GOLD")
# 
# ggplot(votes) +
#   geom_point(aes(x=NO,y=YES,color=GOLD))
# 
# plot(x=votes$"FALSE",y=votes$"TRUE")
# 
# 

###
# https://stackoverflow.com/questions/23891140/r-how-to-visualize-confusion-matrix-using-the-caret-package
###
# https://stats.stackexchange.com/questions/253351/predict-function-tuning-for-random-forest
# 
# draw_confusion_matrix <- function(cm) {
#   
#   layout(matrix(c(1,1,2)))
#   par(mar=c(2,2,2,2))
#   plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
#   title('CONFUSION MATRIX', cex.main=2)
#   
#   # create the matrix 
#   rect(150, 430, 240, 370, col='#3F97D0')
#   text(195, 435, 'FALSE', cex=1.2)
#   rect(250, 430, 340, 370, col='#F7AD50')
#   text(295, 435, 'TRUE', cex=1.2)
#   text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
#   text(245, 450, 'Actual', cex=1.3, font=2)
#   rect(150, 305, 240, 365, col='#F7AD50')
#   rect(250, 305, 340, 365, col='#3F97D0')
#   text(140, 400, 'FALSE', cex=1.2, srt=90)
#   text(140, 335, 'TRUE', cex=1.2, srt=90)
#   
#   # add in the cm results 
#   res <- as.numeric(cm$table)
#   text(195, 400, res[1], cex=1.6, font=2, col='white')
#   text(195, 335, res[2], cex=1.6, font=2, col='white')
#   text(295, 400, res[3], cex=1.6, font=2, col='white')
#   text(295, 335, res[4], cex=1.6, font=2, col='white')
#   
#   # add in the specifics 
#   plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
#   text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
#   text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
#   text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
#   text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
#   text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
#   text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
#   text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
#   text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
#   text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
#   text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
#   
#   # add in the accuracy information 
#   text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
#   text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
#   text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
#   text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
# }  



#Finally, pass in the cm object that we calculated when using caret to create the confusion matrix:

# construct the evaluation dataset
# set.seed(144)
# true_class <- factor(sample(paste0("Class", 1:2), size = 1000, prob = c(.2, .8), replace = TRUE))
# true_class <- sort(true_class)
# class1_probs <- rbeta(sum(true_class == "Class1"), 4, 1)
# class2_probs <- rbeta(sum(true_class == "Class2"), 1, 2.5)
# test_set <- data.frame(obs = true_class,Class1 = c(class1_probs, class2_probs))
# test_set$Class2 <- 1 - test_set$Class1
# test_set$pred <- factor(ifelse(test_set$Class1 >= .5, "Class1", "Class2"))
# # calculate the confusion matrix
# cm <- confusionMatrix(data = test_set$pred, reference = test_set$obs)
# 
# 

#Creo que es hasta aqui

# # https://www.guru99.com/r-random-forest-tutorial.html
# prediction <-predict(fit_rf, data_test)
# confusionMatrix(prediction, data_test$survived)
# 
# 
# library(e1071)
# 
# set.seed(1234)
# # Define the control
# trControl <- trainControl(method = "cv",
#                           number = 10,
#                           search = "grid")
# ## Random forest
# # x <- tDB[,1:(ncol(tDB)-1)]
# # y <- tDB[, ncol(tDB)]
# # rf <- randomForest(x[,seq(5,ncol(x)-1,1)], factor(y), ntree=500)
# data_train <- tDB
# data_train$GOLD <- as.factor(data_train$GOLD)
# # Run the model
# rf_default <- train(GOLD~.,
#                     data = data_train,
#                     method = "rf",
#                     metric = "Accuracy",
#                     trControl = trControl)
# # Print the results
# print(rf_default)
# 
# # test the model with values of mtry from 1 to 10 
# set.seed(1234)
# tuneGrid <- expand.grid(.mtry = c(1: 10))
# rf_mtry <- train(survived~.,
#                  data = data_train,
#                  method = "rf",
#                  metric = "Accuracy",
#                  tuneGrid = tuneGrid,
#                  trControl = trControl,
#                  importance = TRUE,
#                  nodesize = 14,
#                  ntree = 300)
# print(rf_mtry)
# 
# # create a loop to evaluate the different values of maxnodes
# store_maxnode <- list()
# tuneGrid <- expand.grid(.mtry = best_mtry)
# for (maxnodes in c(5: 15)) {
#   set.seed(1234)
#   rf_maxnode <- train(survived~.,
#                       data = data_train,
#                       method = "rf",
#                       metric = "Accuracy",
#                       tuneGrid = tuneGrid,
#                       trControl = trControl,
#                       importance = TRUE,
#                       nodesize = 14,
#                       maxnodes = maxnodes,
#                       ntree = 300)
#   current_iteration <- toString(maxnodes)
#   store_maxnode[[current_iteration]] <- rf_maxnode
# }
# results_mtry <- resamples(store_maxnode)
# summary(results_mtry)
# 
# # tune the number of trees
# store_maxtrees <- list()
# #for (ntree in c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)) {
# for (ntree in c(250)) {
#   set.seed(5678)
#   rf_maxtrees <- train(GOLD~.,
#                        data = data_train,
#                        method = "rf",
#                        metric = "Accuracy",
#                        #tuneGrid = tuneGrid,
#                        #trControl = trControl,
#                        importance = TRUE,
#                        #                       nodesize = 14,
#                        #                       maxnodes = 24,
#                        ntree = ntree)
#   key <- toString(ntree)
#   store_maxtrees[[key]] <- rf_maxtrees
# }
# results_tree <- resamples(store_maxtrees)
# summary(results_tree)
# 
