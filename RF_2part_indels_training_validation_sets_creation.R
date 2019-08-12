# VALIDATION AND TRAINING DATASETS CREATION

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

# CREACIÓN DE GRÁFICA PARA VER CUALES TIENEN EN COMUN CUANTOS TIENEN SOLOS ETC (FANCY VENN DIAGRAM)
u <- DNA_indel_final[, c("T_VAF_pla", "T_VAF_g4m2", "T_VAF_sva", "T_VAF_pin", "T_VAF_str", "T_VAF_lan")]
#Lo de arriba es dataframe solo con valores de VAF de esos programas
u[!is.na(u)] <- 1
u[is.na(u)] <- 0
colnames(u) <- c("platypus","mutect2", "svaba" ,"pindel","strelka", "lancet")
u <- data.frame(cbind(Identifier=DNA_indel_final$mut, u))
u$Identifier <- as.factor(u$Identifier)
u[,2:ncol(u)] <- sapply(u[,2:ncol(u)], function(x) as.numeric(as.character(x)))
upset(u, nsets = 6, sets = c( "platypus","mutect2", "svaba" ,"pindel","strelka", "lancet"), nintersects=NA, sets.bar.color = "#56B4E9", keep.order = F, empty.intersections = "on")

muts2 <- DNA_indel_final

#nrow(goldenT4)
#nrow(muts2)
# TP nrow(muts2[muts2$GOLD,])
muts2$GOLD <- ifelse(muts2$mut %in% nuevo_golden_indel_validated$mut, TRUE, FALSE)
muts2[is.na(muts2)] <- -1
#Faltaría añadir la columna de si esta en el golden o no para tenerlo todo en un único objeto

## Create training and validation series
db <- muts2
db_TP <- TP_DNA
db_TP[is.na(db_TP)] <- -1

db_FP <- FP
db_FP[is.na(db_FP)] <- -1

db_FN <- FN
db_FN[is.na(db_FN)] <- -1

#Una cosa si esta asignación es aleatoria cada vez, no será normal que cada vez de una cosa diferente?
# resulta que todo lo que se modifique se debe hacer antes de serparar los de test y training sin no el rf se lia
# Es decir todo, las eliminaciones de las clases categoricas, las conversiones en factores etc..
#Por otra parte creo que el código es muy redundante y en 3 liuneas podria estar tal y como explican aqui:
#https://stats.stackexchange.com/questions/235764/new-factors-levels-not-present-in-training-data


summary(db_TP)

summary(db_FP)
#Convieritendo en factor la columna de T_PROGS



df_db_TP <- as.data.frame(table(db_TP[,63]))

#Se van a eliminar 17 variantes TP, no pasa nada pq hay muchas y el RF es capaz de gestionar que en el 40% haya cosas que no estan en el 60% 
#además es importante que esta separación se haga después



#Vale ya esta todo claro parecía que no se eliminaban pero ya esta, lo que pasa es que el método de replace solo funciona con characteres
#no con factores, por tanto primero hay que convertirlo a factores, luego ver cuantos tipos hay con str, luego hacer dataframe y sacar aquellos tipos minoritarios
#Luego, reejecutar lo de arriba para que vuelva a la forma de characteres, luego hacer el replace, y luego volver a convertirlo a factor
#Luego comprobar con str que el numero de categorias es 52 en el caso de TP.


replace <- c("lancet,pindel,platypus",
             "pindel,platypus,svaba",
             "gatk4m2,lancet,pindel,platypus",
             "gatk4m2,pindel,platypus,svaba",
             "gatk4m2,lancet,platypus",
             "gatk4m2,lancet,pindel,platypus,svaba",
             "lancet,pindel,svaba") 

db_TP<-db_TP[!(db_TP$T_PROGS %in% replace),]
nrow(db_TP)  #11973
nrow(db_TP2) #11956
nrow(db_FP)


db_TP[,63] <- as.factor(db_TP[,63]) 
db_FP[,63] <- as.factor(db_FP[,63]) 
str(db_TP[,63]) #52 levels esto hay que bajarlo a 52 para que no de error nunca
str(db_FP[,63]) #33 levels



# df_db_TP <- as.data.frame(table(db_TP2[,63]))


# db_FP<-db_FP[!(db_FP$T_PROGS %in% replace),]
# summary(db_TP)
# 
# summary(db_FP)
# 
# db_TP[,63] <- as.factor(db_TP[,63]) 
# db_FP[,63] <- as.factor(db_FP[,63]) 

#Esto tiene sentido, solo TP y FP, las FN no pq no tienes la info de las variables es decir nunca podras aumentar el 
#Número de TP, ya que esto no es un variant caller pero básicamente pq no tienes esas variables.
#Pero y que hay de las TN, las que has identificado como no variante estando en lo correcto, para que las cosas que se
#parezcan a eso retirarles el status de variante?
#Supongo gque lo que quieres es reducir el numero de FP pero sin reducir el de TP, por eso le das ambas métricas



# DATASETS SEPARATION


vector_to_chose <- 1:nrow(db_TP)
number_of_items_to_chose <- nrow(db_TP)*1
training_TP <- sample(vector_to_chose, number_of_items_to_chose, replace = F)

vector_to_chose_FP <- 1:nrow(db_FP)
number_of_items_to_chose_FP <- nrow(db_FP)*1
training_FP <- sample(vector_to_chose_FP, number_of_items_to_chose_FP, replace = F)

length(training_TP) #7183 antes y ahora 7173, NUEVO 7191
length(training_FP) #66antes y ahora 66
training_FN <- sample(1:nrow(db_FN), nrow(db_FN)*0.6, replace = F)
# length(training_FN) #621 antes y ahora 621



#Al parecer el cálculo de los tamaños de los objetos puede deberse a la allocation en memoria de parte de ellos, esto podría explicar
# las diferencias de tamaños entre los mismo objetos con el mismo numero de rows y datos entre R terminal y Rstudio.
#Sin embargo como los datos son los mismos esto no debería alterar el resultado de las funciones aplicadas a posteriori.


ss <- db_TP[training_TP,]
ss$GOLD <- TRUE
nrow(ss) #7183
ncol(ss) #121
object.size(ss) #8659904 bytes

ss2 <-db_FP[training_FP,]
ss2$GOLD <- FALSE
nrow(ss2) #66
ncol(ss2) # 121
object.size(ss2) #103840 bytes

tDB1 <- rbind(ss,ss2)

#Aquí ya hay que hacer cambios para que sea igual que lo que hay en el de SNVs

#Renombramos el mismo objeto para que todo lo de los nombres siga siendo igual
#En principio es quitar todas las normales menos las de DP, el genotipo, y las variables no numéricas excepto t_progs
#Sería como solo dejar las col_number pero aqui no despues, que es como está hecho. 
#N_GT, T_GT, T_Filter_str

col_names <- c("CHROM", "POSITION",  "REF",   "ALT", "N_DP_avg","N_DP_std","T_NPROGS","T_PROGS","T_DP_avg",
               "T_DP_std","T_VAF_avg","T_MBQ_g4m2","T_MMQ_g4m2","T_MQ_str",
               "T_MQ0_str","T_SomaticEVS_str","T_RC_str","T_IC_str","T_IHP_str",
               "T_FETS_lan","T_SB_lan","T_QUAL_pin","T_S1_pin","T_S2_pin","T_GQ_pla",
               "T_GOF_pla","T_MMLQ_pla","T_HP_pla","T_MQ_pla","T_QD_pla","T_LO_sva",
               "T_LR_sva","T_MAPQ_sva","T_LOD_sva","mut", "GOLD")

col_numbers <- which(names(tDB1)%in%col_names)



tDB <- tDB1[,col_numbers]



# replace <- c("lancet,platypus","gatk4m2,lancet,pindel,platypus", "gatk4m2,lancet,pindel,platypus,svaba", "gatk4m2,lancet,platypus", "gatk4m2,pindel,platypus,svaba")
# tDB<-tDB[!(tDB$T_PROGS %in% replace),]
# tDB[,8] <- as.factor(x[,8])

# 7242 + 4846 60%y 40% 
# 7242/12088
# 4846/12088

# tDB1<-tDB1[!(tDB1$T_PROGS %in% replace),]

nrow(tDB) #7249
object.size(tDB)

table(tDB$GOLD)
# FALSE  TRUE 
# 66  7191 
#Esta asignación del GOLD no esta bien hecha, es decir da lo que tiene que dar pero esta mal. Debe ser hecha a nivel
#de TP y FP no aqui, pq da a entender cosas que no son. Todo el script esta hecho de esta forma que es nada reproducible
#y da lugar a dudas

# tDB$GOLD <- ifelse(tDB1$mut %in% nuevo_golden_indel_validated$mut, TRUE, FALSE)

#Ojo que esta tiene las rows del Golden evidentemente
tDB_FN <- db_FN[training_FN,]
#Vaya esta variable ya no se vuelve a usar pq no tiene info de variables de los callers


#Haciendo esto de meterle la GOLD antes cosa que no tiene el archivo original de Romina te ahorras tener que metersela despues
ss_val <- db_TP[-training_TP,]
ss_val$GOLD <- TRUE

ss2_val <-db_FP[-training_FP,]
ss2_val$GOLD <- FALSE

vDB2 <- rbind(ss_val,ss2_val)

nrow(vDB2)
nrow(db_FP[training_FP,])
nrow(db_FP[-training_FP,])
#Lo del db_TP[-training_TP,] es para coger el 40% restante, es decir para excluir todos aquellos que son de training, y solo coger los que
#conformaran los de validacion 

col_numbers <- which(names(vDB2)%in%col_names)
# AQUI vDB tiene 4839 rows
vDB <- vDB2[,col_numbers] #Simplemente este cambio ya hace que varie mucho el resultado, es decir, es evidente que lo otro esta mal, 
#no puedes usar los mismos para validación que de training. eso está en el script inicial de romina de snvs. Con lo cual mis resultados estan bien

# 21 TN   435 FN
# 23 FP  4795 TP

#Y los que se conseguian anteriormente eran:
# 61 TN   435 FN
#  5 TP  4795 TP



# replace <- c("lancet,platypus","gatk4m2,lancet,pindel,platypus", "gatk4m2,lancet,pindel,platypus,svaba", "gatk4m2,lancet,platypus", "gatk4m2,pindel,platypus,svaba")
# vDB<-vDB[!(vDB$T_PROGS %in% replace),]
# vDB1<-vDB1[!(vDB1$T_PROGS %in% replace),]
# vDB[,8] <- as.factor(vDB[,8])


vDB_FN <- db_FN[-training_FN,]
# vDB$GOLD <- ifelse(vDB$mut %in% nuevo_golden_indel_validated$mut, TRUE, FALSE)
###
# https://shiring.github.io/machine_learning/2017/03/31/webinar_code
###
bc_data <- db
## PCA
## features
library(tidyr)
#Esto hace una gráfica en la que se ven para cada variable sus FP y TP para ver si se encuentran variables que distinguen bien
#La forma de elección de las columnas esta hecha fisicamente y por tanto no reproducible, solo hay que coger solumnas con valores no con strings
# OJITO QUE PARA HACER LOS PLOTS SOLO SE DEBEN COGER COLUMNAS CON VALORES NO CON NAs, y ojito en como se conveirten esos NAs
#He introducido arriba un comando para que todas aquellas celdas con NAs les meta un -1 en muts2

#los de avg SUPONGO QUE ES MEJOR NO USARLOS YA QUE LOS NA hay no los entiendo aun su significado
colnames(bc_data)

# Cogemos las variables de los programas excluyendo la DP, VAF y las que tienen strings
# gather(bc_data, x, y, N_DP_avg,N_DP_std,T_NPROGS,T_DP_avg, T_DP_std,T_VAF_avg, T_MBQ_g4m2,T_MMQ_g4m2,T_MQ_str,T_MQ0_str,T_SomaticEVS_str,T_RC_str,T_IC_str,T_IHP_str,T_FETS_lan,T_SB_lan,T_QUAL_pin,T_S1_pin,T_S2_pin,T_GQ_pla,T_GOF_pla,T_MMLQ_pla,T_HP_pla,T_MQ_pla,T_QD_pla,T_LO_sva,T_LR_sva,T_MAPQ_sva,T_LOD_sva) %>%
#   ggplot(aes(x = y, color = GOLD, fill = GOLD)) +
#   geom_density(alpha = 0.3) +
#   facet_wrap( ~ x, scales = "free", ncol = 6)


## training vs validating datasets

# library(dplyr)
# 
# train_data <- tDB
# test_data <- vDB
# rbind(data.frame(group = "train", train_data),
#       data.frame(group = "test", test_data)) %>%
#   gather(x, y, N_DP_avg,N_DP_std,T_NPROGS,T_DP_avg, T_DP_std,T_VAF_avg, T_MBQ_g4m2,T_MMQ_g4m2,T_MQ_str,T_MQ0_str,T_SomaticEVS_str,T_RC_str,T_IC_str,T_IHP_str,T_FETS_lan,T_SB_lan,T_QUAL_pin,T_S1_pin,T_S2_pin,T_GQ_pla,T_GOF_pla,T_MMLQ_pla,T_HP_pla,T_MQ_pla,T_QD_pla,T_LO_sva,T_LR_sva,T_MAPQ_sva,T_LOD_sva) %>%
#   ggplot(aes(x = y, color = group, fill = group)) +
#   geom_density(alpha = 0.3) +
#   facet_wrap( ~ x, scales = "free", ncol = 6)
##############################################################
#Este segmento no esta en el incial pero lo que hay a continuación esta igual, osea que es lo mismo parece.
tDB$GOLD2 <- as.factor(tDB$GOLD)

tDB_nogold <- tDB[,1:(ncol(tDB)-2)]
tDB_only_VC_variables <- tDB_nogold[,seq(5,ncol(tDB_nogold)-1,1)]

tDB_golden_column <- tDB[, ncol(tDB)] # Esto lo que hace es coger la última columna que es la del GOLD y da un objeto que solo tienen 1 columna

tDB_only_VC_variables_and_golden <- tDB_nogold[,seq(5,ncol(tDB_nogold)-1,1)]
tDB_only_VC_variables_and_golden$GOLD2 <- as.factor(tDB$GOLD)




summary(tDB)
str(x)


#Por ejemplo aqui como sabe el RF que columnas no debe utilizar? Ah claro se le indica
#porque aqui le estas dando solo columnas adecuadas, no le das ninguna de posicion etc.. y las ultimas tampoco. 
#Ahora debe estar perfecto pq he limpiado el objeto


## Predict validation series:
x_val <- vDB[,1:(ncol(vDB)-1)]

# vDB2 <- rbind(ss_val,ss2_val)
# col_numbers <- which(names(vDB2)%in%col_names)
# vDB <- vDB2[,col_numbers]

vDB_only_VC_variables <- vDB[,seq(5,ncol(vDB)-1,1)]
vDB$GOLD2 <- as.factor(vDB$GOLD)
vDB_golden_colum <- vDB[, ncol(vDB)]
vDB_VCvariables_golden_only <- vDB_only_VC_variables
vDB_VCvariables_golden_only$GOLD2 <- vDB$GOLD2

colnames(vDB_only_VC_variables)
summary(x_val)

 
#############################################################

## Random forest
#OJITO AQUI QUE COLUMNAS HAY QUE INCLUIR????
x <- tDB[,1:(ncol(tDB)-1)]
object.size(x)

y <- tDB[, ncol(tDB)]
tDB$GOLD2 <- as.factor(tDB$GOLD)
#A continuación hay tb una selección de columnas que tampoco es reproducible... seq(5,ncol(x)-1,1)
#Deben escribirse las columnas que se quieren coger

#rf <- randomForest(x[,seq(5,ncol(x)-1,1)], factor(y), ntree=500)



#A partir de aquí mucho ojito cuando se use x, siempre usarla con x[col_numbers] y siempre ir haciendo pruebas para verificar los formatos con el de SNVs
#lo de la reproducibilidad del RF es algo que luego deberé mirar. 
#Al parecer es algo más complejo ya que https://stats.stackexchange.com/questions/120446/different-results-from-several-passes-of-random-forest-on-same-dataset
#describe el caso de que simplemente al copiar el codigo desde sublime da diferente, hay que analizar en que paso varia. 

#OJO METER T_PROGS AQUI O NO? Romina lo metió en el de SNVS asique pues tb lo meto y ya está.

# col_numbers <- which(names(x)%in%c("N_DP_avg","N_DP_std","T_NPROGS","T_PROGS","T_DP_avg"," T_DP_std","T_VAF_avg"," T_MBQ_g4m2","T_MMQ_g4m2","T_MQ_str","T_MQ0_str","T_SomaticEVS_str","T_RC_str","T_IC_str","T_IHP_str","T_FETS_lan","T_SB_lan","T_QUAL_pin","T_S1_pin","T_S2_pin","T_GQ_pla","T_GOF_pla","T_MMLQ_pla","T_HP_pla","T_MQ_pla","T_QD_pla","T_LO_sva","T_LR_sva","T_MAPQ_sva","T_LOD_sva"))
# object.size(col_numbers) #176 bytes
# object.size(x)
# 
# PRUEBA <- x[,col_numbers] #7249 rows

# PRUEBA <- x[,seq(5,ncol(x)-1,1)]
# #hay algo diferente este este PRUEBA y el otro de SNVS ya que en la función random forest aqui sale lo de NA y en el otro no. 
# summary(PRUEBA) #Gracias a esto he detectado que la diferencia estriba en que T_PROGS está como character aqui y en el otro no
# sapply(PRUEBA, class) # Con esto he detectado que en el otro (SNVs) esa columna es de .factor y aquí estava como character
# PRUEBA[,4] <- as.factor(PRUEBA[,4])
# summary(PRUEBA)

# summary(x)
# x[,8] <- as.factor(x[,8])
# 
# summary(x$T_PROGS)
# unique(x$T_PROGS) #57 Levels
# table(x$T_PROGS)
# sapply(x, class)
# 
# df <- as.data.frame(table(x$T_PROGS))


# set.seed(500)
# rf <- randomForest(x[,seq(5,ncol(x)-1,1)], factor(y), ntree=500)
# object.size(rf)

#Error in randomForest.default(x[, seq(5, ncol(x) - 1, 1)], factor(y),  : 
#Can not handle categorical predictors with more than 53 categories.

#Para solucionar esto,  https://stats.stackexchange.com/questions/157331/random-forest-predictors-have-more-than-53-categories
#O hacer lo de Jeremy Howards, o usar otro modelo, o mezclarlas todas en una categoria, (esto puede ser sucio pero como las variantes encontradas por las categorías menores son solo 1 o 2, pues eso)
#Además lo de juntarlos no se hasta que punto puede romper el modelo. Ya que de moemnto el número de T_NPROGS esta siempre correlacionado con el de T_PROGS
#Ya que cuando hay 4 en uno hay 4 en el otro no siempre los mismos 4 pero siguiendo esa dinámica. Por tanto no se hasta que punto, 
#Estaría el modelo bien ya que entre estas dos columnas hay una correlación más grande que entre todas las otras.

#Otra opción sería eliminar T_NPROGS, ya que la info que este recoje ya esta en la de T_PROGS, si realmente hay correlación entre
#ser encontrado por más programas para ser más veraz pues entonces eso quedará reflejado en la importancia de las variables, ya que esa será la primera, 
#Y los árboles harán la pregunta de OJOOOOOOO CLARO no se que pregunta hacen si habra más de 4 factores en esa columna o si estan unos fáctores específicos.
#¿Cómo trata el RF las columnas categóricas? No las puede tratar como números, y las de los números las puede llegar a tratar como intervalos o no?
#Royo si la variable X tiene valor superior a tal pero inferior a tal por un lado y el resto por otro. OTra cosa, el RF, las decisiones
#que hace son simepre binarias? O puede hacer decisiones de más posibilidades? todas estas preguntas son clave. 

#he pensado que por motivos de simplicidad podría elminarlas directamente, ya que como no todas las categorías estan representadas en este dataset
#sin eliminar nada ya se va a dar el caso de que nuevas ejecuciones de los programas encuentran combinaciones de variant callers
#no contempladas en el modelo por tanto el modelo pienso que ya debe estar preparado de antemano (por si solo) para lidiar con esto.
#Es por ello por lo que pienso que eliminar esas rows sería lo menos perjudicial. 
#Me tocaría eliminar 5 varintes que eso, pertenecen a las categorías de T_PROGS que menos variantes han encontrado, (1 y 2) por tanto
#la probabilidad de encontrar esa combinación de variant callers en la vida real es menor. Además como ya he mencionado anteriormente
#hay combinaciones no representadas en este dataset y no pasa nada.

#Entonces para eliminarlas podria decirle que elimine aquellas rows donde encunetre que en la columna de T_NPROGS coinciden esos valores


str(x)

# x1<-x
# 
# replace <- c("gatk4m2,lancet,pindel,platypus", "gatk4m2,lancet,pindel,platypus,svaba", "gatk4m2,lancet,platypus", "gatk4m2,pindel,platypus,svaba")
# x1<-x1[!(x1$T_PROGS %in% replace),]

# y1<-y
# y1<-y1[!(y1$T_PROGS %in% replace),]
