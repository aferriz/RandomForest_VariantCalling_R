# TCGA INDEL DATA PREPARATION FOR RF




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


#install.packages("export")
#install.packages("rgl")
#devtools::install_github('araastat/reprtree')
#install.packages("openssl")




## Read inputs
# TCGA datasets indels
#Read only 10 rows of file to check columns and so on
#combineFile_10 <- "/home/yo/Indels_RF_TCGA_10rows.tab"
#DNA1_10 <- fread(combineFile_10, header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))
#DNA1_10 <- read.table("/home/yo/Indels_RF_TCGA_10rows.tab", header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))



combineFile <- "/media/computer/HDD/MOUNT/INDELS_UCEC_CESC.tab"
DNA1 <- fread(combineFile, header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))
#Diferent numbers of rows obtained
#DNA1 <- read.table("/home/yo/INDELS_UCEC_CESC.tab", header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))

fileGoldenT4 <- "/media/computer/HDD/MOUNT/Golden_all"
goldenT4 <- fread(fileGoldenT4,header=TRUE)


table(goldenT4$Variant_Type)

a12 <- goldenT4[goldenT4$Variant_Type=="DEL"]
table(a12$Tumor_Sample_Barcode)
#goldenT4 <- read.table("/home/yo/Golden_all", header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))
colnames(goldenT4)




file_otro_golden_norm <- "/media/computer/HDD/MOUNT/other_golden_norm"
otro_golden_norm <- fread(file_otro_golden_norm,header=FALSE)
#otro_golden_norm <- read.table("/home/yo/other_golden_norm",header=FALSE)


colnames(otro_golden_norm) <- c("Chromosome"         ,         "Start_Position"     ,         "End_Position"      ,          "Variant_Classification"   ,   "Variant_Type"   ,            
                                "Reference_Allele"          ,  "Tumor_Seq_Allele1"        ,   "Tumor_Seq_Allele2"       ,    "Tumor_Sample_Barcode"   ,
                                "mutval_status_rna"         ,  "mutval_status_wex"       ,    "mutval_status_targeted" )                        

aa <- colnames(otro_golden_norm)


#bbb <- otro_golden_norm.join(golden_DEL_complex, otro_golden_norm.Tumor_Sample_Barcode == golden_DEL_complex.Tumor_Sample_Barcode, 'outer')


#joinedDF = df1.join(df2.select('A', 'E'), ['A'])
#neww <- merge.data.frame(otro_golden_norm, golden_DEL_complex, by="Tumor_Sample_Barcode" , all.x=TRUE )






DNA1 <- as.data.frame(DNA1)

#Here, I think that you havae to merge colums with _ CHROM and POS in order to create a new column to compare later
#Due to the fact that several sample are being used, it is also required to merge using T_SAMPLE
DNA1$mut <- paste0(DNA1$CHROM,"_",DNA1$POSITION,"_",DNA1$REF,"_",DNA1$ALT,"_",DNA1$T_SAMPLE)

# There is a challengue, and it is that INS, and DEL columns of the golden are not represented like the variant caller.tab
# I have thougth that it may be easier to transform the VC.tab columns to the format of the golden than doing the reverse option
# This way I would have all of them in the same format



goldenT4$mut <- paste0(goldenT4$Chromosome,"_",goldenT4$Start_position,"_",goldenT4$Reference_Allele,"_",goldenT4$Tumor_Seq_Allele2)

#### Processing starts to achieve same formats

#Separe indels and snv from the VC.tab file
DNA1_INS <- DNA1[which(DNA1[,4]> DNA1[,3]),]
DNA1_DEL <- DNA1[which(DNA1[,3]> DNA1[,4]),]

DNA1_INS <- DNA1[which(nchar(DNA1[,4])> nchar(DNA1[,3])),]
DNA1_DEL <- DNA1[which(nchar(DNA1[,4])< nchar(DNA1[,3])),]
DNA1_UM <- DNA1[which(nchar(DNA1[,4]) == nchar(DNA1[,3])),]

aa<-as.data.frame(DNA1_INS$REF)
bb<-as.data.frame(DNA1_DEL$REF)

summary(aa, maxsum= 299)
summary(bb, maxsum= 299) # You can check here the lack of normalization




golden_INS <- goldenT4[which(goldenT4[,5]== "INS"),]
golden_DEL <- goldenT4[which(goldenT4[,5]== "DEL"),]

prueba36 <- golden_INS[which(golden_INS[,13] == "unvalidated_powered"),] #52 variantes 1 no normalizada
prueba37 <- golden_DEL[which(golden_DEL[,13] == "unvalidated_powered"),] #218 variantes 1 no normalizada
#Inicialmente solo hay 270 variantes N. Con lo de normalización creo que me cargo la de INS no normalizada y quedan 269



#OJO LOS QUE NO SON NI DE UN TIPO NI DE OTRO PARA LO DE ARRIBA!! LOS QUE SON IGUALES!!
DNA1_UM <- DNA1[which(nchar(DNA1[,4]) == nchar(DNA1[,3])),]

GG<-as.data.frame(golden_INS$Reference_Allele)
VB<-as.data.frame(golden_DEL$Tumor_Seq_Allele2)

summary(GG, maxsum= 299)
summary(VB, maxsum= 299) #Esto da la prueba de que no están normalizados

golden_DEL_complex <- golden_DEL[which(golden_DEL[,8]!="-"),] # Aqui hay 1 unvalidated
#osea que realmente los números estarían casi bien, pq serían todo POSITIVES!!!! Pero no es así.
#LOS NUMEROS NO ESTAN CASI BIEN PQ NO ESTA HECHO LO DE LOS TN, POR TANTO LO QUE SI ESTA BIEN ESQUE DE TN HAY 271 +1 = 272 Y MAS 1 INS = 273
#PERO QUE DIRECTAMENTE SE ELIMINAN, ES DECIR TAL Y COMO SE HAYA CALCULADO LA PRECION DE LOS NPROGS HAY QUE CALCULAR LA PRECISION DEL RANDOM FOREST
#pero tener estas cifras en cuenta. 

#Entonces eso de momento 273 los N que tenemos, pero usé estos para calcular las métricas de NPROGS??




# resulta que de INS hay 2 compleas solo (es decir susceptibles de ser normalizadas, y encima es una validated y otra unvalidated)


golden_DEL_complex_sorted <- golden_DEL_complex[with(golden_DEL_complex, order(as.numeric(as.character(Chromosome))), as.numeric(as.character(Start_Position))), ]


#write.table(golden_DEL_complex, "/media/computer/HDD/MOUNT/golden_DEL_complex.txt", append = FALSE, sep = "\t", dec = ".",
#row.names = FALSE, col.names = FALSE, quote=FALSE)




# Here do you need to select N_ and T_? or only T_? The colum of T_SAMPLE only has Tumoral TCGA_IDs
# Indels variant callers: Mutec2, Lancet, Strelka, Pindel, Platypus, Svaba


# Feature selection, which variables from the programs are going to be used? Do they have to be the most differential ones?

# selectCols <- c("T_SAMPLE","CHROM","POSITION","REF","ALT","N_GT", "T_NPROGS","T_PROGS","mut",
#                 "N_DP_avg","N_DP_std", "T_DP_avg","T_DP_std","T_VAF_avg",
#                 "T_MBQ_g4m2","T_MMQ_g4m2",
#                 "T_FETS_lan","T_SB_lan",
#                 "T_MQ_str","T_MQ0_str","T_SomaticEVS_str"
# )

#List of all variables:

#"CHROM","POSITION","REF","ALT","N_GT", "T_NPROGS","T_PROGS",
# "N_DP_avg","N_DP_std","N_DP_REF_avg",           "N_DP_ALT_avg",           "N_VAF_avg",  "N_DP_REF_g4m2",          "N_DP_ALT_g4m2",          "N_VAF_g4m2",
# "N_MBQ_g4m2", "N_MMQ_g4m2", "N_IN_PON_g4m2",          "N_FILTER_g4m2",          "N_DP_REF_str",           "N_DP_ALT_str",          
# "N_VAF_str",  "N_MQ_str",   "N_MQ0_str",  "N_SomaticEVS_str",       "N_RU_str",   "N_RC_str",  
# "N_IC_str",   "N_IHP_str",  "N_FILTER_str",           "N_DP_REF_lan",           "N_DP_ALT_lan",           "N_VAF_lan", 
# "N_FETS_lan", "N_SB_lan",   "N_FILTER_lan",          "N_DP_REF_pin",           "N_DP_ALT_pin",           "N_VAF_pin",  "N_QUAL_pin", "N_S1_pin",              
# "N_S2_pin",         "N_FILTER_pin",     "N_DP_REF_pla",     "N_DP_ALT_pla",     "N_VAF_pla",        "N_GQ_pla",        
# "N_GOF_pla",        "N_MMLQ_pla",       "N_HP_pla",         "N_MQ_pla",         "N_QD_pla",         "N_FILTER_pla",    
# "N_DP_REF_sva",     "N_DP_ALT_sva",     "N_VAF_sva",        "N_LO_sva",         "N_LR_sva",         "N_MAPQ_sva",      
# "N_LOD_sva",        "N_REPSEQ_sva",     "N_FILTER_sva",     "T_GT",             "T_NPROGS",         "T_PROGS",         
# "T_DP_avg",         "T_DP_std",         "T_DP_REF_avg",     "T_DP_ALT_avg",     "T_VAF_avg",        "T_DP_REF_g4m2",   
# "T_DP_ALT_g4m2",    "T_VAF_g4m2",       "T_MBQ_g4m2",       "T_MMQ_g4m2",       "T_IN_PON_g4m2",    "T_FILTER_g4m2",   
# "T_DP_REF_str",     "T_DP_ALT_str",     "T_VAF_str",        "T_MQ_str",         "T_MQ0_str",        "T_SomaticEVS_str",
# "T_RU_str",         "T_RC_str",         "T_IC_str",         "T_IHP_str",        "T_FILTER_str",     "T_DP_REF_lan",    
# "T_DP_ALT_lan",     "T_VAF_lan",        "T_FETS_lan",       "T_SB_lan",         "T_FILTER_lan",     "T_DP_REF_pin",     "T_DP_ALT_pin",     "T_VAF_pin",       
# "T_QUAL_pin",       "T_S1_pin",         "T_S2_pin",         "T_FILTER_pin",     "T_DP_REF_pla",     "T_DP_ALT_pla",    
# "T_VAF_pla",        "T_GQ_pla",         "T_GOF_pla",        "T_MMLQ_pla",       "T_HP_pla",         "T_MQ_pla",        
# "T_QD_pla",         "T_FILTER_pla",     "T_DP_REF_sva",     "T_DP_ALT_sva",     "T_VAF_sva",        "T_LO_sva",        
# "T_LR_sva",         "T_MAPQ_sva",       "T_LOD_sva",        "T_REPSEQ_sva",     "T_FILTER_sva",     "T_SAMPLE" 

#The mut variable is added afterwards so leave it there

selectCols <- c("mut","CHROM","POSITION","REF","ALT","N_GT", "T_NPROGS","T_PROGS",
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
                "T_LR_sva",         "T_MAPQ_sva",       "T_LOD_sva",        "T_REPSEQ_sva",     "T_FILTER_sva",     "T_SAMPLE"
)





muts2 <- DNA1[,names(DNA1) %in% selectCols]

# DNA1_10 <- as.data.frame(DNA1_10)
# muts2_10 <- DNA1_10[,names(DNA1_10) %in% selectCols]

muts2$T_PROGS <- as.factor(muts2$T_PROGS)


muts2$GOLD <- ifelse(muts2$mut %in% goldenT4$mut, TRUE, FALSE)

muts2[is.na(muts2)] <- -1
head(muts2)
sapply(muts2, class)
sapply(muts2, typeof)
table(muts2$T_PROGS)













#FINALLY!!
otro_golden_norm$Matched_Norm_Sample_Barcode <- golden_DEL_complex$Matched_Norm_Sample_Barcode[match(otro_golden_norm$Tumor_Sample_Barcode, golden_DEL_complex$Tumor_Sample_Barcode)]
#Esto para solo quedarte con las que no estan en golde_DEL_complex pq luego voy a adicionarlas normalizadas
golden_DEL_complex$mut <- NULL
xxxxx <- anti_join(goldenT4, golden_DEL_complex)

zxzx <- data.frame(otro_golden_norm[otro_golden_norm$Variant_Type =="DEL",])
yyy <- data.frame(otro_golden_norm[otro_golden_norm$Variant_Type =="SNV",])
zzzzz <- data.frame(otro_golden_norm[otro_golden_norm$Variant_Type =="INS",])
zzzzz$Tumor_Seq_Allele1 <- "-"
zzzzz$Reference_Allele <- "-"


zxzx$Start_Position <- zxzx$Start_Position+1
zxzx$Tumor_Seq_Allele2 <- "-"

zxzx$Reference_Allele <- sub('.', '', zxzx$Reference_Allele)


xxxxx$mut <- NULL

nuevo_golden <- rbind(xxxxx,zxzx)
nuevo_golden <- rbind(nuevo_golden,yyy)
nuevo_golden <- rbind(nuevo_golden,zzzzz) 
nuevo_golden2 <- data.frame(lapply(nuevo_golden, function(x) {
  gsub("SNP", "SNV", x)
}))


#Antes de paste, separar por SNV e indels, y convertir con lo del -.

#CON ESTO YA ESTARIA EL GOLDEN LISTO!!!!!!!!
nuevo_golden2$mut <- paste0(nuevo_golden2$Chromosome,"_",nuevo_golden2$Start_Position,"_",nuevo_golden2$Reference_Allele,"_",nuevo_golden2$Tumor_Seq_Allele2,"_",nuevo_golden2$Tumor_Sample_Barcode)


#write.table(nuevo_golden2, "/media/computer/HDD/MOUNT/nuevo_golden2.txt", append = FALSE, sep = "\t", dec = ".",
#  row.names = FALSE, col.names = FALSE, quote=FALSE)


nuevo_golden_ins <- nuevo_golden2[which(nuevo_golden2[,5] == "INS" ),]
nuevo_golden_del <- nuevo_golden2[which(nuevo_golden2[,5] == "DEL" ),]
nuevo_golden_indel <- rbind(nuevo_golden_ins,nuevo_golden_del)


rownames(nuevo_golden_indel) <- c()

nuevo_golden_indel[13301,13] = factor("unvalidated_powered")
nuevo_golden_indel[13302,13] = factor("unvalidated_powered")

# Ahora sí. Por tanto para medir el recall y specificity pillar solo los unvalidated powered, y para recall solo los diferentes de unvalidated powered, que son validated powered
#NOS QUEDAMOS POR AQUI 3/7/2019 ESTO RECALCULAR MÉTRICAS TENIENDO EN CUENTA LO DE ARRIBA


#En este golden, todo es validated_powered excepto las que pone unvalidated_powered y las dos del sample sin tag TCGA-A5-A0GB-10A-01W-A062-09
#Posiciones 237024373 y 237024398 (Estas son unvalidated, es decir que no deben estar y si no hago nada se evaluan como que si que deben estar)


prueba34 <- nuevo_golden_indel[which(nuevo_golden_indel[,13] == "unvalidated_powered"),] #269 variantes

# write.table(nuevo_golden_indel, "/media/computer/HDD/MOUNT/nuevo_golden_indel.txt", append = FALSE, sep = "\t", dec = ".",
#             row.names = FALSE, col.names = FALSE, quote=FALSE)



#Ahora a preparar los .tab, hay que pasarlos a -. Y no estan normalizados, por tanto...
#1) pillar todas las del que tienen mas de 1 nt en ALT.
#2) quitarlas de esa tabla
#3) hacer dummy con ellas y normalizarlas
#4) pegar debajo en la misma tabla
#5) ordenar como pueda, ya que el tab no esta en vcf y no se puede usar vcf sort, pero con un for se puede hacer fácil


#1)
Bad_DEL <- DNA1_DEL[which(nchar(DNA1_DEL[,4])> 1),]

#2)
DNA1_DEL_good <- anti_join(DNA1_DEL, Bad_DEL)

#3) Esto lo hago en python, sería sacar del tab la info necesaria para el VCF, y correr pre.py. Pero ese paso de sacar la info necesaria, 
#En principio no se la forma de guardar la info para cada variante y que al normalizarse se conserve en cada nueva variante que se crea, no se si eso es
#automático o no, Y de todas maneras tendra que construir un VCF con todas estas features, metiendolas en la parte de sample, separados por dos puntos


# Count number of occurences in R
table(Bad_DEL$T_SAMPLE)
# Con esto se ve que las Bad_DEL corresponden a muchos samples no solo a unos cuantos y que por tanto no es una gran perdida eliminarlas teneinedo en cuenta que cada sample tiene muchas más
# Por tanto no se utilizan y ya esta luego siempre se puede volver a estas. Por ser eficaz.

DNA1_DEL_good_10 <- head(DNA1_DEL_good)
#En principio es fácil  de estas cambiar el ALT por un - y eliminar un nt del REF, por otra parte, sumarle una posición a POS


DNA1_DEL_good_modified <- DNA1_DEL_good
DNA1_DEL_good_modified$POSITION <- DNA1_DEL_good$POSITION +1
DNA1_DEL_good_modified$ALT <- ("-")
DNA1_DEL_good_modified$REF <- sub('.', '', DNA1_DEL_good_modified$REF)
DNA1_DEL_good_modified$mut <- paste0(DNA1_DEL_good_modified$CHROM,"_",DNA1_DEL_good_modified$POSITION,"_",DNA1_DEL_good_modified$REF,"_",DNA1_DEL_good_modified$ALT,"_",DNA1_DEL_good_modified$T_SAMPLE)


#esto deberia arreglarlo y ordenarlo para que estuviese todo reploducidble y bien en caso de que pase algo 


selectCols <- c("mut","CHROM","POSITION","REF","ALT","N_GT", "T_NPROGS","T_PROGS",
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
                "T_LR_sva",         "T_MAPQ_sva",       "T_LOD_sva",        "T_REPSEQ_sva",     "T_FILTER_sva",     "T_SAMPLE"
)

table(DNA1_INS$REF)
# A     C     G     T 
# 8615 11322 11440 10036 

DNA1_INS_good_modified <- DNA1_INS
DNA1_INS_good_modified$POSITION <- DNA1_INS_good_modified$POSITION
DNA1_INS_good_modified$REF <- ("-")
DNA1_INS_good_modified$ALT <- sub('.', '', DNA1_INS_good_modified$ALT)
DNA1_INS_good_modified$mut <- paste0(DNA1_INS_good_modified$CHROM,"_",DNA1_INS_good_modified$POSITION,"_",DNA1_INS_good_modified$REF,"_",DNA1_INS_good_modified$ALT,"_",DNA1_INS_good_modified$T_SAMPLE)






DNA_indel <- rbind(DNA1_INS_good_modified, DNA1_DEL_good_modified)

DNA_indel_final <- DNA_indel[,names(DNA_indel) %in% selectCols]
#Ahora ya si que si que se puede comenzar

# write.table(DNA_indel_final, "/media/computer/HDD/MOUNT/DNA_indel_final.txt",append = FALSE, sep = "\t", dec = ".",
#             row.names = FALSE, col.names = FALSE, quote=FALSE)

# Ahora ya tenemos los dos archivos preparados para su comparación y puedo calcular las métricas para NPROGS y para cada programa por separado

# Convendria ordenar el script antes de hacer esto o lo ordeno justo antes de pasar a hacer el RF?

DNA_indel_final_10 <- head(DNA_indel_final)


DNA_indel_final$chr_pos <- paste0(DNA_indel_final$CHROM,"_",DNA_indel_final$POSITION,"_",DNA_indel_final$T_SAMPLE)


# T_NPROGS al menos de 1 (es decir de 1 a 5, todas las rows)
# En principio es hacer cuantos de la tabla de nuevo_golden_indel$mut estan presentes en esta tabla?


#OJOOOOOOOOO
#Falta hacer subset a los targeted validated powered, ya que aqui hay tanto TP como TN, ojito!!!!!



nuevo_golden_indel # 13378 MUTACIONES



#Ahora meter a todos los guiones del golden en la columna de targeted validated_powered
nuevo_golden_indel[13] <- data.frame(lapply(nuevo_golden_indel[13], function(x) {gsub("-", "validated_powered", x)}))



nuevo_golden_indel_validated <- nuevo_golden_indel[nuevo_golden_indel$mutval_status_targeted == "validated_powered",] # 13107 indels antes 13009
#Pues supongo que ahora es recalcular todo de nuevo ejecutar todo y au

# Se crean los dos tipos de TP, con datos del golden o bien con datos de los variant callers, tb dejo un objeto TP para no cambiar el resto que es como el TP_golden

TP_golden <- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% DNA_indel_final$mut,]
TP <- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% DNA_indel_final$mut,]
TP_DNA <- DNA_indel_final[DNA_indel_final$mut %in% nuevo_golden_indel_validated$mut,]

13009 - 1036
1-0.5904

(0.5904*111)/0.4096
FN <- anti_join(nuevo_golden_indel_validated, TP) # 1036 -- Ahora 1086

Recall_otraforma <- nrow(TP)/(nrow(FN)+nrow(TP)) *100 # 91.71435


Recall_metric <- nrow(TP)/nrow(nuevo_golden_indel_validated) *100
Recall_metric # 92.03628     Ahora = 91.71435


# Vale pues sigo calculando el resto de Recalls para las otras posibilidades. 
# NPROGS igual a 2 o más
progs_2_3_4_5_6 <- c("2","3","4","5","6")
indel_2_or_more <- DNA_indel_final[DNA_indel_final$T_NPROGS %in% progs_2_3_4_5_6,]
golden_in_indel_2_or_more<- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% indel_2_or_more$mut,]
Recall_2_or_more <- nrow(golden_in_indel_2_or_more)/nrow(nuevo_golden_indel_validated) *100
Recall_2_or_more # 83.59597   Ahora 83.28374
FN_2_or_more <- anti_join(nuevo_golden_indel_validated, golden_in_indel_2_or_more) #2191
nrow(FN_2_or_more)
TN_2_or_more <- 
  
  # NPROGS igual a 3 o más
  progs_3_4_5_6 <- c("3","4","5","6")
indel_3_or_more <- DNA_indel_final[DNA_indel_final$T_NPROGS %in% progs_3_4_5_6,]
golden_in_indel_3_or_more<- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% indel_3_or_more$mut,]
Recall_3_or_more <- nrow(golden_in_indel_3_or_more)/nrow(nuevo_golden_indel_validated) *100
Recall_3_or_more # 75.10185 Ahora 74.81498
FN_3_or_more <- anti_join(nuevo_golden_indel_validated, golden_in_indel_3_or_more) #3301
nrow(FN_3_or_more)


# NPROGS igual a 4 o más
progs_4_5_6 <- c("4","5","6")
indel_4_or_more <- DNA_indel_final[DNA_indel_final$T_NPROGS %in% progs_4_5_6,]
golden_in_indel_4_or_more<- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% indel_4_or_more$mut,]
Recall_4_or_more <- nrow(golden_in_indel_4_or_more)/nrow(nuevo_golden_indel_validated) *100
Recall_4_or_more # 64.97041 Ahora  64.65248
FN_4_or_more <- anti_join(nuevo_golden_indel_validated, golden_in_indel_4_or_more) #4633
nrow(FN_4_or_more)


# NPROGS igual a 5 o más
progs_5_6 <- c("5","6")
indel_5_or_more <- DNA_indel_final[DNA_indel_final$T_NPROGS %in% progs_5_6,]
golden_in_indel_5_or_more<- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% indel_5_or_more$mut,]
Recall_5_or_more <- nrow(golden_in_indel_5_or_more)/nrow(nuevo_golden_indel_validated) *100
Recall_5_or_more # 49.99616   Ahora 49.74441
FN_5_or_more <- anti_join(nuevo_golden_indel_validated, golden_in_indel_5_or_more) #6587
nrow(FN_5_or_more)

# NPROGS igual a 6
indel_6 <- DNA_indel_final[DNA_indel_final$T_NPROGS %in% "6",]
golden_in_indel_6<- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% indel_6$mut,]
Recall_6 <- nrow(golden_in_indel_6)/nrow(nuevo_golden_indel_validated) *100
Recall_6 # 24.337 Ahora 24.2237
FN_6 <- anti_join(nuevo_golden_indel_validated, golden_in_indel_6) #9932
nrow(FN_6)


#Recall para los programas sueltos: gatk4m2,lancet,pindel,platypus,strelka,svaba

indel_gatk4m2 <- DNA_indel_final[DNA_indel_final$T_PROGS %like% "gatk4m2",]
table(indel_gatk4m2$T_PROGS) #Check - perfect
golden_in_indel_gatk4m2<- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% indel_gatk4m2$mut,]
Recall_gatk4m2 <- nrow(golden_in_indel_gatk4m2)/nrow(nuevo_golden_indel_validated) *100
Recall_gatk4m2 # 70.06688    69.76425
FN_gatk4m2 <- anti_join(nuevo_golden_indel_validated, golden_in_indel_gatk4m2) #3963
nrow(FN_gatk4m2)


indel_lancet <- DNA_indel_final[DNA_indel_final$T_PROGS %like% "lancet",]
table(indel_lancet$T_PROGS) #Check - perfect
golden_in_indel_lancet<- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% indel_lancet$mut,]
Recall_lancet <- nrow(golden_in_indel_lancet)/nrow(nuevo_golden_indel_validated) *100
Recall_lancet # 48.78161   48.58473
FN_lancet <- anti_join(nuevo_golden_indel_validated, golden_in_indel_lancet) #6739
nrow(FN_lancet)

indel_pindel <- DNA_indel_final[DNA_indel_final$T_PROGS %like% "pindel",]
table(indel_pindel$T_PROGS) #Check - perfect
golden_in_indel_pindel<- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% indel_pindel$mut,]
Recall_pindel <- nrow(golden_in_indel_pindel)/nrow(nuevo_golden_indel_validated) *100
Recall_pindel # 65.42394  65.10262
FN_pindel <- anti_join(nuevo_golden_indel_validated, golden_in_indel_pindel) #4574
nrow(FN_pindel)

indel_platypus <- DNA_indel_final[DNA_indel_final$T_PROGS %like% "platypus",]
table(indel_platypus$T_PROGS) #Check - perfect
golden_in_indel_platypus<- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% indel_platypus$mut,]
Recall_platypus <- nrow(golden_in_indel_platypus)/nrow(nuevo_golden_indel_validated) *100
Recall_platypus # 51.41056   51.1635
FN_platypus <- anti_join(nuevo_golden_indel_validated, golden_in_indel_platypus) #6401
nrow(FN_platypus)

indel_strelka <- DNA_indel_final[DNA_indel_final$T_PROGS %like% "strelka",]
table(indel_strelka$T_PROGS) #Check - perfect
golden_in_indel_strelka<- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% indel_strelka$mut,]
Recall_strelka <- nrow(golden_in_indel_strelka)/nrow(nuevo_golden_indel_validated) *100
Recall_strelka # 86.4171   86.11429
FN_strelka <- anti_join(nuevo_golden_indel_validated, golden_in_indel_strelka) #1820
nrow(FN_strelka)

indel_svaba <- DNA_indel_final[DNA_indel_final$T_PROGS %like% "svaba",]
table(indel_svaba$T_PROGS) #Check - perfect
golden_in_indel_svaba<- nuevo_golden_indel_validated[nuevo_golden_indel_validated$mut %in% indel_svaba$mut,]
Recall_svaba <- nrow(golden_in_indel_svaba)/nrow(nuevo_golden_indel_validated) *100
Recall_svaba # 67.93758   67.70428
FN_svaba <- anti_join(nuevo_golden_indel_validated, golden_in_indel_svaba) #4233
nrow(FN_svaba)

row_names <- c("1 or more","2 or more","3 or more","4 or more","5 or more","6","gatk4m2","lancet","pindel","platypus","strelka","svaba")
recall_values <- c(Recall_metric,Recall_2_or_more,Recall_3_or_more, Recall_4_or_more,Recall_5_or_more,Recall_6,Recall_gatk4m2,Recall_lancet,Recall_pindel,Recall_platypus,Recall_strelka,Recall_svaba)
METRICS <- data.table(row_names,recall_values)
METRICS

#PRECISION
#Aquí solo se pueden evaluar aquellos del golden que son unvalidated_powered que son muy poquitos
#Además el approach es distinto, aquí lo que hay que hacer es: De todas las posiciones donde no debe haber sido encontrada ninguna mutación (solo las posiciones CHROM y POS, cuantas han sido recogidas en el .tab? eso dirá el número de FP)

nuevo_golden_indel_unvalidated <- nuevo_golden_indel[nuevo_golden_indel$mutval_status_targeted == "unvalidated_powered",] # 269 indels
#A esta habria que sumarle 269 +1 que es normalizado que no esta con tag, y de inserciones tb, otra más, y eso da 271




#Ahora hay que crear otras columnas solo basadas en CHROM y POS, ya que se supone que es esto lo que vamos a comparar. 
nuevo_golden_indel_unvalidated$chr_pos <- paste0(nuevo_golden_indel_unvalidated$Chromosome,"_",nuevo_golden_indel_unvalidated$Start_Position,"_",nuevo_golden_indel_unvalidated$Tumor_Sample_Barcode)
# nuevo_golden_indel_unvalidated serían los N
DNA_indel_final$chr_pos <- paste0(DNA_indel_final$CHROM,"_",DNA_indel_final$POSITION,"_",DNA_indel_final$T_SAMPLE)

FP <- DNA_indel_final[DNA_indel_final$mut %in% nuevo_golden_indel_unvalidated$mut,] #110 ahora sí ya que los TN tb son de cada sample -- Ahora 112

FP2 <- nuevo_golden_indel_unvalidated[nuevo_golden_indel_unvalidated$mut %in% DNA_indel_final$mut,] #110 ahora sí ya que los TN tb son de cada sample -- Ahora 112


TN <- anti_join(nuevo_golden_indel_unvalidated, FP) # 161
#OJO que la forma correcta de calcular los TN no es la de arriba, es esta: COMPROBACIÓN DE QUE ESTÁ BIEN
TN_correcta <- nuevo_golden_indel_unvalidated[!(nuevo_golden_indel_unvalidated$mut %in% DNA_indel_final$mut),] #161


#ojito que las métricas estan totalmente biased por el numero de P Y N, ya que si el numero de N, es muy bajo el posible número de FP será super bajo tb.
#Lo mismo con specificity va mejor aqui

Specificity_metric <- nrow(TN) / (nrow(TN) + nrow(FP)) *100
Specificity_metric # 59.04059

Precision_metric <- nrow(TP)/(nrow(TP) +nrow(FP)) *100
Precision_metric # 99.08143
#Esta disparidad ocurre pq de todos los P que hemos detectado, solo podemos decir sobre unos pocos, si son verdaderos o falsos, el resto no lo sabemos y quedan en un limbo.

nrow(DNA_indel_final)
# De los 232029 que has detectado como P, solo puedes evaluar el estatus de 13378 (nrow(nuevo_golden_indel) (real mut + real places without mut))

# Entonces en este caso es mejor usar Specificity, o checkear que todo sea correcto y no haya que hacer ninguna operación para poder tratar con estos datos
# Esto es importante y puede suponer una nueva fuente de aprendizaje


#NPROGS igual a 2 o más
progs_2_3_4_5_6 <- c("2","3","4","5","6")
indel_2_or_more <- DNA_indel_final[DNA_indel_final$T_NPROGS %in% progs_2_3_4_5_6,]
FP_2 <- indel_2_or_more[indel_2_or_more$mut %in% nuevo_golden_indel_unvalidated$mut,] #79
TN_2 <- nuevo_golden_indel_unvalidated[!(nuevo_golden_indel_unvalidated$mut %in% indel_2_or_more$mut),] #192
nrow(TN_2)
Specificity_metric_2_bien <- nrow(TN_2) / (nrow(TN_2)+nrow(FP_2)) *100 # 79.7048
Specificity_metric_2_bien
Specificity_metric_2_or_more <- nrow(TN) / (nrow(TN) + nrow(FP_2)) *100
Specificity_metric_2_or_more # 67.22689

Precision_metric_2_or_more <- nrow(TP)/(nrow(TP) +nrow(FP_2)) *100
Precision_metric_2_or_more # 99.35275



# NPROGS igual a 3 o más
progs_3_4_5_6 <- c("3","4","5","6")
indel_3_or_more <- DNA_indel_final[DNA_indel_final$T_NPROGS %in% progs_3_4_5_6,]
FP_3 <- indel_3_or_more[indel_3_or_more$mut %in% nuevo_golden_indel_unvalidated$mut,] # 55
TN_3 <- nuevo_golden_indel_unvalidated[!(nuevo_golden_indel_unvalidated$mut %in% indel_3_or_more$mut),] #216
nrow(TN_3)
Specificity_metric_3_bien <- nrow(TN_3) / (nrow(TN_3)+nrow(FP_3)) *100 # 79.7048
Specificity_metric_3_bien
Specificity_metric_3_or_more <- nrow(TN) / (nrow(TN) + nrow(FP_3)) *100
Specificity_metric_3_or_more # 74.4186

Precision_metric_3_or_more <- nrow(TP)/(nrow(TP) +nrow(FP_3)) *100
Precision_metric_3_or_more # 99.54273



# NPROGS igual a 4 o más
progs_4_5_6 <- c("4","5","6")
indel_4_or_more <- DNA_indel_final[DNA_indel_final$T_NPROGS %in% progs_4_5_6,]
FP_4 <- indel_4_or_more[indel_4_or_more$mut %in% nuevo_golden_indel_unvalidated$mut,] # 42
TN_4 <- nuevo_golden_indel_unvalidated[!(nuevo_golden_indel_unvalidated$mut %in% indel_4_or_more$mut),] #229
nrow(TN_4)
Specificity_metric_4_bien <- nrow(TN_4) / (nrow(TN_4)+nrow(FP_4)) *100 # 84.50185
Specificity_metric_4_bien

Specificity_metric_4_or_more <- nrow(TN) / (nrow(TN) + nrow(FP_4)) *100
Specificity_metric_4_or_more # 79.20792

Precision_metric_4_or_more <- nrow(TP)/(nrow(TP) +nrow(FP_4)) *100
Precision_metric_4_or_more # 99.65044

# NPROGS igual a 5 o más
progs_5_6 <- c("5","6")
indel_5_or_more <- DNA_indel_final[DNA_indel_final$T_NPROGS %in% progs_5_6,]
FP_5 <- indel_5_or_more[indel_5_or_more$mut %in% nuevo_golden_indel_unvalidated$mut,] # 26
TN_5 <- nuevo_golden_indel_unvalidated[!(nuevo_golden_indel_unvalidated$mut %in% indel_5_or_more$mut),] #245
nrow(TN_5)

Specificity_metric_5_bien <- nrow(TN_5) / (nrow(TN_5)+nrow(FP_5)) *100 # 90.4059
Specificity_metric_5_bien


Specificity_metric_5_or_more <- nrow(TN) / (nrow(TN) + nrow(FP_5)) *100
Specificity_metric_5_or_more # 86.02151

Precision_metric_5_or_more <- nrow(TP)/(nrow(TP) +nrow(FP_5)) *100
Precision_metric_5_or_more # 99.78332


# NPROGS igual a 6
indel_6 <- DNA_indel_final[DNA_indel_final$T_NPROGS %in% "6",]
FP_6 <- indel_6[indel_6$mut %in% nuevo_golden_indel_unvalidated$mut,] # 5
TN_6 <- nuevo_golden_indel_unvalidated[!(nuevo_golden_indel_unvalidated$mut %in% indel_6$mut),] #266
nrow(TN_6)

Specificity_metric_6_bien <- nrow(TN_6) / (nrow(TN_6)+nrow(FP_6)) *100 # 98.15498


# Specificity_metric_6 <- nrow(TN) / (nrow(TN) + nrow(FP_6)) *100
# Specificity_metric_6 # 96.9697

Precision_metric_6 <- nrow(TP)/(nrow(TP) +nrow(FP_6)) *100
Precision_metric_6 # 99.95826


#Precision y specificity para los programas sueltos: gatk4m2,lancet,pindel,platypus,strelka,svaba

indel_gatk4m2 <- DNA_indel_final[DNA_indel_final$T_PROGS %like% "gatk4m2",]
table(indel_gatk4m2$T_PROGS) #Check - perfect
FP_gatk4m2 <- indel_gatk4m2[indel_gatk4m2$mut %in% nuevo_golden_indel_unvalidated$mut,] # 64
Specificity_metric_gatk4m2 <- nrow(TN) / (nrow(TN) + nrow(FP_gatk4m2)) *100
Specificity_metric_gatk4m2 # 71.11111
Precision_metric_gatk4m2 <- nrow(TP)/(nrow(TP) +nrow(FP_gatk4m2)) *100
Precision_metric_gatk4m2 # 99.46004
TN_gatk4m2 <- nuevo_golden_indel_unvalidated[!(nuevo_golden_indel_unvalidated$mut %in% indel_gatk4m2$mut),] # 207
nrow(TN_gatk4m2)
Specificity_metric_gatk4m2_bien <- nrow(TN_gatk4m2) / (nrow(TN_gatk4m2)+nrow(FP_gatk4m2)) *100 # 75.82418
Specificity_metric_gatk4m2_bien

indel_lancet <- DNA_indel_final[DNA_indel_final$T_PROGS %like% "lancet",]
table(indel_lancet$T_PROGS) #Check - perfect
FP_lancet <- indel_lancet[indel_lancet$mut %in% nuevo_golden_indel_unvalidated$mut,] #33
Specificity_metric_lancet <- nrow(TN) / (nrow(TN) + nrow(FP_lancet)) *100
Specificity_metric_lancet # 82.90155
Precision_metric_lancet <- nrow(TP)/(nrow(TP) +nrow(FP_lancet)) *100
Precision_metric_lancet # 99.72514
TN_lancet <- nuevo_golden_indel_unvalidated[!(nuevo_golden_indel_unvalidated$mut %in% indel_lancet$mut),] # 238
nrow(TN_lancet)
Specificity_metric_lancet_bien <- nrow(TN_lancet) / (nrow(TN_lancet)+nrow(FP_lancet)) *100 # 87.82288
Specificity_metric_lancet_bien




indel_pindel <- DNA_indel_final[DNA_indel_final$T_PROGS %like% "pindel",]
table(indel_pindel$T_PROGS) #Check - perfect
FP_pindel <- indel_pindel[indel_pindel$mut %in% nuevo_golden_indel_unvalidated$mut,] # 52
TN_pindel <- nuevo_golden_indel_unvalidated[!(nuevo_golden_indel_unvalidated$mut %in% indel_pindel$mut),] # 219
nrow(TN_pindel)
Specificity_metric_pindel_bien <- nrow(TN_pindel) / (nrow(TN_pindel)+nrow(FP_pindel)) *100 # 80.81181
Specificity_metric_pindel_bien
Specificity_metric_pindel <- nrow(TN) / (nrow(TN) + nrow(FP_pindel)) *100
Specificity_metric_pindel # 75.4717
Precision_metric_pindel <- nrow(TP)/(nrow(TP) +nrow(FP_pindel)) *100
Precision_metric_pindel # 99.5675771.11111




indel_platypus <- DNA_indel_final[DNA_indel_final$T_PROGS %like% "platypus",]
table(indel_platypus$T_PROGS) #Check - perfect
FP_platypus <- indel_platypus[indel_platypus$mut %in% nuevo_golden_indel_unvalidated$mut,] # 36
TN_platypus <- nuevo_golden_indel_unvalidated[!(nuevo_golden_indel_unvalidated$mut %in% indel_platypus$mut),] # 235
nrow(TN_platypus)
Specificity_metric_platypus_bien <- nrow(TN_platypus) / (nrow(TN_platypus)+nrow(FP_platypus)) *100 # 86.71587
Specificity_metric_platypus_bien
Specificity_metric_platypus <- nrow(TN) / (nrow(TN) + nrow(FP_platypus)) *100
Specificity_metric_platypus # 81.63265
Precision_metric_platypus <- nrow(TP)/(nrow(TP) +nrow(FP_platypus)) *100
Precision_metric_platypus # 99.70022

indel_strelka <- DNA_indel_final[DNA_indel_final$T_PROGS %like% "strelka",]
table(indel_strelka$T_PROGS) #Check - perfect
FP_strelka <- indel_strelka[indel_strelka$mut %in% nuevo_golden_indel_unvalidated$mut,] # 88
TN_strelka <- nuevo_golden_indel_unvalidated[!(nuevo_golden_indel_unvalidated$mut %in% indel_strelka$mut),] # 183
nrow(TN_strelka)
Specificity_metric_strelka_bien <- nrow(TN_strelka) / (nrow(TN_strelka)+nrow(FP_strelka)) *100 # 67.52768
Specificity_metric_strelka_bien
Specificity_metric_strelka <- nrow(TN) / (nrow(TN) + nrow(FP_strelka)) *100
Specificity_metric_strelka # 64.51613
Precision_metric_strelka <- nrow(TP)/(nrow(TP) +nrow(FP_strelka)) *100
Precision_metric_strelka # 99.27038

indel_svaba <- DNA_indel_final[DNA_indel_final$T_PROGS %like% "svaba",]
table(indel_svaba$T_PROGS) #Check - perfect
FP_svaba <- indel_svaba[indel_svaba$mut %in% nuevo_golden_indel_unvalidated$mut,] # 44
TN_svaba <- nuevo_golden_indel_unvalidated[!(nuevo_golden_indel_unvalidated$mut %in% indel_svaba$mut),] # 227
nrow(TN_svaba)
Specificity_metric_svaba_bien <- nrow(TN_svaba) / (nrow(TN_svaba)+nrow(FP_svaba)) *100 # 83.76384
Specificity_metric_svaba_bien
Specificity_metric_svaba <- nrow(TN) / (nrow(TN) + nrow(FP_svaba)) *100
Specificity_metric_svaba # 78.81773
Precision_metric_svaba <- nrow(TP)/(nrow(TP) +nrow(FP_svaba)) *100
Precision_metric_svaba # 99.64214

#ACCURACY
accuracy_1 <- 100*(nrow(TN)+nrow(TP))/(nrow(TN)+nrow(FN)+nrow(TP)+nrow(FP))
accuracy_2 <- 100*(nrow(TN_2)+nrow(golden_in_indel_2_or_more))/(nrow(TN_2)+nrow(FN_2_or_more)+nrow(golden_in_indel_2_or_more)+nrow(FP_2))
accuracy_3 <- 100*(nrow(TN_3)+nrow(golden_in_indel_3_or_more))/(nrow(TN_3)+nrow(FN_3_or_more)+nrow(golden_in_indel_3_or_more)+nrow(FP_3))
accuracy_4 <- 100*(nrow(TN_4)+nrow(golden_in_indel_4_or_more))/(nrow(TN_4)+nrow(FN_4_or_more)+nrow(golden_in_indel_4_or_more)+nrow(FP_4))
accuracy_5 <- 100*(nrow(TN_5)+nrow(golden_in_indel_5_or_more))/(nrow(TN_5)+nrow(FN_5_or_more)+nrow(golden_in_indel_5_or_more)+nrow(FP_5))
accuracy_6 <- 100*(nrow(TN_6)+nrow(golden_in_indel_6))/(nrow(TN_6)+nrow(FN_6)+nrow(golden_in_indel_6)+nrow(FP_6))

accuracy_gatk4m2 <- 100*(nrow(TN_gatk4m2)+nrow(golden_in_indel_gatk4m2))/(nrow(TN_gatk4m2)+nrow(FN_gatk4m2)+nrow(golden_in_indel_gatk4m2)+nrow(FP_gatk4m2))
accuracy_lancet <- 100*(nrow(TN_lancet)+nrow(golden_in_indel_lancet))/(nrow(TN_lancet)+nrow(FN_lancet)+nrow(golden_in_indel_lancet)+nrow(FP_lancet))
accuracy_pindel <- 100*(nrow(TN_pindel)+nrow(golden_in_indel_pindel))/(nrow(TN_pindel)+nrow(FN_pindel)+nrow(golden_in_indel_pindel)+nrow(FP_pindel))
accuracy_platypus <- 100*(nrow(TN_platypus)+nrow(golden_in_indel_platypus))/(nrow(TN_platypus)+nrow(FN_platypus)+nrow(golden_in_indel_platypus)+nrow(FP_platypus))
accuracy_strelka <- 100*(nrow(TN_strelka)+nrow(golden_in_indel_strelka))/(nrow(TN_strelka)+nrow(FN_strelka)+nrow(golden_in_indel_strelka)+nrow(FP_strelka))
accuracy_svaba <- 100*(nrow(TN_svaba)+nrow(golden_in_indel_svaba))/(nrow(TN_svaba)+nrow(FN_svaba)+nrow(golden_in_indel_svaba)+nrow(FP_svaba))

row_names <- c("1 or more","2 or more","3 or more","4 or more","5 or more","6","gatk4m2","lancet","pindel","platypus","strelka","svaba")
specificity_values <- c(Specificity_metric,Specificity_metric_2_bien,Specificity_metric_3_bien, Specificity_metric_4_bien,Specificity_metric_5_bien,Specificity_metric_6_bien,Specificity_metric_gatk4m2_bien,Specificity_metric_lancet_bien,Specificity_metric_pindel_bien,Specificity_metric_platypus_bien,Specificity_metric_strelka_bien,Specificity_metric_svaba_bien)
precision_values <- c(Precision_metric,Precision_metric_2_or_more,Precision_metric_3_or_more, Precision_metric_4_or_more,Precision_metric_5_or_more,Precision_metric_6,Precision_metric_gatk4m2,Precision_metric_lancet,Precision_metric_pindel,Precision_metric_platypus,Precision_metric_strelka,Precision_metric_svaba)

#Nos hemos quedado aqui!!!!!!!!

FP_values <- c(nrow(FP), nrow(FP_2),nrow(FP_3),nrow(FP_4),nrow(FP_5),nrow(FP_6),nrow(FP_gatk4m2),nrow(FP_lancet),nrow(FP_pindel),nrow(FP_platypus),nrow(FP_strelka),nrow(FP_svaba))
TP_values <- c(nrow(TP),nrow(golden_in_indel_2_or_more),nrow(golden_in_indel_3_or_more),nrow(golden_in_indel_4_or_more),nrow(golden_in_indel_5_or_more),nrow(golden_in_indel_6),nrow(golden_in_indel_gatk4m2),nrow(golden_in_indel_lancet),nrow(golden_in_indel_pindel),nrow(golden_in_indel_platypus),nrow(golden_in_indel_strelka),nrow(golden_in_indel_svaba))
FN_values <- c(nrow(FN), nrow(FN_2_or_more),nrow(FN_3_or_more),nrow(FN_4_or_more),nrow(FN_5_or_more),nrow(FN_6),nrow(FN_gatk4m2),nrow(FN_lancet),nrow(FN_pindel),nrow(FN_platypus),nrow(FN_strelka),nrow(FN_svaba))
TN_values <- c(nrow(TN), nrow(TN_2),nrow(TN_3),nrow(TN_4),nrow(TN_5),nrow(TN_6),nrow(TN_gatk4m2),nrow(TN_lancet),nrow(TN_pindel),nrow(TN_platypus),nrow(TN_strelka),nrow(TN_svaba))
accuracy_values <- c(accuracy_1,accuracy_2,accuracy_3,accuracy_4,accuracy_5,accuracy_6, accuracy_gatk4m2,accuracy_lancet,accuracy_pindel,accuracy_platypus,accuracy_strelka,accuracy_svaba)
golden_data <- c("Number of variants", "13107", "Number of reference positions", "271")
METRICS <- data.table(row_names,recall_values,specificity_values,precision_values,FP_values,TP_values,FN_values, TN_values,accuracy_values)
METRICS
golden_data

# "Number of variants"            "13107"   (TP)                      "Number of reference positions" "271" (N
# De estos N, solo podemos usar los que sean FP, ya que solo esos tienen la info necesaria para el modelo (110)
# Sin embargo se pueden usar el resto de TN para las métricas, que estarán presentes en los programas sueltos y NPROGS y en el resultado del RF
# El RF solo tendrá poder para pasar variantes de P a N, con suerte las que pase sean de FP, y con más suerte pase la gran mayoria de FP a N.
# Sin embargo tb puede pasar TP a N y eso baja el recall. Por otra parte 
# Con lo cual, el número de FP jamás podrá aumentar con un RF. Pero si que podrá aumentar el número de TP, (disminuyendo el de FP, )
# 
# row_names     recall_values specificity_values precision_values FP_values TP_values FN_values TN_values accuracy_values
#  1 or more      91.71435           59.40959         99.09323       110     12021      1086       161        91.05995
#  2 or more      83.28374           70.84871         99.34711        79     10916      2191       192        83.03184
#  3 or more      74.81498           79.70480         99.54455        55      9806      3301       216        74.91404
#  4 or more      64.65248           84.50185         99.65183        42      8474      4633       229        65.05457
#  5 or more      49.74441           90.40590         99.78418        26      6520      6587       245        50.56810
#          6      24.22370           98.15498         99.95842         5      3175      9932       266        25.72133
#    gatk4m2      69.76425           76.38376         99.47042        64      9144      3963       207        69.89834
#     lancet      48.58473           87.82288         99.72623        33      6368      6739       238        49.37958
#     pindel      65.10262           80.81181         99.56929        52      8533      4574       219        65.42084
#   platypus      51.16350           86.71587         99.70142        36      6706      6401       235        51.88369
#    strelka      86.11429           67.52768         99.27327        88     11287      1820       183        85.73778
#      svaba      67.70428           83.76384         99.63531        44      8874      4233       227        68.02960

