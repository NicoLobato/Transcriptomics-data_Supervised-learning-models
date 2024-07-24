###########################
### Data preparation ###
##########################

## Loading libraries

library("data.table")
library("edgeR")
library("tidyverse")

## Reading the data

data <- fread("data/allcountsSalmon.csv", sep = ",")
data <- data[0:10000,]
dim(data) # 206574     26
head(data[,0:10])

genes_names <- data$target_id # guardamos los identificadores de las variables originales
samples_names <- names(data)[3:ncol(data)] # guardamos los de las muestras

data <- as.matrix(data[,-c("V1","target_id")])
rownames(data) <- genes_names
rm(genes_names)

# Antes de comenzar con el análisis, es altamente recomendable realizar un filtrado muy común en transcriptómica que consiste en descartar aquellos genes con un conteo por millón (CPM) inferior a 1 en menos de 3 muestras
# The matrix of read counts is filtered so that I keep genes where cpm>=1 in at least 3 samples (CPM = counts per million)
keep=rowSums(cpm(data)>=1)>=3
length(keep)
data <- data[keep,]
rm(keep)
dim(data) # 44881    24
# Vemos que se reduce bastante el número de genes

data <- t(data) %>%
  as.data.table(.)

# Creamos la estructura de los datos a usar: 1º columna: identificadores de las muestras; 2º columna: tipo de tratamiento (Control o Sequía), tipo de púa y tipo de porta; resto de columnas: variables
treatment_type <- ifelse(endsWith(samples_names, "C"), "Control", "Drought")
stem_type <- ifelse(grepl("PuS", samples_names), "Sensitive", "Tolerant")
root_type <- ifelse(grepl("PoS", samples_names), "Sensitive", "Tolerant")

#type <- paste(treatment_type, stem_type, root_type, sep="_")
#data$treatment_stem_root <- type
prepared_data <- cbind(Sample_ID = samples_names, treatment = treatment_type,
              stem = stem_type, root = root_type, data)

save(prepared_data, file = "data/prepared_data.RData")
rm(list = ls()) # liberamos el espacio para agilizar la ejecución del código
gc()
