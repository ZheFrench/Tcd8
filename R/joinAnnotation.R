rm(list=ls())
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-m", "--matrix"), type="character",  help="A big matrix without annotation", metavar="PATH2MATRIX"),
  make_option(c("-a", "--annotation"), type="character",  help="File with IDs and Annotations", metavar="PATH2ANNOTATION")
)

# clinicData.tsv  
# NanoString.normalised.tsv
# TIS.Score.tsv
# /data/villemin/data/Tcd8/

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = TRUE);
opt  <- arguments$options
args <- arguments$args

print("> OPTS : ")
print("> ARGS : ")
print(args)

# jour/mois/année
base.dir <- "/data/villemin/data/Tcd8/"
dataframe.nanostring <-fread(glue("{base.dir}NanoString.normalised.tsv"),data.table = F)

rownames(dataframe.nanostring) <- dataframe.nanostring$V1
dataframe.nanostring$V1 <- NULL

#dataframe.annotation <-fread(glue("/data/villemin/data/Tcd8/clinicData.tsv"),data.table = F)
#head(dataframe.annotation)

#dataframe.final <-  merge(x = dataframe.matrix, y = dataframe.annotation,  by.x = "genes", by.y = "p_val")
#head(dataframe.final)

dataframe.annotation <-fread(glue("{base.dir}clinicData.tsv"),data.table = F)
#head(dataframe.annotation)
Genes <- unlist(rownames(dataframe.nanostring))

Ids <- unlist(colnames(dataframe.nanostring))

# Loose rownames in transposition
dataframe.nanostring.transposed <- as.data.frame(transpose(dataframe.nanostring))

# transpose all but the first column (name)
colnames(dataframe.nanostring.transposed) <- Genes
rownames(dataframe.nanostring.transposed) <- Ids
dataframe.nanostring.transposed$Id <- Ids

dataframe.nanostring.transposed.annotated <- inner_join(x = dataframe.nanostring.transposed, y = dataframe.annotation %>% select("Id","Best.Response","Responder.1","Responder.2")  , by = "Id")
dataframe.nanostring.transposed.annotated <- dataframe.nanostring.transposed.annotated %>% select(c("Id","Best.Response","Responder.1","Responder.2") , everything())



dataframe.nanostring.annotated <- transpose(dataframe.nanostring.transposed.annotated)

rownames(dataframe.nanostring.annotated) <- c (c("Id","Best.Response","Responder.1","Responder.2" ) , Genes)

write.table(dataframe.nanostring.annotated,file = glue("{base.dir}/heatmap.annotated.tsv"),quote=F,row.names=T,col.names=F,sep="\t")
