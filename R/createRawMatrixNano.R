#################################################################
#
# date: June 01, 2021
# platform: Ubuntu 10.04
# R.version : 4.0.3
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# description.R
# Usage : 
# 
# description.R
# 
# Description : 
#
# 
#
#################################################################

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
#

#CodeClass,Name,Accession,Count
#Endogenous,CCNO,NM_021147.4,20

base.dir <- "/data/villemin/data/Tcd8/raw_nanostring/CLEAN/"
final.file <- "/data/villemin/data/Tcd8/NanoString.raw.tsv"

files.list <- list.files(glue("{base.dir}"),pattern="(*).clean.csv$")

for (file in files.list){
  
  print(file)
  asbolutepath2file <- glue("{base.dir}/{file}")
  print(asbolutepath2file)
  
  if (!exists("dataset") ){
    dataset <- fread(asbolutepath2file,data.table=F)
    ID <- str_extract(file,pattern="_C\\d+_|_H\\d+_")
    ID <- str_sub(ID, 2, -2)
    print(ID)# Positive Negative Housekeeping
    dataset <- dataset %>% filter(CodeClass =="Endogenous") %>%  select(Name,Count) %>% rename( !!ID := Count )
    print(head(dataset))
   
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <- fread(asbolutepath2file,data.table=F)
    ID <- str_extract(file,pattern="_C\\d+_|_H\\d+_")
    ID <- str_sub(ID, 2, -2)
    print(ID)
    temp_dataset <- temp_dataset %>% filter(CodeClass =="Endogenous") %>%  select(Name,Count) %>% rename( !!ID := Count )
    dataset <- inner_join(dataset,temp_dataset, by =c("Name"),keep=FALSE)

    rm(temp_dataset)
  }
  
}


clean.dataset <- dataset[grep(".y$", names(dataset), invert = TRUE)]
colnames(clean.dataset) <- gsub(".x",'', colnames(clean.dataset), fixed=TRUE)
head(clean.dataset)

write.table(clean.dataset,file=final.file,quote=F,row.names=F,sep="\t")

