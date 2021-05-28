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
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(viridis))

git config --global https.proxy https://proxy.company.com:8888

base.dir <- "/data/villemin/data/Tcd8/"

dir.create(glue("{base.dir}/plots"), showWarnings = F)

tis.score <- fread(glue("{base.dir}TIS_Score.tsv"),data.table=F)
head(tis.score)

density.plot <- tis.score %>% 
  ggplot(aes(x = TIS_Score, fill = "grey")) + 
   geom_density(alpha = 0.2) + 
  theme_ipsum() +
  ylab("TIS Score Density") 


png(file=glue("{base.dir}/plots/TIS_density.png"),width=900,height=800)
density.plot
dev.off()

#write.table(dataset.nes,file=final.file.nes,quote=F,row.names=F,sep="\t")
