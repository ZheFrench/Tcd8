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
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(reshape2))
#suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(gtsummary))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(survminer))
#suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(gt))
#suppressPackageStartupMessages(library(atable))
#suppressPackageStartupMessages(library(tableone))
#suppressPackageStartupMessages(library(table1))
#suppressPackageStartupMessages(library(formattable))
#suppressPackageStartupMessages(library(flextable))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(clusterProfiler))
# http://rstudio-pubs-static.s3.amazonaws.com/250023_a4a1795bc6db4421ad178bd8520b1197.html
# https://www.themillerlab.io/post/survival_analysis/
# https://github.com/LucoLab/Survival/blob/main/SurvivalV4.R
suppressPackageStartupMessages(library(edgeR))

option_list = list(
  make_option(c("-a", "--annotation"), type="character",  help="File with IDs and Annotations", metavar="PATH2ANNOTATION")
)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = TRUE);
opt  <- arguments$options
args <- arguments$args

print("> OPTS : ")
print("> ARGS : ")
print(args)
set.seed(12)

base.dir <- "/data/villemin/data/Tcd8/"

dir.create(glue("{base.dir}/plots"), showWarnings = F)
dir.create(glue("{base.dir}/DE"), showWarnings = F)

dataframe.Annotation <- fread(opt$annotation,data.table = F)
dataframe.tis.score  <- fread(glue("{base.dir}TIS.Score.tsv"),data.table = F)
dataframe.nanostring <- fread(glue("{base.dir}NanoString.normalised.tsv"),data.table = F)
#dataframe.nanostring <- fread(glue("{base.dir}NanoString.raw.32.tsv"),data.table = F)
#dataframe.nanostring <- fread(glue("{base.dir}NanoString.normalised.wtoutliers.tsv"),data.table = F)

head(dataframe.Annotation)

colnames(dataframe.Annotation)[1]<-"Id"
dim(dataframe.Annotation$Patient)
genes <- dataframe.nanostring$V1

dim(dataframe.nanostring)

# Round Dataframe Values
dataframe.nanostring[,-1] <-round(dataframe.nanostring[,-1],0)
head(dataframe.nanostring)


ids.matched <- dataframe.Annotation$Id[dataframe.Annotation$Id %in% colnames(dataframe.nanostring)]
print(dataframe.Annotation$Id)
print(ids.matched)

length(ids.matched)
matrix.nanostring <- dataframe.nanostring %>% select (c(-V1)) %>% select(c(unlist(ids.matched))) %>%  data.matrix( )
dim(matrix.nanostring)
#matrix.nanostring <- dataframe.nanostring %>% select (c(-Names)) %>% data.matrix( )

rownames(matrix.nanostring) <- genes

#head(matrix.nanostring)
#good <- rowSums( apply(counts,1,function(x) x>= 5 ) ) >= 3
#counts <- counts[good,]
#print(glue("# Filtered Genes {length(rownames(counts))}"))

dge <- DGEList(matrix.nanostring) #,genes=rownames(matrix.nanostring)

print("DGE samples.")


#dge <- calcNormFactors(dge)

dge$samples$Id <- rownames(dge$samples)
dim(dge$samples)#31 4 

print ("Annotation")
dim(dataframe.Annotation)#33 3

#atient group CD8Plus.Total.Ex
#1     C43   low                0
#2     C16   low                2
#3     H85   low               23
#4     C90  high              568
#5     C04  high              608
#6     H21  high              613
dge$samples <- dge$samples %>% select (c(-group))
dge$samples  <- inner_join (x = dge$samples ,y = dataframe.Annotation,by ="Id")

dge$samples$group <-  glue("condition_{dge$samples$group}")
print(dge$samples)


cond1 = "high"
cond2 = "low"

print(dge$samples )
# names need to have specific format to go into model matrix meaning see help(names)

# Aucun effet ici???
dge$samples$group <- relevel(factor(dge$samples$group),ref = glue("condition_{cond2}"))
# was each time comparing to high
comp     <- glue("condition_{cond1}-condition_{cond2}")

de.design <- model.matrix(~0 + dge$samples$group)

colnames(de.design) <- gsub("^dge\\$samples\\$group","",colnames(de.design))
print("de.design")
#levels(de.design)
print("dge$samples$group")
#dge$samples$group
#Levels: condition_low condition_high
#                Contrasts
#Levels           condition_high-condition_low
#  condition_low                            -1
#  condition_high                            1
#Levels: condition_low condition_high
#                Contrasts
#Levels           condition_high-condition_low
#  condition_low                            -1
#  condition_high                            1
cm     <- makeContrasts(contrasts = comp,levels = dge$samples$group)
cm
dge    <- estimateDisp(dge, de.design,robust=T)

fit.y  <- glmFit(dge, de.design)
lrt    <- glmLRT(fit.y,contrast = cm)

de.design

#######################################################################################################
###################################      FILES          ###############################################
#######################################################################################################
filename<- basename(opt$annotation)
processed.dirname <- dirname(opt$annotation)

dir.create(glue("{processed.dirname}/DE"), showWarnings = F)

print("Writing differential results...")

result <- as.data.frame( topTags(lrt, adjust.method="BH",n=Inf, sort.by="PValue", p.value=1))
result <- cbind(genes = rownames(result), result, row.names = NULL)

result     <- cbind(result,dge$counts[result$genes,])

write.table(result,file = glue("{processed.dirname}/DE/{filename}_{cond1}_{cond2}-differential.tsv"),quote=F,row.names=F,sep="\t")

result_up    <- subset(result, ( logFC >= 1 & FDR < 0.05 ))
result_up <- result_up[order(abs(result_up$logFC),decreasing = TRUE),]
rownames(result_up)

write.table(result_up,file = glue("{processed.dirname}/DE/{filename}_{cond1}_{cond2}-differential-up.tsv"),quote=F,row.names=F,sep="\t")

result_down <- subset(result, ( logFC <= -1 & FDR < 0.05 ))
result_down <- result_down[order(abs(result_down$logFC),decreasing = TRUE),]
rownames(result_down)

write.table(result_down,file = glue("{processed.dirname}/DE/{filename}_{cond1}_{cond2}-differential-down.tsv"),quote=F,row.names=F,sep="\t")

summary <- as.data.frame(summary(dt<-decideTestsDGE(lrt, adjust.method="BH",p.value = 0.05,lfc = 1)))
colnames(summary)[1] <- "FC"
colnames(summary)[2] <- "Analyse"
write.table(summary ,file = glue("{processed.dirname}/DE/{filename}_{cond1}_{cond2}-summary.tsv"),quote=F,row.names=F,sep="\t")
