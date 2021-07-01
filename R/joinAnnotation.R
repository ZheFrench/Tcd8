rm(list=ls())
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))


base.dir <- "/data/villemin/data/Tcd8/"


option_list = list(
  make_option(c("-a", "--annotation"), type="character",  help="File with IDs and Annotations", metavar="PATH2ANNOTATION"),
  make_option(c("-e", "--expression"), type="character",  default = glue("{base.dir}NanoString.normalised.tsv"),help="File with IDs and Annotations", metavar="PATH2ANNOTATION")

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

print(opt$expression)
dataframe.clinical <-fread(glue("{base.dir}clinicData.tsv"),data.table = F)
dataframe.nanostring <-fread(opt$expression,data.table = F)
dataframe.tis.score  <- fread(glue("{base.dir}TIS.Score.tsv"),data.table=F)

rownames(dataframe.nanostring) <- dataframe.nanostring$V1
dataframe.nanostring$V1 <- NULL

dataframe.annotation <-fread(opt$annotation,data.table = F)

filname.annotation <- basename(opt$annotation)

dataframe.final <-  merge(x = dataframe.annotation, y = dataframe.clinical,  by.x = "Patient", by.y = "Id")
dim(dataframe.final)


Genes <- unlist(rownames(dataframe.nanostring))

Ids <- unlist(colnames(dataframe.nanostring))

# Loose rownames in transposition
dataframe.nanostring.transposed <- as.data.frame(transpose(dataframe.nanostring))

# transpose all but the first column (name)
colnames(dataframe.nanostring.transposed) <- Genes
rownames(dataframe.nanostring.transposed) <- Ids
dataframe.nanostring.transposed$Patient <- Ids
# Patient is a column. All other columns are genes.

# Subset what you want to keep
dataframe.final <- dataframe.final %>% select("Patient","Best.Response","Responder.1","Responder.2","Histo","Tabaco" ,"Nb.Pa", "Sexe","Age.First.Immunotherapy","group") 
# We merge dataframe nanostring.transposed with dataframe.final (clinical data+ Group high low)
dataframe.nanostring.transposed.annotated <- inner_join(x = dataframe.nanostring.transposed, y =dataframe.final ,  by= "Patient")
dim(dataframe.nanostring.transposed.annotated)

dataframe.nanostring.transposed.annotated <- dataframe.nanostring.transposed.annotated %>% select(c("Patient","Best.Response","Responder.1","Responder.2","Histo","Tabaco" ,"Nb.Pa", "Sexe","Age.First.Immunotherapy","group") , everything())

dataframe.nanostring.annotated <- transpose(dataframe.nanostring.transposed.annotated)

rownames(dataframe.nanostring.annotated) <- c (c("Patient","Best.Response","Responder.1","Responder.2","Histo","Tabaco" ,"Nb.Pa", "Sexe","Age.First.Immunotherapy","group" ) , Genes)
dim(dataframe.nanostring.annotated)

# Change base.dir to be in same directory of file processessed
processed.dirname <- dirname(opt$annotation)
# don't forget you have 62 patients with clinical data, 32 expression...
write.table(dataframe.nanostring.annotated,file = glue("{processed.dirname}/heatmap.{filname.annotation}.annotated.tsv"),quote=F,row.names=T,col.names=F,sep="\t")

png(file=glue("{processed.dirname}/{filname.annotation}.Boxplot_Age.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$group),], aes(factor(group),Age.First.Immunotherapy, fill=factor(group))) + 
stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE,label.x = 1.5) + geom_boxplot(outlier.shape=NA) + #label.x = 2.45
labs(x = "",fill = "Group :",y = "")  + geom_jitter( position=position_jitter(0.2))+ 
scale_fill_manual(values= c( "low"= "#00BFC4", "high"= "#F8766D"))+ theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum() 
dev.off()

png(file=glue("{processed.dirname}/{filname.annotation}.Barplot_Best.Response.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$group),] )+ geom_bar( aes( factor(group),fill = factor(Best.Response)) ,color="black")  + 
labs(x = "Group",fill = "Best.Response :",y = "") + theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum()  + scale_fill_brewer(palette="Blues")
dev.off()

png(file=glue("{processed.dirname}/{filname.annotation}.Boxplot_Nb.Pa.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$group),], aes(factor(group),Nb.Pa, fill=factor(group))) + geom_boxplot() + 
stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE,label.x = 1.5) + geom_boxplot(outlier.shape=NA) + #label.x = 2.45
labs(x = "",fill = "Group :",y = "")  + geom_jitter( position=position_jitter(0.2))+ 
scale_fill_manual(values= c( "low"= "#00BFC4", "high"= "#F8766D"))+ theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum() 
dev.off()

png(file=glue("{processed.dirname}/{filname.annotation}.Barplot_Histo.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$group),] )+ geom_bar( aes( factor(group),fill = factor(Histo)),color="black"   )   + 
labs(x = "Group",fill = "Histo :",y = "") + theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum()  + scale_fill_brewer(palette="Blues")
dev.off()

png(file=glue("{processed.dirname}/{filname.annotation}.Barplot_Sexe.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$group),] )+ geom_bar( aes( factor(group) ,fill = factor(Sexe)) ,color="black" )  + 
labs(x = "Group",fill = "Sexe :",y = "") + theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum()  + scale_fill_brewer(palette="Blues")
dev.off()

png(file=glue("{processed.dirname}/{filname.annotation}.Barplot_Tabaco.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$group),] )+ geom_bar( aes( factor(group),fill = factor(Tabaco)) ,color="black")  + 
labs(x = "Group",fill = "Tabaco :",y = "") + theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum()  + scale_fill_brewer(palette="Blues")
dev.off()

names(dataframe.tis.score)[1] <- "Patient"
dataframe.final <-  inner_join(x = dataframe.final , y = dataframe.tis.score ,  by.x = "Patient",by.y = "Id")
dim(dataframe.final)
png(file=glue("{processed.dirname}/{filname.annotation}.Boxplot_TIS_Score.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$group),], aes(factor(group),TIS_Score, fill=factor(group))) +
stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE,label.x = 1.5) + geom_boxplot(outlier.shape=NA) + #label.x = 2.45
labs(x = "",fill = "Group :",y = "")  + geom_jitter( position=position_jitter(0.2))+ 
scale_fill_manual(values= c( "low"= "#00BFC4", "high"= "#F8766D"))+ theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum() 
dev.off()
