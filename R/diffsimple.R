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
suppressPackageStartupMessages(library(gtsummary))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(stats))

#CodeClass,Name,Accession,Count
#Endogenous,CCNO,NM_021147.4,20

base.dir <- "/data/villemin/data/Tcd8/"
nanostring.file   <- "/data/villemin/data/Tcd8/NanoString.raw.tsv"

dataframe.Annotation <- fread(glue("{base.dir}annotation.tsv"),data.table=F)
dataframe.nanostring <- fread(nanostring.file,data.table=F)


# Round value if needed
dataframe.nanostring[,-1] <-round(dataframe.nanostring[,-1],0)


genes                          <- dataframe.nanostring$Name
dataframe.nanostring           <- dataframe.nanostring %>% select (c(-Name)) # %>% data.matrix( )
rownames(dataframe.nanostring) <- genes
Ids <- unlist(colnames(dataframe.nanostring))


# Transformation if needed
y   <- DGEList(counts = dataframe.nanostring)
cpm <- cpm(y$counts, log = FALSE)
dataframe.nanostring <- as.data.frame(cpm)

# Loose rownames in transposition
dataframe.nanostring.transposed <- transpose(dataframe.nanostring)

# transpose all but the first column (name)
colnames(dataframe.nanostring.transposed) <- rownames(dataframe.nanostring)
rownames(dataframe.nanostring.transposed) <- Ids
dataframe.nanostring.transposed$Id        <- Ids

dataframe.nanostring.transposed.annotated <- inner_join(x = dataframe.nanostring.transposed, y = dataframe.Annotation  , by = "Id")
gene="ITGA6"
png(file=glue("{base.dir}/plots/Boxplot_{gene}.png"),width=400,height=500)
ggplot(dataframe.nanostring.transposed.annotated, aes(factor(Group),get(gene), fill=factor(Group) )) + geom_boxplot(outlier.shape=NA) +
stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE,label.x = 1.5) + #label.x = 2.45
labs(x = "",fill = "Group :",y = "")  + geom_jitter( position=position_jitter(0.2))+ 
scale_fill_manual(values= c( "1"= "#00BFC4", "2"= "#F8766D"))+ theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum() 
dev.off()

print ("Let's dot it")


test <- dataframe.nanostring.transposed.annotated %>% select(-Id) %>% tbl_summary(by = Group) %>% add_p()

nanostring.tested <- data.frame (Name  = test$table_body$variable,p.value = test$table_body$p.value) 


nanostring.tested$p.adjust <- p.adjust(nanostring.tested$p.value, method = "BH", n = length(nanostring.tested$p.value))

# Remerged...and filter by p-value
dataframe.nanostring$Name <- rownames(dataframe.nanostring)
dataframe.nanostring.final<- left_join(x = dataframe.nanostring, y = nanostring.tested  , by = "Name")

rownames(dataframe.nanostring.final) <- dataframe.nanostring$Name 

dataframe.nanostring.final <- dataframe.nanostring.final %>% select(-Name)

# Write to file
write.table(dataframe.nanostring.final,file = glue("{base.dir}/all-diffsimple.tsv"),quote=F,row.names=T,sep="\t")

dataframe.nanostring.final.filtered.genes <- rownames(dataframe.nanostring.final[!is.na(dataframe.nanostring.final$p.value) & dataframe.nanostring.final$p.adjust < 0.05 ,])

head(dataframe.nanostring.final.filtered.genes)

x <- subset(dataframe.nanostring, rownames(dataframe.nanostring) %in% dataframe.nanostring.final.filtered.genes)

write.table(x,file = glue("{base.dir}/filtred-diffsimple.tsv"),quote=F,row.names=T,col.names=T,sep="\t")
