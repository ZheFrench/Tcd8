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

option_list = list(
  make_option(c("-a", "--annotation"), type="character",  help="File with IDs and Annotations", metavar="PATH2ANNOTATION")
)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = TRUE);
opt  <- arguments$options
args <- arguments$args



base.dir <- "/data/villemin/data/Tcd8/"
nanostring.file   <- "/data/villemin/data/Tcd8/NanoString.normalised.tsv"
#nanostring.file   <- "/data/villemin/data/Tcd8/NanoString.normalised.wtoutliers.tsv"
#nanostring.file   <- "/data/villemin/data/Tcd8/NanoString.raw.32.tsv"


filename<- basename(opt$annotation)

dataframe.Annotation <- fread(opt$annotation,data.table=F)
dataframe.nanostring <- fread(nanostring.file,data.table=F)
#head(dataframe.Annotation)

# Trick to remove value of signal to make low and high group
dataframe.Annotation <- dataframe.Annotation[1:(length(dataframe.Annotation)-1)]

# Round value if needed
dataframe.nanostring[,-1] <-round(dataframe.nanostring[,-1],0)


genes                          <- dataframe.nanostring[,1]
head(genes)
dataframe.nanostring           <- dataframe.nanostring[,-1] # %>% data.matrix( )
rownames(dataframe.nanostring) <- genes
Ids <- unlist(colnames(dataframe.nanostring))

# Transformation if needed in CPM
#y   <- DGEList(counts = dataframe.nanostring)
#cpm <- cpm(y$counts, log = FALSE)
#dataframe.nanostring <- as.data.frame(cpm)

# Loose rownames in transposition
dataframe.nanostring.transposed <- transpose(dataframe.nanostring)

# transpose all but the first column (name)
colnames(dataframe.nanostring.transposed) <- rownames(dataframe.nanostring)
rownames(dataframe.nanostring.transposed) <- Ids
dataframe.nanostring.transposed$Patient    <- Ids

dataframe.nanostring.transposed.annotated <- inner_join(x = dataframe.nanostring.transposed, y = dataframe.Annotation  , by = "Patient")
# From here you get Patient Group ang genes names as column

gene="ITGA6"

# Plot gene value for low and high
png(file=glue("{base.dir}/plots/Boxplot_{gene}.png"),width=400,height=500)
ggplot(dataframe.nanostring.transposed.annotated, aes(factor(group),get(gene), fill=factor(group) )) + geom_boxplot(outlier.shape=NA) +
stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE,label.x = 1.5) + #label.x = 2.45
labs(x = "",fill = "Group :",y = "")  + geom_jitter( position=position_jitter(0.2))+ 
scale_fill_manual(values= c( "low"= "#00BFC4", "high"= "#F8766D"))+ theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum() 
dev.off()

print ("Let's dot it")
processed.dirname <- dirname(opt$annotation)

test <- dataframe.nanostring.transposed.annotated %>% select(-Patient) %>% tbl_summary(by = group) %>% add_p()

nanostring.tested <- data.frame (Name  = test$table_body$variable,p.value = test$table_body$p.value) 
nanostring.tested <- nanostring.tested[!is.na(nanostring.tested$p.value),] 
dim(nanostring.tested)

nanostring.tested$p.adjust <- p.adjust(nanostring.tested$p.value, method = "BH", n = length(nanostring.tested$p.value))
dim(nanostring.tested) #770 3
dim(dataframe.nanostring) #770 35

# Remerged...and filter by p-value
dataframe.nanostring$Name <- rownames(dataframe.nanostring)

dataframe.nanostring.final<- inner_join(x = nanostring.tested  , y = dataframe.nanostring  , by = "Name")
dim(dataframe.nanostring.final) #779 35
rownames(dataframe.nanostring.final) <- dataframe.nanostring.final$Name
# Write to file

dataframe.nanostring.final <- dataframe.nanostring.final  %>% select(Name, everything()) 
write.table(dataframe.nanostring.final,file = glue("{processed.dirname}/{filename}-all-diffsimple.tsv"),quote=F,row.names=F,sep="\t")

dataframe.nanostring.final.filtered.genes <- rownames(dataframe.nanostring.final[dataframe.nanostring.final$p.value < 0.05 ,])
head(dataframe.nanostring.final.filtered.genes)

x <- subset(dataframe.nanostring,rownames(dataframe.nanostring) %in% dataframe.nanostring.final.filtered.genes)

x <- x  %>% select(Name, everything())  #%>% select(-p.value,-p.adjust)

names(x)[1] <- ""
write.table(x,file = glue("{processed.dirname}/{filename}-filtred-diffsimple.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
