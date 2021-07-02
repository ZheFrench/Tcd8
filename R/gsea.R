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
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(gtsummary))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(clusterProfiler))


option_list = list(
  make_option(c("-d", "--differential"), type="character",  help="File with IDs and Annotations", metavar="PATH2ANNOTATION"),
  make_option(c("-a", "--annotation"), type="character",  help="Annotation to test", metavar="PATH2ANNOTATION")

)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = TRUE);
opt  <- arguments$options
args <- arguments$args


dirname    <- dirname(opt$differential)
filename   <- basename(opt$differential)
annotation <- opt$annotation

base.dir <- dirname

dir.create(glue("{base.dir}/fgsea"), showWarnings = F)

asbolutepath2file <- opt$differential


if (annotation == "C5"){file.gmt <- "/data/villemin/annotation/gsea/MSigDB/c5.go.v7.2.symbols.gmt" }
if (annotation == "C2"){file.gmt <- "/data/villemin/annotation/gsea/MSigDB/c2.all.v7.2.symbols.gmt" }
if (annotation == "C3"){file.gmt <- "/data/villemin/annotation/gsea/MSigDB/c3.tft.gtrd.v7.2.symbols.gmt" }
if (annotation == "C6"){file.gmt <- "/data/villemin/annotation/gsea/MSigDB/c6.all.v7.2.symbols.gmt" }
if (annotation == "H"){file.gmt  <- "/data/villemin/annotation/gsea/MSigDB/h.all.v7.2.symbols.gmt" }
if (annotation == "C7"){file.gmt <- "/data/villemin/annotation/gsea/MSigDB/c7.all.v7.2.symbols.gmt" }
if (annotation == "C8"){file.gmt <- "/data/villemin/annotation/gsea/MSigDB/c8.all.v7.2.symbols.gmt" }


base <-sub('\\..[^\\.]*$', '', basename(file.gmt) )
print(base)

h.All <- gmtPathways(file.gmt) # 6226
h.All.bis <- read.gmt(file.gmt)

final.file.padj <- glue("{base.dir}/fgsea/{filename}-{annotation}-fgsea.padj.txt")
final.file.nes  <- glue("{base.dir}/fgsea/{filename}-{annotation}-fgsea.nes.txt")

dataframe.expression <- fread(asbolutepath2file,data.table=F)
dataframe.expression <- subset(dataframe.expression,select=c(genes,logFC))

dataframe.expression <- dataframe.expression[order(dataframe.expression$logFC),]

# ClusterProfiler need a decreasing order...
dataframe.expression.decreasing <- dataframe.expression[order(dataframe.expression$logFC, decreasing = TRUE),]
ranks_decreasing <- deframe(dataframe.expression.decreasing)

ranks <- deframe(dataframe.expression)

#not using fgseaMultilevel...a really tiny difference in NES score due to the fact it use fgsea.10/325 5 /250
# But genes used are the same so I used it to plot with heatplot function of clusterProfiler the leading edge genes contained in the signature I am interested in
egmt2 <- GSEA(ranks_decreasing, TERM2GENE = h.All.bis, by = "fgsea",verbose=TRUE ,nPermSimple = 10000 ,minGSSize  = 10, maxGSSize  = 325 , eps = 0,  pvalueCutoff = 1)

# You dont need the whole object to be written
write.table(egmt2, file=glue("{base.dir}/fgsea/{filename}-{annotation}-full-gsea-clusterprofiler.tsv"),quote=F,row.names=F,sep="\t")

# Yeah I do it again (I know.Don't say a fucking word moron.)
fgseaRes     <- fgseaMultilevel(pathways=h.All, stats=ranks,eps=0, nPermSimple = 10000 ,minSize  = 10, maxSize  = 325)

for (pathway in names(h.All)){
    if (pathway %in% c("HALLMARK_INFLAMMATORY_RESPONSE" )){
      #print(pathway)
      #print(h.All[[pathway]])
      p<- plotEnrichment(h.All[[pathway]], ranks) + labs(title=pathway)

      png(file=glue("{base.dir}/fgsea/{filename}-{annotation}-{pathway}.png"))
      print(p)
      dev.off()

    }  
    
    
    if (pathway %in% c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" )){
      #print(pathway)
      #print(h.All[[pathway]])
      p<- plotEnrichment(h.All[[pathway]], ranks) + labs(title=pathway)

      png(file=glue("{base.dir}/fgsea/{filename}-{annotation}-{pathway}.png"))
      print(p)
      dev.off()

    }  
}


fgseaResTidy <- fgseaRes %>% as_tibble()

fwrite( filter(fgseaResTidy ,padj <= 1), file=glue("{base.dir}/fgsea/{filename}-{annotation}-full-fgsea.tsv"), sep="\t", sep2=c("", "/", ""))

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]

topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]


topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

png(file=glue("{base.dir}/fgsea/{filename}-{annotation}-top-global.png"),width=900)
plotGseaTable(h.All[topPathways], ranks, fgseaRes, gseaParam=0.5) 
dev.off()


dataset.padj <-  fgseaResTidy %>% select(-leadingEdge, -pval,-log2err,-size, -ES,-NES) # %>% rename( !!file := padj)
dataset.nes  <-  fgseaResTidy %>% select(-leadingEdge, -pval,-log2err,-size, -ES,-padj) #%>% rename(!!file :=  NES)
    
write.table(dataset.padj,file=final.file.padj,quote=F,row.names=F,sep="\t")
write.table(dataset.nes,file=final.file.nes,quote=F,row.names=F,sep="\t")
