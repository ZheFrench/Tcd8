#################################################################
#
# date: June 01, 2021
# platform: Ubuntu 10.04
# R.version : 4.0.3
# author: Villemin Jean-Philippeall_of
# Institute : IRCM
#
# description.R
# Usage : all_of
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

# http://rstudio-pubs-static.s3.amazonaws.com/250023_a4a1795bc6db4421ad178bd8520b1197.html
# https://www.themillerlab.io/post/survival_analysis/
# https://github.com/LucoLab/Survival/blob/main/SurvivalV4.R
base.dir <- "/data/villemin/data/Tcd8/"

option_list = list(
  make_option(c("-f", "--file"), type="character",  help="Path to annotated file", metavar="PATH2MATRIX"),
  make_option(c("-v", "--variable"), type="character",  help="Variable to check", metavar="VARIABLE"),
  make_option(c("-a", "--analysis"), type="character",  help="Analyse OS or PFS", metavar="ANALYSE"),
  make_option(c("-c", "--cutoff"), type="character",  help="Median or Tertile", metavar="CUTOFF")

)

# clinicData.tsv  
# NanoString.normalised.tsv
# TIS.Score.tsv
# /data/villemin/data/Tcd8/ TIS.Score.tsv

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = TRUE);
opt  <- arguments$options
args <- arguments$args

print("===> SURVIVAL SCRIPT: ")

print("> OPTS : ")
print(opt)


dataframe.Annotation <- fread(opt$file,data.table=F)
variable <- opt$variable
dataframe.clinical   <- fread(glue("{base.dir}clinicData.tsv"),data.table=F)
dataframe.tis.score  <- fread(glue("{base.dir}TIS.Score.tsv"),data.table=F)

dim(dataframe.Annotation)
dataframe.clinical.selected <- dataframe.clinical %>% select(Sexe,Histo,Best.Response,Tabaco,Nb.Pa)

tbl_summary(dataframe.clinical.selected) %>%  bold_labels() %>% as_gt() %>%  gtsave( "stats.html", inline_css = TRUE, path = glue("{base.dir}/")) 


dataframe.clinical$Date.Death.Or.Last.Contact <- mdy(dataframe.clinical$Date.Death.Or.Last.Contact)
dataframe.clinical$Date.First.Immunotherapy   <- mdy(dataframe.clinical$Date.First.Immunotherapy)
dataframe.clinical$Date.Biopsy                <- mdy(dataframe.clinical$Date.Biopsy)

names(dataframe.Annotation)[1] <-"Id"
names(dataframe.Annotation)[2] <-"Group"

#head(dataframe.Annotation)

# C09 DOUBLE
dataframe.final <-  left_join(x = dataframe.Annotation , y = dataframe.clinical ,  by = "Id")
#dataframe.final <-  left_join(x = dataframe.final , y = dataframe.tis.score ,  by = "Id")
print(dataframe.final)
dataframe.final <- mutate(dataframe.final,Death.binary=ifelse(Death.binary=="oui", 1,0))

dataframe.final <- dataframe.final[!is.na(dataframe.final$Group), ]
dataframe.final <- dataframe.final[!is.na(dataframe.final$Date.Biopsy), ]
#dataframe.final <- dataframe.final[dataframe.final$Histo=="ADK", ]
if (opt$analysis == "PFS"){
 dataframe.final <- dataframe.final[!is.na(dataframe.final$Date.Progression.or.Endpoint), ] 
}

dataframe.final$Date.Progression.or.Endpoint   <- mdy(dataframe.final$Date.Progression.or.Endpoint)

X.months <- 18 
type.survival.analysis = opt$analysis

#############################################
# OS -     # Limit to 6month , 12 months etc
#############################################
if(FALSE){
    dataframe.final$temp <-  dataframe.final$Date.First.Immunotherapy %m+% months(X.months) 

    dataframe.final[dataframe.final$Date.Death.Or.Last.Contact > dataframe.final$temp & dataframe.final$Death.binary == 1, ]$Death.binary <- 0

    dataframe.final[dataframe.final$Date.Death.Or.Last.Contact > dataframe.final$temp , ]$Date.Death.Or.Last.Contact <- dataframe.final[dataframe.final$Date.Death.Or.Last.Contact >  dataframe.final$temp , ]$temp 

    dataframe.final$temp  <- NULL
}

dataframe.final$Diff.days   <- interval(dataframe.final$Date.First.Immunotherapy , dataframe.final$Date.Death.Or.Last.Contact) %>% as.numeric('days')
dataframe.final$Diff.months <- interval(dataframe.final$Date.First.Immunotherapy , dataframe.final$Date.Death.Or.Last.Contact)  %>% as.numeric('months')
dataframe.final$Diff.years  <- interval(dataframe.final$Date.First.Immunotherapy , dataframe.final$Date.Death.Or.Last.Contact) %>% as.numeric('years')

#########################
######## PFS  ##########
########################

if( opt$analysis == "PFS" ){
    
    type.survival.analysis = opt$analysis
    # Si date de progression est antérieure à la date d'immunothérapie ya un prob C44 / C90 -> la ou ya 2 dates l'une avt therapie et l'autre 13 pr le mois
    dataframe.final <- dataframe.final[dataframe.final$Date.Progression.or.Endpoint > dataframe.final$Date.First.Immunotherapy,]
    # Limit to 6month , 12 months etc
    ###################################
    if(FALSE){
        dataframe.final$temp <-  dataframe.final$Date.First.Immunotherapy %m+% months(X.months) 

        dataframe.final[dataframe.final$Date.Progression.or.Endpoint > dataframe.final$temp & dataframe.final$Death.binary == 1, ]$Death.binary <- 0
        dataframe.final[dataframe.final$Date.Progression.or.Endpoint >  dataframe.final$temp , ]$Date.Progression.or.Endpoint <- dataframe.final[dataframe.final$Date.Progression.or.Endpoint >  dataframe.final$temp , ]$temp 

        dataframe.final$temp  <- NULL
    }

    # PFS
    dataframe.final$Diff.days   <- interval(dataframe.final$Date.First.Immunotherapy , dataframe.final$Date.Progression.or.Endpoint) %>% as.numeric('days')
    
    print(dataframe.final)
    dataframe.final[dataframe.final$Date.Death.Or.Last.Contact !=  dataframe.final$Date.Progression.or.Endpoint,]$Death.binary <- 1
}

dataframe.final

dataframe.final
dim(dataframe.final)

#dataframe.final$Group <- dataframe.final$Responder.2
#dataframe.final <- mutate(dataframe.final,Group=ifelse(Responder.2=="Oui", 1,2))


surv_object <- Surv(time = dataframe.final$Diff.days, event = dataframe.final$Death.binary)

fit1        <- survfit(surv_object ~ Group, data = dataframe.final)

test = surv_pvalue(fit1, dataframe.final)


fit.coxph          <- coxph(surv_object ~ Group, data = dataframe.final)
confidenceInterval <- summary(fit.coxph)$conf.int

print(test)
print(test$pval)
print(summary(fit.coxph))

pvc             <- coef(summary(fit.coxph))[,5]
hr              <- coef(summary(fit.coxph))[,2]
low_confidence  <-confidenceInterval[,3]
high_confidence <-confidenceInterval[,4]

print(low_confidence)	
print(high_confidence)
print(hr)
print(pvc)

tobepasted = c("HR = ",round(hr,digits = 2),"[",round(low_confidence,digits = 2),"-",round(high_confidence,digits = 2),"]")
title <- paste(tobepasted, collapse = "",sep = "")
 
mypalette = c("#F8766D" , "#00BFC4" ) # high -> red alphabetic orders
#values= c( "high"= "#00BFC4", "low"= "#F8766D")
# "high"= "#00BFC4", "low"= "#F8766D"
# Année diagnostic ??
file      = glue("{base.dir}/plots/survival/{opt$cutoff}_{opt$analysis}/{variable}_Surv.png")
table(dataframe.final$Group)
# use ggpar to change the font of one or a list of ggplots at once
group1 <- sprintf(glue("Group High (n=%s)"), table(dataframe.final$Group)[1])
group2 <- sprintf(glue("Group Low (n=%s)"), table(dataframe.final$Group)[2])

table(dataframe.final$Group)
    write.table(dataframe.final,file = glue("/data/villemin/data/Tcd8/plots/survival/{opt$cutoff}_{opt$analysis}/{variable}_survival.csv"), row.names =  FALSE, quote = F,sep = "\t")

dataframe.final$Group

if (X.months< 12){
ggsurv <- 	ggsurvplot(fit1, data = dataframe.final, pval = TRUE ,   legend.labs = c(group1, group2), legend.title  = "Group :",title = title ,xlab="Months",palette = mypalette) 
} else
{ ggsurv <- 	ggsurvplot(fit1, data = dataframe.final, pval = TRUE ,   legend.labs = c(group1, group2), legend.title  = "Group :",title = title ,xlab="Years",palette = mypalette,break.x.by = 365.25,xscale= 365.25)
}

legend <- glue("{type.survival.analysis}: {variable}")
ggsurv$plot <- ggsurv$plot +  ggplot2::annotate("text", x = 500, y = 1,label = legend, size = 2)

png(file,width = 12,height = 10,units = 'cm', res = 300)

print({	ggpubr::ggpar(ggsurv, 
                    legend = "top",
                    font.legend = list(size = 12, color = "black", face = "bold"),
                    font.x = c(12 , "bold", "grey"),           # font for x axises in the plot, the table and censor part
                    font.y = c(12, "bold", "grey"),       # font for y axises in the plot, the table and censor part
                    font.tickslab = c(12, "plain", "black") ,  
        
                    ggtheme = theme_ipsum()  
)   })
dev.off()
#00BFC4 -> bleu
#F8766D -> rouge



png(file=glue("{base.dir}/plots/survival/{opt$cutoff}_{opt$analysis}/{variable}_Boxplot_Age.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),], aes(factor(Group),Age.First.Immunotherapy, fill=factor(Group))) + 
stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE,label.x = 1.5) + geom_boxplot(outlier.shape=NA) + #label.x = 2.45
labs(x = "",fill = "Group :",y = "")  + geom_jitter( position=position_jitter(0.2))+ 
scale_fill_manual(values= c("high"= "#F8766D", "low"= "#00BFC4"))+ theme(legend.position="right",axis.text.y = element_text(size=14)) + theme_ipsum() 
dev.off()

#png(file=glue("{base.dir}/plots/survival/{variable}_Boxplot_TIS_Score.png"),width=400,height=500)
#ggplot(dataframe.final[!is.na(dataframe.final$Group),], aes(factor(Group),TIS_Score, fill=factor(Group))) +
#stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE,label.x = 1.5) + geom_boxplot(outlier.shape=NA) + #label.x = 2.45
#labs(x = "",fill = "Group :",y = "")  + geom_jitter( position=position_jitter(0.2))+ 
#scale_fill_manual(values= c( "high"= "#00BFC4", "low"= "#F8766D"))+ theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum() 
#dev.off()

png(file=glue("{base.dir}/plots/survival/{opt$cutoff}_{opt$analysis}/Boxplot_Nb.Pa.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),], aes(factor(Group),Nb.Pa, fill=factor(Group))) + geom_boxplot() + 
stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE,label.x = 1.5) + geom_boxplot(outlier.shape=NA) + #label.x = 2.45
labs(x = "",fill = "Group :",y = "")  + geom_jitter( position=position_jitter(0.2))+ 
scale_fill_manual(values= c("high"= "#F8766D", "low"= "#00BFC4"))+ theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum() 
dev.off()

png(file=glue("{base.dir}/plots/survival/{opt$cutoff}_{opt$analysis}/Barplot_Histo.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),] )+ geom_bar( aes( factor(Group),fill = factor(Histo)),color="black"   )   + 
labs(x = "Group",fill = "Histo :",y = "") + theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum()  + scale_fill_brewer(palette="Blues")
dev.off()

png(file=glue("{base.dir}/plots/survival/{opt$cutoff}_{opt$analysis}/Barplot_Sexe.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),] )+ geom_bar( aes( factor(Group) ,fill = factor(Sexe)) ,color="black" )  + 
labs(x = "Group",fill = "Sexe :",y = "") + theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum()  + scale_fill_brewer(palette="Blues")
dev.off()

png(file=glue("{base.dir}/plots/survival/{opt$cutoff}_{opt$analysis}/Barplot_Best.Response.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),] )+ geom_bar( aes( factor(Group),fill = factor(Best.Response)) ,color="black")  + 
labs(x = "Group",fill = "Best.Response :",y = "") + theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum()  + scale_fill_brewer(palette="Blues")
dev.off()

png(file=glue("{base.dir}/plots/survival/{opt$cutoff}_{opt$analysis}/Barplot_Tabaco.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),] )+ geom_bar( aes( factor(Group),fill = factor(Tabaco)) ,color="black")  + 
labs(x = "Group",fill = "Tabaco :",y = "") + theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum()  + scale_fill_brewer(palette="Blues")
dev.off()


df <- data.frame(NAME=character(),PVAL=double(),HR=double(),PVALHR=double(),CILOW=double(),CIHIGH=double())

df    <- rbind(df, data.frame(NAME = variable, PVAL = test$pval,HR=hr,PVALHR=pvc,CILOW=low_confidence,CIHIGH=high_confidence))

if ( file.exists(glue("{base.dir}/plots/survival/{opt$cutoff}_{opt$analysis}/wholeSurv_{opt$cutoff}_{opt$analysis}.csv"))) {
       write.table(df,file=glue("{base.dir}/plots/survival/{opt$cutoff}_{opt$analysis}/wholeSurv_{opt$cutoff}_{opt$analysis}.csv"),quote=FALSE,row.names=FALSE,col.names = TRUE, append = TRUE)
} else {
    
        write.table(df,file=glue("{base.dir}/plots/survival/{opt$cutoff}_{opt$analysis}/wholeSurv_{opt$cutoff}_{opt$analysis}.csv"),quote=FALSE,row.names=FALSE,col.names = FALSE, append = TRUE)
    }





stop()#, label.y = 0.95,

head(dataframe.tis.score)
head(dataframe.clinical)
head(dataframe.nanostring)

density.plot <- dataframe.tis.score %>% 
  ggplot(aes(x = TIS_Score, fill = "grey")) + 
   geom_density(alpha = 0.2) + 
  theme_ipsum() +
  ylab("TIS Score Density") 

png(file=glue("{base.dir}/plots/TIS_density.png"),width=900,height=800)
density.plot
dev.off()
