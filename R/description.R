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

# http://rstudio-pubs-static.s3.amazonaws.com/250023_a4a1795bc6db4421ad178bd8520b1197.html
# https://www.themillerlab.io/post/survival_analysis/
# https://github.com/LucoLab/Survival/blob/main/SurvivalV4.R
base.dir <- "/data/villemin/data/Tcd8/"

dir.create(glue("{base.dir}/plots"), showWarnings = F)

dataframe.Annotation <- fread(glue("{base.dir}annotation.tsv"),data.table=F)

dataframe.tis.score  <- fread(glue("{base.dir}TIS.Score.tsv"),data.table=F)
dataframe.clinical   <- fread(glue("{base.dir}clinicData.tsv"),data.table=F)
dataframe.nanostring <- fread(glue("{base.dir}NanoString.normalised.tsv"),data.table=F)


dataframe.clinical.selected <- dataframe.clinical %>% select(Sexe,Histo,Best.Response,Tabaco,Nb.Pa)

tbl_summary(dataframe.clinical.selected) %>%  bold_labels() %>% as_gt() %>%  gtsave( "stats.html", inline_css = TRUE, path = glue("{base.dir}/")) 


dataframe.clinical$Date.Death.Or.Last.Contact <- mdy(dataframe.clinical$Date.Death.Or.Last.Contact)
dataframe.clinical$Date.First.Immunotherapy   <- mdy(dataframe.clinical$Date.First.Immunotherapy)
dataframe.clinical$Date.Biopsy                <- mdy(dataframe.clinical$Date.Biopsy)


# C09 DOUBLE
dataframe.final <-  left_join(x = dataframe.Annotation , y = dataframe.clinical ,  by = "Id")
dataframe.final <-  left_join(x = dataframe.final , y = dataframe.tis.score ,  by = "Id")

dataframe.final <- mutate(dataframe.final,Death.binary=ifelse(Death.binary=="oui", 1,0))

dataframe.final <- dataframe.final[!is.na(dataframe.final$Group), ]
dataframe.final <- dataframe.final[!is.na(dataframe.final$Date.Biopsy), ]
#dataframe.final <- dataframe.final[dataframe.final$Histo=="ADK", ]

dataframe.final <- dataframe.final[!is.na(dataframe.final$Date.Progression.or.Endpoint), ]
dataframe.final$Date.Progression.or.Endpoint   <- mdy(dataframe.final$Date.Progression.or.Endpoint)

dataframe.final

X.months <- 18 

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

if(TRUE){
    #
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
    dataframe.final[dataframe.final$Date.Death.Or.Last.Contact !=  dataframe.final$Date.Progression.or.Endpoint,]$Death.binary <- 1
}

dataframe.final
dataframe.final$Group <- dataframe.final$Responder.2
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

mypalette = c("#00BFC4", "#F8766D")
# Année diagnostic ??
file      = glue("{base.dir}/plots/Surv.png")

# use ggpar to change the font of one or a list of ggplots at once
group1 <- sprintf(glue("Group 1 (n=%s)"), table(dataframe.final$Group)[1])
group2 <- sprintf(glue("Group 2 (n=%s)"), table(dataframe.final$Group)[2])

table(dataframe.final$Group)

dataframe.final$Group

if (X.months< 12){
ggsurv <- 	ggsurvplot(fit1, data = dataframe.final, pval = TRUE ,   legend.labs = c(group1, group2), legend.title  = "Group :",title = title ,xlab="Months",palette = mypalette)  } else
{ ggsurv <- 	ggsurvplot(fit1, data = dataframe.final, pval = TRUE ,   legend.labs = c(group1, group2), legend.title  = "Group :",title = title ,xlab="Years",palette = mypalette,break.x.by = 365.25,xscale= 365.25)}

png(file,width = 10,height = 10,units = 'cm', res = 300)


print({	ggpubr::ggpar(ggsurv, 
                    legend = "top",
                    font.legend = list(size = 12, color = "black", face = "bold"),
                    font.x = c(12 , "bold", "grey"),           # font for x axises in the plot, the table and censor part
                    font.y = c(12, "bold", "grey"),       # font for y axises in the plot, the table and censor part
                    font.tickslab = c(12, "plain", "black") ,  
                 
                    ggtheme = theme_ipsum() 
) })
dev.off()
#00BFC4 -> bleu
#F8766D -> rouge
png(file=glue("{base.dir}/plots/Boxplot_Age.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),], aes(factor(Group),Age.First.Immunotherapy, fill=factor(Group))) + 
stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE,label.x = 1.5) + geom_boxplot(outlier.shape=NA) + #label.x = 2.45
labs(x = "",fill = "Group :",y = "")  + geom_jitter( position=position_jitter(0.2))+ 
scale_fill_manual(values= c( "1"= "#00BFC4", "2"= "#F8766D"))+ theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum() 
dev.off()

png(file=glue("{base.dir}/plots/Boxplot_TIS_Score.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),], aes(factor(Group),TIS_Score, fill=factor(Group))) +
stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE,label.x = 1.5) + geom_boxplot(outlier.shape=NA) + #label.x = 2.45
labs(x = "",fill = "Group :",y = "")  + geom_jitter( position=position_jitter(0.2))+ 
scale_fill_manual(values= c( "Oui"= "#00BFC4", "Non"= "#F8766D"))+ theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum() 
dev.off()

png(file=glue("{base.dir}/plots/Boxplot_Nb.Pa.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),], aes(factor(Group),Nb.Pa, fill=factor(Group))) + geom_boxplot() + 
stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE,label.x = 1.5) + geom_boxplot(outlier.shape=NA) + #label.x = 2.45
labs(x = "",fill = "Group :",y = "")  + geom_jitter( position=position_jitter(0.2))+ 
scale_fill_manual(values= c( "1"= "#00BFC4", "2"= "#F8766D"))+ theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum() 
dev.off()

png(file=glue("{base.dir}/plots/Barplot_Histo.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),] )+ geom_bar( aes( factor(Group),fill = factor(Histo)),color="black"   )   + 
labs(x = "Group",fill = "Histo :",y = "") + theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum()  + scale_fill_brewer(palette="Blues")
dev.off()

png(file=glue("{base.dir}/plots/Barplot_Sexe.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),] )+ geom_bar( aes( factor(Group) ,fill = factor(Sexe)) ,color="black" )  + 
labs(x = "Group",fill = "Sexe :",y = "") + theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum()  + scale_fill_brewer(palette="Blues")
dev.off()

png(file=glue("{base.dir}/plots/Barplot_Best.Response.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),] )+ geom_bar( aes( factor(Group),fill = factor(Best.Response)) ,color="black")  + 
labs(x = "Group",fill = "Best.Response :",y = "") + theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum()  + scale_fill_brewer(palette="Blues")
dev.off()

png(file=glue("{base.dir}/plots/Barplot_Tabaco.png"),width=400,height=500)
ggplot(dataframe.final[!is.na(dataframe.final$Group),] )+ geom_bar( aes( factor(Group),fill = factor(Tabaco)) ,color="black")  + 
labs(x = "Group",fill = "Tabaco :",y = "") + theme(legend.position="right"   ,axis.text.y = element_text(size=14)) + theme_ipsum()  + scale_fill_brewer(palette="Blues")
dev.off()



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
