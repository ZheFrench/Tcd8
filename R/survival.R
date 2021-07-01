suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggcorrplot))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rstatix))
library(reshape2)
dataframe.clinical   <- fread(glue("/data/villemin/data/Tcd8/clinicData.tsv"),data.table=F)
TIS <-fread(glue("/data/villemin/data/Tcd8/TIS.Score.tsv"),data.table = F)
names(TIS)[1] <- "Patient"

option_list = list(
  make_option(c("-a", "--analysis"), type="character",default="OS",  help="Analyse", metavar="ANALYSE"),
  make_option(c("-c", "--cutoff"), type="character", default="tertile" ,  help="Cutoff", metavar="CUTOFF")

)


# clinicData.tsv  
# NanoString.normalised.tsv
# TIS.Score.tsv
# /data/villemin/data/Tcd8/ TIS.Score.tsv

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = TRUE);
opt   <- arguments$options
args  <- arguments$args

print("> OPTS : ")
print("> ARGS : ")
print(args)

analysis = opt$analysis
cutoff   = opt$cutoff 
print(analysis)
print(cutoff)
dir.create(glue("/data/villemin/data/Tcd8/plots/survival/{cutoff}_{analysis}/"), showWarnings = F)

# jour/mois/année
base.dir <- "/data/villemin/data/Tcd8/"

set.seed(122)

percentage <- c("percentage.Exhausted.Total",
"percentage.Trm.Stroma",
"percentage.Trm.Tumeur","percentage.Trm.Total")
#"percentage.Trm.Total"
#"density.Trm.Total"
density <- c("density.Exhausted.Total",
"density.Trm.Stroma",
"density.Trm.Tumeur","density.Trm.Total"
)

for (file in c(percentage)){ 
    
    item <-  deparse(substitute(percentage))

    print(file)
    dataframe <-fread(glue("{base.dir}raw_infiltration/{file}.tsv"),data.table = F)
    
    dataframe  <- na.omit(dataframe) 

    rownames(dataframe) <- dataframe$Patient
    dataframe$Patient <- NULL

    dataframe.Correlation <- cor(dataframe)
    head(round(dataframe.Correlation,2))

    png(file=glue("{base.dir}/plots/{file}.png"),width=1000,height=1000)
    print(ggcorrplot(dataframe.Correlation,title = file, lab = TRUE, insig = "blank", lab_size = 8,type = "upper", hc.order = TRUE) + theme_ipsum(base_size = 18) + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)))
    dev.off()
    
   dataframe$Patient    <-  rownames(dataframe)  
   rownames(dataframe)  <- NULL
    
   if (!exists("dataset.merged") ){
    dataset.merged <- dataframe
    next
  }
    
   if (exists("dataset.merged")){

        temp_dataset   <- dataframe
        dataset.merged <- inner_join(dataset.merged,temp_dataset, by = c("Patient"),keep = FALSE)

        rm(temp_dataset)
    }
}

# When you want to integrate the tis
#dataset.merged <- inner_join(dataset.merged,TIS, by = c("Patient"),keep = FALSE)

rownames(dataset.merged) <- dataset.merged$Patient
dataset.merged$Patient <- NULL

dataframe.Correlation <- cor(dataset.merged)

write.table(cor_gather(dataframe.Correlation) ,file = glue("/data/villemin/data/Tcd8/plots/{item}_correlation.csv"), row.names =  FALSE,col.names = TRUE, quote = F,sep = "\t")

#dataframe.Correlation
data <- ggcorrplot(dataframe.Correlation,title =item, lab = TRUE, insig = "blank", lab_size = 8,type = "upper", hc.order = TRUE) 
# To get data reorderered by clustering for coloration
pg <- ggplot_build(data)

#head(pg$Var1)
y <- ifelse(grepl("Stroma",pg$layout$panel_params[[1]]$y$get_labels()
), "blue", "red")
y <- ifelse(grepl("Total",pg$layout$panel_params[[1]]$y$get_labels()
), "black",y)
y <- ifelse(grepl("TIS",pg$layout$panel_params[[1]]$y$get_labels()
), "green",y)
y <- ifelse(grepl("Ex",pg$layout$panel_params[[1]]$y$get_labels()
), "purple",y)

x <- ifelse(grepl("Stroma",pg$layout$panel_params[[1]]$x$get_labels()
), "blue", "red")
x <- ifelse(grepl("Total",pg$layout$panel_params[[1]]$x$get_labels()
), "black",x)
x <- ifelse(grepl("TIS",pg$layout$panel_params[[1]]$x$get_labels()
), "green",x)
x <- ifelse(grepl("Ex",pg$layout$panel_params[[1]]$x$get_labels()
), "purple",x)

png(file = glue("{base.dir}/plots/Corr_global_{item}.png"),width = 2200,height = 2000)
print(ggcorrplot(dataframe.Correlation,title = item, lab = TRUE, insig = "blank", lab_size = 8,type = "upper", hc.order = TRUE) + theme_ipsum(base_size = 18) + theme(axis.text.y = element_text(colour = y) ,axis.text.x = element_text(angle = 30, vjust = 1, hjust=1,colour = x)))
dev.off()

dataset.merged$Patient <- rownames(dataset.merged)
rownames(dataset.merged) <- NULL
names(dataframe.clinical)[1] <- "Patient"

# TO BE SURE YOU HAVE THE CLINICAL DATA, AND COMPUTE TRUE QUANTILE 
dataset.merged <-  inner_join(x = dataset.merged , y = dataframe.clinical ,  by = "Patient")
dataset.merged <- dataset.merged %>% select(-Date.Biopsy,-Date.Birth,-Death.binary ,Date.Progression.or.Endpoint , -Date.Death.Or.Last.Contact , -Date.First.Immunotherapy ,-Tabaco ,-Nb.Pa ,-Sexe ,-Histo ,-Best.Response ,-Responder.1 ,-Responder.2 ,-Age.First.Immunotherapy)
dim(dataset.merged)

# You check that there is a date for PFS if not there will be bug in number of observations between boxplot and survival
if (analysis == "PFS"){
    dataset.merged <- dataset.merged[!is.na(dataset.merged$Date.Progression.or.Endpoint), ] 
}


count = 0
cutoffs <-  data.frame("VAR"=character(),"FIRST"=double(), "MEDIAN"=double(), "TERTILE"=double(),"LOW"=integer(),"HIGH"=integer(),"EFFECTIF"=integer()) 


for (feature in  names(dataset.merged) ) {
    
    if (feature == "Patient"){next}
    if (feature == "Date.Progression.or.Endpoint"){next}
    print('======================>>>')
    print(feature)

    # I know I do that several times......
    meltData <- melt(dataset.merged)
    names(meltData)[names(meltData) == "variable"] <- "feature"

    png(file=glue("{base.dir}/plots/survival/BoxplotDescriptif.png"),width = 1300,height = 1000)
    print(ggplot(meltData, aes(factor(feature), value)) +
    geom_boxplot(outlier.shape=NA)  +                                              
    labs(x = "")  + geom_jitter( position=position_jitter(0.2)) +
    theme(legend.position="right"   ,  axis.text.x = element_text(angle = 90,size=14, vjust = 1, hjust=1),axis.text.y = element_text(size=14)  )  + coord_flip() )
    dev.off()
    
    first   <- quantile(dataset.merged[[feature]], .30,na.rm = TRUE) #get 1st  tertile
    last    <- quantile(dataset.merged[[feature]], .70,na.rm = TRUE) #get last tertile .70
    median  <- quantile(dataset.merged[[feature]], .50,na.rm = TRUE) #get last tertile .70

    first   <- unname(first)
    last    <- unname(last)
    median  <- unname(median)
    
	if (first == last){next}
    
    print(first)
    print(last)
    #print(median)
    print(class(first))
    print(class(last))
    dataset.merged$group <- "MIDDLE"
    count = count + 1
    #print(dataset.merged[feature])
    print(class(dataset.merged))
    if (cutoff == "tertile" ) { 
        print("TERTILE ")
        dataset.merged$group[dataset.merged[[feature]] <= first] <- "low" 
        dataset.merged$group[dataset.merged[[feature]] >= last]  <- "high"
    }
    
    #print(dataset.merged[,c("Patient","group",all_of(feature))] %>% arrange(dataset.merged[[feature]]))
    
    if (cutoff == "median" ) { 
        print("MEDIAN ")
        dataset.merged$group[dataset.merged[[feature]] < as.numeric(median)] <- "low" 
        dataset.merged$group[dataset.merged[[feature]] > as.numeric(median)] <- "high"
    }
    
    subset <- dataset.merged %>% select(Patient, group  , all_of(feature)) %>% arrange(dataset.merged[[feature]])
    #write.table(subset  ,file = glue("/data/villemin/data/Tcd8/plots/survival/{cutoff}_{analysis}/{feature}_subset.csv"), row.names =  FALSE, quote = F,sep = "\t")
   
    subset <- subset[subset$group != "MIDDLE",]
    
    print(dim(subset))
    low_count <- nrow(subset[subset$group == "low",])
    print (low_count)

    high_count <- nrow(subset[subset$group == "high",])
    print (high_count)
  
    cutoffs.2add <- c(feature,first,median,last,low_count,high_count,dim(dataset.merged)[1])
    cutoffs <- rbind(cutoffs, cutoffs.2add)

    xlabs <- paste(levels(subset$group),"\n( N = ",table(subset$group)," )",sep="")
    png(file=glue("{base.dir}/plots/survival/{cutoff}_{analysis}/_{feature}_Boxplot.png"),width = 400,height = 500)
    print(ggplot(subset, aes(factor(group),get(feature), fill=factor(group))) +
    stat_compare_means(label="p.signif",method = "wilcox.test", paired = FALSE,label.x = 1.5) +     geom_boxplot(outlier.shape=NA) + scale_x_discrete(labels=xlabs)  +                                              
    labs(title = feature,x = "",fill = "Group :",y = "")  + geom_jitter( position=position_jitter(0.2)) +
    scale_fill_manual(values= c( "high"= "#F8766D", "low"= "#00BFC4"))+ theme(legend.position="right"   ,axis.text.y = element_text(size=14)  ) + theme_ipsum(  plot_title_size = 10) )
    dev.off()
    
    write.table(subset %>% select(Patient, group, all_of(feature)) ,file = glue("/data/villemin/data/Tcd8/plots/survival/{cutoff}_{analysis}/{feature}_annotation.csv"), row.names =  FALSE, quote = F,sep = "\t")

    print(glue("Rscript survival_panel.R -f /data/villemin/data/Tcd8/plots/survival/{cutoff}_{analysis}/{feature}_annotation.csv -v {feature}  -a {analysis} -c {cutoff} "))
    system(glue("Rscript survival_panel.R -f /data/villemin/data/Tcd8/plots/survival/{cutoff}_{analysis}/{feature}_annotation.csv -v {feature} -a {analysis} -c {cutoff}  "))
     #if (count == 2 ) { break() } 
}

cutoffs <- data.frame(setNames(cutoffs, c("VAR", "FIRST", "MEDIAN","TERTILE","LOW","HIGH","EFFECTIF")))

head(cutoffs)
write.table(cutoffs ,file = glue("/data/villemin/data/Tcd8/plots/survival/{item}_{analysis}_{cutoff}_cutoffs.csv"), col.names= TRUE,row.names =  FALSE, quote = F,sep = "\t")

print("===> ")
print(glue("{count} features were tested."))