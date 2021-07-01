#!/bin/bash
#################################################################
#
#date: July 01, 2021
#platform: Ubuntu 18.04
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# doit.sh
# Usage :
#
# ./doit.sh
#
# Description :
#
# Do everything we have done until now for survival analysis & Co.
#
################

DIR=/data/villemin/data/Tcd8/raw_nanostring

###################################################
#Start
##################################################
# This must be run two times. One variable is hard coded and need to be change manually each time it runs.
# Density / Percentage
#Rscript  /data/villemin/code/Tcd8/R/diver.R -a "OS" -c "median"
#Rscript  /data/villemin/code/Tcd8/R/diver.R -a "PFS" -c "median"

#Rscript /data/villemin/code/Tcd8/R/diver.R -a "OS" -c "tertile"
#Rscript  /data/villemin/code/Tcd8/R/diver.R -a "PFS" -c "tertile"

#CD8Plus_TCF1Plus.Stroma.Trm
VAR=CD8Plus_CD49aPlus_CD103Plus_TCF1Moins.Tumeur.Trm
###################################################
#Start
###################################################
echo ${VAR}_annotation.csv
echo ''
#=>>>>>>>>>>>>>>>>>>>>>>>>
echo "Add clinical data & expression for group of patients define for survival analysis"
Rscript /data/villemin/code/Tcd8/R/joinAnnotation.R -a /data/villemin/data/Tcd8/plots/survival/DENSITY/tertile_OS/${VAR}_annotation.csv
# This create /data/villemin/data/Tcd8/plots/survival/DENSITY/tertile_OS/heatmap.CD8Plus_CD49aPlus_CD103Plus_TCF1Moins.Tumeur.Trm_annotation.csv.annotated.tsv

Rscript /data/villemin/code/Tcd8/R/joinAnnotation.R -a /data/villemin/data/Tcd8/plots/survival/DENSITY/median_OS/${VAR}_annotation.csv
#/data/villemin/data/Tcd8/plots/survival/DENSITY/median_OS/heatmap.CD8Plus_CD49aPlus_CD103Plus_TCF1Moins.Tumeur.Trm_annotation.csv.annotated.tsv
exit 0
#=>>>>>>>>>>>>>>>>>>>>>>>>
echo "Wilcoxon of normalised gene expression values between groups"
Rscript /data/villemin/code/Tcd8/R/diffsimple.R -a /data/villemin/data/Tcd8/plots/survival/DENSITY/median_OS/${VAR}_annotation.csv
# This create /data/villemin/data/Tcd8/heatmap.CD8Plus_CD49aPlus_CD103Plus_TCF1Moins.Tumeur.Trm_annotation.csv.annotated.tsv

Rscript /data/villemin/code/Tcd8/R/diffsimple.R -a /data/villemin/data/Tcd8/plots/survival/DENSITY/tertile_OS/${VAR}_annotation.csv
#CD8Plus_CD49aPlus_CD103Moins.Total.Trm_annotation.csv-filtred-diffsimple.tsv

#=>>>>>>>>>>>>>>>>>>>>>>>>
echo "Re-Add annotation on filtered subset of data"
# Add Annotation to filtered stuffs
Rscript /data/villemin/code/Tcd8/R/joinAnnotation.R -a /data/villemin/data/Tcd8/plots/survival/DENSITY/median_OS/${VAR}_annotation.csv -e   /data/villemin/data/Tcd8/plots/survival/DENSITY/median_OS/${VAR}_annotation.csv-filtred-diffsimple.tsv

Rscript /data/villemin/code/Tcd8/R/joinAnnotation.R -a /data/villemin/data/Tcd8/plots/survival/DENSITY/tertile_OS/${VAR}_annotation.csv -e   /data/villemin/data/Tcd8/plots/survival/DENSITY/tertile_OS/${VAR}_annotation.csv-filtred-diffsimple.tsv

#=>>>>>>>>>>>>>>>>>>>>>>>>
echo "Differential on raw data using EdgeR"
# Do differential expression based on raw data
Rscript /data/villemin/code/Tcd8/R/diffDE.R -a /data/villemin/data/Tcd8/plots/survival/DENSITY/median_OS/${VAR}_annotation.csv

Rscript /data/villemin/code/Tcd8/R/diffDE.R -a /data/villemin/data/Tcd8/plots/survival/DENSITY/tertile_OS/${VAR}_annotation.csv

#=>>>>>>>>>>>>>>>>>>>>>>>>
echo "GSEA"
# Gsea on previous output
Rscript /data/villemin/code/Tcd8/R/gsea.R -d /data/villemin/data/Tcd8/plots/survival/DENSITY/median_OS/DE/${VAR}annotation.csv_high_low-differential.tsv

Rscript /data/villemin/code/Tcd8/R/gsea.R -d /data/villemin/data/Tcd8/plots/survival/DENSITY/tertile_OS/DE/${VAR}annotation.csv_high_low-differential.tsv

#############################
#End
#############################
exit 0
###################################################
#Command to pull the raw data into a matrix.
##################################################
#find /data/villemin/data/Tcd8/raw_nanostring/ -name "*.RCC" | parallel -j 62 "awk 'NR > 26 && NR < 812' {1} > {.}.clean.csv"  

###################################################
# Old stuffs to filter heatmap using a custom gene list.
##################################################

#head -n=5 /data/villemin/data/Tcd8/heatmap.CD8Plus_CD49aPlus_CD103Plus_TCF1Moins.Tumeur.Trm_annotation.csv.annotated.tsv > heatmap.tsv
#grep -f signature.TRM.Mami-Chouaib.txt /data/villemin/data/Tcd8/heatmap.CD8Plus_CD49aPlus_CD103Plus_TCF1Moins.Tumeur.Trm_annotation.csv.annotated.tsv >> heatmap.tsv
