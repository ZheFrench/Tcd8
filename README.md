### Workflow analysis
_________________

You will find bash and R scripts in two distincts directories.  

1. `doit.sh` will call the different bash script of the worflow.

2. `survival.R` will execute logrank test for different cutoff depending parameters. Hard coded , the set of variables. Need to be changed manually by user.  It also called another R script called survival_panel.R.

```shell
# Density / Percentage
Rscript  /data/villemin/code/Tcd8/R/diver.R -a "OS" -c "median"
Rscript  /data/villemin/code/Tcd8/R/diver.R -a "PFS" -c "median"

Rscript /data/villemin/code/Tcd8/R/diver.R -a "OS" -c "tertile"
Rscript  /data/villemin/code/Tcd8/R/diver.R -a "PFS" -c "tertile"
```

3. `joinAnnotation.R` will join the output of script 1 with clinical data and expression matrix in order to visualize that in Morpheus in order to create a Matrix.
It also plot and test other values relatives to both groups.( like sexe, age, tisscore etc...)

```shell
Rscript /data/villemin/code/Tcd8/R/joinAnnotation.R -a /data/villemin/data/Tcd8/plots/survival/DENSITY/tertile_OS/CD8Plus_CD49aPlus_CD103Plus_TCF1Moins.Tumeur.Trm_annotation.csv
```

4. `diffsimple.R` will apply a simple man withney with bonferrony between groups for expression using normalised data...

```shell
Rscript /data/villemin/code/Tcd8/R/diffsimple.R -a /data/villemin/data/Tcd8/plots/survival/DENSITY/tertile_OS/CD8Plus_CD49aPlus_CD103Moins.Total.Trm_annotation.csv
```

5.  `diffsimpleDE.R` applied tools of rnaseq for differential expression... should be done on raw data....

```shell
# Do differential expression based on raw data
Rscript /data/villemin/code/Tcd8/R/diffDE.R -a /data/villemin/data/Tcd8/plots/survival/DENSITY/tertile_OS/CD8Plus_CD49aPlus_CD103Moins.Total.Trm_annotation.csv
```

6. `gsea.R` Gene set enrichment analysis.

```shell
Rscript /data/villemin/code/Tcd8/R/gsea.R -d /data/villemin/data/Tcd8/DE/CD8Plus_CD49aPlus_CD103Moins.Total.Trm_annotation.csv_high_low-differential.tsv
# Gsea on output
```