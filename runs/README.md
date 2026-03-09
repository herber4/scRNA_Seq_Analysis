
### Create multiconfig.csv

```
cd ../
vi multiconfig.csv
```


### Run Cell Ranger

```
cd ../runs
sbatch slurm_submit.sh
```

### Setup Seurat

Install Seurat within R version used on openondemand

```
ml R/4.1.2
R
install.packages('Seurat')
install.packages('ggpubr')
install.packages('patchwork')
BiocManager::install("celldex")
BiocManager::install("SingleR")
BiocManager::install("ensembldb")
BiocManager::install("scrapper")
q()
y
```




