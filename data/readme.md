# scRNA Seq analysis of 10X universal 3' Universal scRNA-Seq Data

### set up project dir

Set up working dir, fetch and unzip cellranger bins

```
cd /to/dir/home/
mkdir data bin scripts ref runs
cd bin/
curl -o cellranger-10.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=1773020233&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=cxl0cvIpifksytmQnYWEsM7dIMYr0Jq7pojfi4XCNHZsayzDggBn9xXk4kb8kaiCsBW95yFfl~gPIt8Wz0T3sQ1-z9wjHo-QtBrhLN6iQohB50bwm7qAoezDqfBA6onB4hdPUqlRHWhpPzPh0rLVc15QPsL8Zq4omcMw0oFNnGi2JuPRw5f3oXULlkvbpIoilUsYdkRw05l8wGSlzp70CcgZlDrt9VNSWhzx923TERyUzI03zYL1Na3CLNAF7GuxIT1E4NriEnLaygC8~YbiZ4klN8RN7Z7GL4dNgfOI2yqt5Uyx~WIrDjdLXwCVya-abPUIsYQUiXcFYpNunsbr2g__"
tar -xvzf cellranger-10.0.0.tar.gz 
```

Add path to bashrc

```
cd ~
vi .bash_profile
export PATH:/pbmc/bin/cellranger-10.0.0:$PATH
source .bash_profile
```
Check install

```
cellranger testrun --id=check_install
```

### Setup refs

Download and install the refs

```
cd ../ref
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -xvzf refdata-gex-GRCh38-2024-A.tar.gz
```
Download vdj ref

```
curl -O "https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0.tar.gz"
tar -xvzf refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0.tar.gz 
```

### Fetch and unzip data

```
cd data/
wget https://cf.10xgenomics.com/samples/cell-exp/6.1.0/3p_Citrate_CPT/3p_Citrate_CPT_fastqs.tar
tar -xf 3p_Citrate_CPT_fastqs.tar
```

### Create multiconfig.csv

```
cd ../
vi multiconfig.csv
```


### Create conda env
```
conda activate
conda create --name pbmc
```

Install conda packages

```

```
