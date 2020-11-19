# TTT (Tissue Transcriptomics Toolbox)

Spatial Transcriptomics (ST) is a sequencing technology that combines RNA-Seq 
with spatial barcodes. This enables analysis of transcriptomic data in
a spatial context. TTT contains tools for preprocessing, normalization,
clustering and further analysis as well as visualisation methods 
for ST data.

## Installation
You can install this package in one of the following ways:
- with BiocManager:
```
install.packages(c("BiocManager", "remotes"))
BiocManager::install("tmirus/TTT")
```
- with devtools:
```
install.packages("devtools")
devtools::install_github("https://github.com/tmirus/TTT")
```
- manually:
```
git clone https://github.com/tmirus/TTT
cd TTT
R CMD INSTALL ./
```

The last option requires manual installation of all dependencies:
- ggplot2
- EBImage
- genlasso 
- biomaRt 
- Rtsne 
- multtest 
- cluster 
- gridExtra 
- RColorBrewer
- doParallel
- topGO
- foreach
- uwot
- igraph
- FNN
- GO.db
- sctransform
- factoextra