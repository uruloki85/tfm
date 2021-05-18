if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

# locate package installation paths
.libPaths()

# install in local folder, not system's
install.packages(c("mgcv","spatial"))
# Then, manually remove them from the system folder

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")


install.packages("ggplot2")
install.packages("ggrepel")


install.packages("stringr")   
# install.packages("splitstackshape")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")


BiocManager::install("cBioPortalData")
BiocManager::install("AnVIL")

# install in local folder, not system's
install.packages(c("lattice","survival"))
# Then, manually remove them from the system folder


install.packages("UpSetR")
install.packages("survminer")
install.packages("pheatmap")
install.packages("reshape2")
install.packages("randomForest")
