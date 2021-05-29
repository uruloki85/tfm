if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install(version = "3.12")
BiocManager::install(version = "3.13")

BiocManager::install("biomaRt")

BiocManager::install("GEOquery")

BiocManager::install("org.Hs.eg.db")

BiocManager::install("cBioPortalData")
BiocManager::install("AnVIL")

# locate package installation paths
.libPaths()

# install in local folder, not system's
install.packages(c("mgcv","spatial"))
# Then, manually remove them from the system folder

# remove.packages("ggplot2")

install.packages("ggplot2")
install.packages("ggrepel")

install.packages("stringr")   
# install.packages("splitstackshape")


# install in local folder, not system's
install.packages(c("lattice","survival"))
# Then, manually remove them from the system folder

install.packages("UpSetR")
install.packages("survminer")
# install.packages("pheatmap")
install.packages("reshape2")
# install.packages("randomForest")


deps <- tools::package_dependencies("ggplot2", recursive = TRUE)$ggplot2
for (dep in deps)
  try(install.packages(dep))

