require(GEOquery)
require(Biobase)

# Read from file system
# data_path <- "/home/sabela/OneDrive/Personal/UOC/2020-21_02/TFM/data/treatment_vs_outcome/original"
# fileName <- "GSE112282_series_matrix.txt.gz"
# file_path <- paste(data_path, fileName, sep="/")
# file_path
# gse <- getGEO(filename=file_path,
#               GSEMatrix =TRUE, 
#               getGPL = FALSE)



harmonize_genes <- function(seriesName) {
  
  gse <- getGEO(seriesName, GSEMatrix=TRUE, getGPL = TRUE)
  gse <- gse[[1]]
  # show(gse)
  # phenoData(gse)
  # show(pData(phenoData(gse[[1]]))[1:5,c(1,6,8)])
  
  # phenoData(gse[[1]])
  
  # phenoData(gse)
  # gse@phenoData
  # featureData(gse)
  # gse@featureData@data
  
  # Convert feature IDs to gene IDs
  featureData <- as.data.frame(gse@featureData@data)
  cols1 <- c('ID','GB_ACC')
  cols2 <- c('ID','GB_LIST')
  cols3 <- c('ID','Accession')
  cols <- ifelse("GB_LIST" %in% colnames(featureData), cols2,
                 ifelse("Accession" %in% colnames(featureData), cols3,
                        cols1))
  print(cols)
  featureData[,cols]
  
  ID_mapping <- featureData[,cols]
  
  ## exprs get the expression levels as a data frame and get the distribution
  summary(exprs(gse))
  
  exprs(gse) <- log2(exprs(gse))
  boxplot(exprs(gse),outline=FALSE)
  
  # Get expression data
  exp_data <- gse@assayData[["exprs"]]
  harmonized_exp_data <- cbind(ID_mapping, exp_data)
  
  
  
  # Save as CSV
  write.table(harmonized_exp_data,file=paste(seriesName,"exp_data_harmonized.csv", sep="_")) 
  write_csv(harmonized_exp_data,file=paste(seriesName,"exp_data_harmonized.csv", sep="_")) 
  # read.table("test.txt",header=TRUE,row.names=1) # says first column are rownames
  
  return(0)
}

# head(Meta(gse))
# names of all the GSM objects contained in the GSE
# names(GSMList(gse))

series_treatment_vs_outcome <- c("GSE37645", "GSE45757", "GSE112282", "GSE14426")
series_gene_exp_vs_outcome <- c("GSE28735", "GSE62452", "GSE21501", "GSE62165", "GSE57495")

harmonize_genes(series_treatment_vs_outcome[0])

for(i in series_treatment_vs_outcome) {
  harmonize_genes(i)
}

for(i in series_gene_exp_vs_outcome) {
  harmonize_genes(i)
}
