library(dplyr)
library(readr)
require(GEOquery)
library(limma)
library(ggplot2)

########################
# Treatment vs Outcome #
########################
# Done: GSE112282, GSE45757
# New: GSE71729, GSE56560, GSE55643, GSE77858, GSE11838, GSE16515

# Weird result: GSE60980

# Too few DEG: GSE14426

# Weird boxplot: GSE17891
# Discarded: GSE37645

##############################
# Gene expression vs Outcome #
##############################
# Done: GSE28735, GSE62452, GSE21501, GSE15471

# Weird result: GSE62165
# Discarded: GSE57495


###############################################################################
##                              Get the data                                 ##
###############################################################################

seriesName <- "GSE16515"
do_volcano_plots <- TRUE

gse <- getGEO(seriesName, GSEMatrix=TRUE, getGPL=TRUE)
gse <- gse[[1]]
# show(gse)

######################
## Explore the data ##
######################

## exprs get the expression levels as a data frame and get the distribution
summary(exprs(gse))
# exprs(gse)

# if data is log2, will be between 0 and 16
if (seriesName == "GSE37645" || seriesName == "GSE45757" 
    || seriesName == "GSE112282" || seriesName ==  "GSE14426"
    || seriesName == "GSE89396" #|| seriesName == "GSE60980"
    || seriesName == "GSE56560" || seriesName == "GSE55643"
    || seriesName == "GSE11838" || seriesName == "GSE16515") {
  exprs(gse) <- log2(exprs(gse))
}

# Verify data has been normalized
boxplot(exprs(gse), 
        outline=FALSE, 
        las=2, 
        boxwex=0.6,
        cex.axis = 0.7)

###############################################################################
##                              Metadata                                     ##
###############################################################################

################################################
## Get the relevant metadata from the samples ##
################################################

# Get sample info
sampleInfo <- pData(gse)
# View(sampleInfo)

# Select some columns

# Treatment vs outcome
if(seriesName == "GSE112282") {
  sampleInfo <- dplyr::select(sampleInfo, 
                       "cell line:ch1", 
                       "replicate info:ch1", 
                       "treatment:ch1")
  sampleInfo <- dplyr::rename(sampleInfo, 
                       line="cell line:ch1", 
                       replicate="replicate info:ch1", 
                       treatment="treatment:ch1")
}
# if (seriesName == "GSE37645") {
#   sampleInfo <- dplyr::select(sampleInfo, 
#                        "characteristics_ch1.3")
#   sampleInfo <- dplyr::rename(sampleInfo, 
#                        sensitive="characteristics_ch1.3")
# }
if (seriesName == "GSE45757") {
  sampleInfo <- dplyr::select(sampleInfo, 
                       "treated with:ch1", 
                       "cell line:ch1")
  sampleInfo <- dplyr::rename(sampleInfo, 
                       treated="treated with:ch1", 
                       line="cell line:ch1")
  # Fix typo: replace Mpanc96 with MPanc96
  sampleInfo$line[sampleInfo$line == "Mpanc-96"] <- "MPanc-96"
}
if (seriesName == "GSE14426") {
  # Select only samples for 24h and 168h
  library("stringr")  
  sampleInfoSubset <- sampleInfo[str_detect(sampleInfo$title, "24hr|168hr"), ] 
  sampleInfo <- dplyr::select(sampleInfoSubset, 
                       "source_name_ch1")
  sampleInfo <- dplyr::rename(sampleInfo, 
                       source="source_name_ch1")
  
  samples_to_keep <- row.names(sampleInfo)
}
if (seriesName == "GSE71729") {
  sampleInfo <- dplyr::select(sampleInfo, 
                              "tissue type:ch2")
  sampleInfo <- rename(sampleInfo,
                       tissue="tissue type:ch2")
  sampleInfo <- na.omit(sampleInfo, "tissue")
  
  samples_to_keep <- row.names(sampleInfo)
}
if (seriesName == "GSE60980") {
  sampleInfo <- dplyr::select(sampleInfo, 
                              "tissue:ch1")
  sampleInfo <- rename(sampleInfo,
                       tissue="tissue:ch1")
}
if (seriesName == "GSE56560") {
  sampleInfo <- dplyr::select(sampleInfo,
                              "tissue:ch1")
  sampleInfo <- rename(sampleInfo,
                       tissue="tissue:ch1")
}
if (seriesName == "GSE55643") {
  sampleInfo <- dplyr::select(sampleInfo,
                              "tissue:ch1")
  sampleInfo <- rename(sampleInfo,
                       tissue="tissue:ch1")
}
if (seriesName == "GSE77858") {
  sampleInfo <- dplyr::select(sampleInfo,
                              "morphology:ch2")
  sampleInfo <- rename(sampleInfo,
                       tissue="morphology:ch2")
  # Fix typo
  sampleInfo[sampleInfo == "Panreatitis"] <- "Pancreatitis"
  # unique(sampleInfo$tissue)
}
if (seriesName == "GSE11838") {
  sampleInfo <- dplyr::select(sampleInfo,
                              "source_name_ch1")
  sampleInfo <- rename(sampleInfo,
                       tissue="source_name_ch1")
  # unique(sampleInfo$tissue)
}
if (seriesName == "GSE16515") {
  sampleInfo <- dplyr::select(sampleInfo,
                              "tissue:ch1")
  sampleInfo <- rename(sampleInfo,
                       tissue="tissue:ch1")
  # unique(sampleInfo$tissue)
}

# Gene expression vs outcome
if (seriesName == "GSE28735" || seriesName == "GSE62452") {
  sampleInfo <- dplyr::select(sampleInfo, "tissue:ch1")
  sampleInfo <- dplyr::rename(sampleInfo, 
                       tissue="tissue:ch1")
}
if (seriesName == "GSE21501") {
  sampleInfo <- dplyr::select(sampleInfo, 
                       "characteristics_ch2.5", 
                       "characteristics_ch2.6")
  sampleInfo <- dplyr::rename(sampleInfo, 
                       risk="characteristics_ch2.5", 
                       risk2="characteristics_ch2.6")
  
  # Information is misplaced in these samples
  sampleInfo["GSM536946","risk"] <- sampleInfo["GSM536946","risk2"]
  sampleInfo["GSM536892","risk"] <- sampleInfo["GSM536892","risk2"]
  sampleInfo <- dplyr::select(sampleInfo, risk)
  
  # Remove samples with empty value
  sampleInfo[sampleInfo == ""] <- NA
  sampleInfo <- na.omit(sampleInfo, "risk")
  
  samples_to_keep <- row.names(sampleInfo)
  # samples_to_keep
  
  length(unique(sampleInfo$risk)) # 2 unique values
}
if (seriesName == "GSE62165") {
  # sampleInfo <- select(sampleInfo,
  #                      "grouped stage:ch1")
  # sampleInfo <- rename(sampleInfo,
  #                      stage="grouped stage:ch1")
  sampleInfo <- select(sampleInfo,
                       "Stage:ch1")
  sampleInfo <- rename(sampleInfo,
                       stage="Stage:ch1")
  # Remove samples with NA value
  sampleInfo[sampleInfo == "NA"] <- NA
  sampleInfo <- na.omit(sampleInfo, "stage")
  
  # Convert to 4 groups: Early (1,2) & Advanced (3,4) stages
  # sampleInfo$stage[sampleInfo$stage == "1a"] <- "Early"
  # sampleInfo$stage[sampleInfo$stage == "1b"] <- "Early"
  # sampleInfo$stage[sampleInfo$stage == "2a"] <- "Early"
  # sampleInfo$stage[sampleInfo$stage == "2b"] <- "Early"
  # sampleInfo$stage[sampleInfo$stage == "3"] <- "Advanced"
  # sampleInfo$stage[sampleInfo$stage == "4"] <- "Advanced"
  
  # sampleInfo$stage[sampleInfo$stage == "1a"] <- "1"
  # sampleInfo$stage[sampleInfo$stage == "1b"] <- "1"
  # sampleInfo$stage[sampleInfo$stage == "2a"] <- "2"
  # sampleInfo$stage[sampleInfo$stage == "2b"] <- "2"
  # sampleInfo$stage[sampleInfo$stage == "3"] <- "3"
  # sampleInfo$stage[sampleInfo$stage == "4"] <- "4"
  
  sampleInfo$stage[sampleInfo$stage == "1a"] <- "1"
  sampleInfo$stage[sampleInfo$stage == "1b"] <- "1"
  sampleInfo$stage[sampleInfo$stage == "2a"] <- "2"
  sampleInfo$stage[sampleInfo$stage == "2b"] <- "2"
  sampleInfo$stage[sampleInfo$stage == "3"] <- "3"
  sampleInfo$stage[sampleInfo$stage == "4"] <- "3"
  
  samples_to_keep <- row.names(sampleInfo)
  # samples_to_keep
  
  length(unique(sampleInfo$stage)) # 2 unique values
}
if (seriesName == "GSE15471") {
  sampleInfo <- dplyr::select(sampleInfo, "sample:ch1")
  sampleInfo <- rename(sampleInfo, 
                       tissue="sample:ch1")
}
# if (seriesName == "GSE57495") {
#   # sampleInfo <- select(sampleInfo, "overall survival (month):ch1")
#   # sampleInfo <- rename(sampleInfo, 
#   #                      OS="overall survival (month):ch1")
#   # sampleInfo$OS <- as.numeric(sampleInfo$OS)
#   # # Low risk -> median of median survival of 35 months
#   # # High risk -> median of median survival of 15 months
#   # sampleInfo$risk[sampleInfo$OS <= 23] <- "High"
#   # sampleInfo$risk[sampleInfo$OS > 23] <- "Low"
#   # sampleInfo <- select(sampleInfo, risk)
#   
#   sampleInfo <- select(sampleInfo, "Stage:ch1")
#   sampleInfo <- rename(sampleInfo, 
#                        stage="Stage:ch1")
#   sampleInfo$stage[sampleInfo$stage == "1B"] <- "1"
#   sampleInfo$stage[sampleInfo$stage == "2A"] <- "2"
#   sampleInfo$stage[sampleInfo$stage == "2B"] <- "2"
# }

#################################################
## Get the relevant metadata from the platform ##
#################################################

features <- fData(gse)
# View(features)

# We keep the probe ID and the gene accession
if(seriesName == "GSE112282" || seriesName == "GSE37645" 
   || seriesName == "GSE45757" || seriesName == "GSE21501"
   || seriesName == "GSE89396" || seriesName == "GSE15471"
   || seriesName == "GSE60980" || seriesName == "GSE55643"
   || seriesName == "GSE77858" || seriesName == "GSE11838"
   || seriesName == "GSE16515") {
  features <- dplyr::select(features, ID, GB_ACC)
}
if (seriesName == "GSE14426") {
  features <- dplyr::select(features, ID, Accession)
}
if (seriesName == "GSE28735" || seriesName == "GSE62452"
    || seriesName == "GSE62165" #|| seriesName == "GSE57495"
    || seriesName == "GSE56560") {
  features <- dplyr::select(features, ID, GB_LIST)
}

if (seriesName == "GSE71729") {
  features <- dplyr::select(features, ID)
}

# full_output <- cbind(features,exprs(gse))
# write_csv(full_output,file=paste(seriesName,"exp_data_harmonized.csv", sep="_")) 

###############################################################################
##                                Design                                     ##
###############################################################################

############################################
## Define the design for the DEG analysis ##
############################################
if(seriesName == "GSE112282") {
  design_colnames <- c("BET","BETMEK","MEK","VEHICLE","COLO201",
                       "HPAFII","NCIH510","RKO","Replicate2")
  design <- model.matrix(~0+sampleInfo$treatment
                         +sampleInfo$line
                         +sampleInfo$replicate)
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(BET - VEHICLE, 
                             BETMEK - VEHICLE, 
                             MEK - VEHICLE, 
                             levels=design)
}
# if (seriesName == "GSE37645") {
#   design <- model.matrix(~0+sampleInfo$sensitive)
#   colnames(design) <- c("NonSensitive","Sensitive")
#   contrasts <- makeContrasts(Sensitive - NonSensitive, levels=design)
# }
if (seriesName == "GSE45757") {
  design_colnames <- c("Treated","Untreated","Capan2","CFPAC1","COLO357",
                       "HPAFII","Hs766T","L33","L36pl","L36sl","MIAPaCa2",
                       "MPanc96","Panc1","Panc1005","Panc203","Panc213",
                       "Panc327","Panc504","Panc603","Panc813","PL45",
                       "SU8686","SW1990")
  design <- model.matrix(~0+sampleInfo$treated
                         +sampleInfo$line)
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(Untreated - Treated, 
                             levels=design)
}
if (seriesName == "GSE14426") {
  design <- model.matrix(~0+sampleInfo$source)
  design_colnames <- c("ATRA168h","ATRA24h","Vehicle168h","Vehicle24h")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(Vehicle168h - ATRA168h, 
                             Vehicle24h - ATRA24h,
                             levels=design)
}
if (seriesName == "GSE89396") {
  design <- model.matrix(~0+sampleInfo$treatment)
  colnames(design) <- c("treated","untreated")
  contrasts <- makeContrasts(untreated - treated,
                             levels=design)
}
if (seriesName == "GSE28735" || seriesName == "GSE62452") {
  design <- model.matrix(~0+sampleInfo$tissue)
  design_colnames <- c("NonTumour","Tumour")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(Tumour - NonTumour, 
                             levels=design)
}
if (seriesName == "GSE21501") {
  design <- model.matrix(~0+sampleInfo$risk)
  design_colnames <- c("HighRisk","LowRisk")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(LowRisk - HighRisk, 
                             levels=design)
}
if (seriesName == "GSE62165") {
  design <- model.matrix(~0+sampleInfo$stage)
  # design_colnames <- c("g1a","g1b","g2a","g2b","g3","g4")
  # colnames(design) <- design_colnames
  # contrasts <- makeContrasts("g1b - g1a",
  #                            "g2a - g1a",
  #                            "g2b - g1a",
  #                            "g3 - g1a",
  #                            "g4 - g1a",
  #                            "g2a - g1b",
  #                            "g2b - g1b",
  #                            "g3 - g1b",
  #                            "g4 - g1b",
  #                            "g2b - g2a",
  #                            "g3 - g2a",
  #                            "g4 - g2a",
  #                            "g3 - g2b",
  #                            "g4 - g2b",
  #                            "g4 - g3",
  #                            levels=design)
  # design_colnames <- c("g1","g2","g3","g4")
  # colnames(design) <- design_colnames
  # contrasts <- makeContrasts("g3 - g1",
  #                            "g3 - g2",
  #                            "g4 - g1",
  #                            "g4 - g2",
  #                            levels=design)
  # design_colnames <- c("Advanced","Early","LNM")
  # colnames(design) <- design_colnames
  # contrasts <- makeContrasts(Early - Advanced,
  #                            LNM - Advanced,
  #                            LNM - Early,
  #                            levels=design)
  # design_colnames <- c("Advanced","Early")
  # colnames(design) <- design_colnames
  # contrasts <- makeContrasts(Early - Advanced,
  #                            levels=design)
  design_colnames <- c("g1","g2","g3")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts("g2 - g1",
                             "g3 - g1",
                             "g3 - g2",
                             levels=design)
}
if (seriesName == "GSE15471") {
  design <- model.matrix(~0+sampleInfo$tissue)
  design_colnames <- c("Normal","Tumor")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(Tumor - Normal, 
                             levels=design)
}
if (seriesName == "GSE71729") {
  design <- model.matrix(~0+sampleInfo$tissue)
  design_colnames <- c("Metastasis", "Normal", "Primary")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(Normal - Metastasis,
                             Primary - Metastasis,
                             Primary - Normal,
                             levels=design)
}
if (seriesName == "GSE56560") {
  design <- model.matrix(~0+sampleInfo$tissue)
  design_colnames <- c("Adjacent", "Normal","PDAC")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(PDAC - Adjacent,
                             PDAC - Normal,
                             levels=design)
}
if (seriesName == "GSE55643") {
  design <- model.matrix(~0+sampleInfo$tissue)
  design_colnames <- c("Normal","Tumor")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(Tumor - Normal,
                             levels=design)
}
if (seriesName == "GSE77858") {
  design <- model.matrix(~0+sampleInfo$tissue)
  design_colnames <- c("Normal","Pancreatitis","Tumor")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(Tumor - Normal,
                             Tumor - Pancreatitis,
                             levels=design)
}
if (seriesName == "GSE11838") {
  design <- model.matrix(~0+sampleInfo$tissue)
  design_colnames <- c("Normal","Tumor")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(Tumor - Normal,
                             levels=design)
}
if (seriesName == "GSE16515") {
  design <- model.matrix(~0+sampleInfo$tissue)
  design_colnames <- c("Normal","Tumor")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(Tumor - Normal,
                             levels=design)
}
# if (seriesName = "GSE60980") {
#   design <- model.matrix(~0+sampleInfo$tissue)
#   design_colnames <- c("Normal", "AmpullaIntestinal", "AmpullaPanc", "Bile", "Duodenal","PDAC")
#   colnames(design) <- design_colnames
#   contrasts <- makeContrasts(PDAC - Normal,
#                              PDAC - AmpullaIntestinal,
#                              PDAC - AmpullaPanc,
#                              PDAC - Bile,
#                              PDAC - Duodenal,
#                              levels=design)
# }
# if(seriesName == "GSE57495") {
#   # design <- model.matrix(~0+sampleInfo$risk)
#   # design_colnames <- c("HighRisk","LowRisk")
#   # colnames(design) <- design_colnames
#   # contrasts <- makeContrasts(LowRisk - HighRisk, 
#   #                            levels=design)
#   design <- model.matrix(~0+sampleInfo$stage)
#   # design_colnames <- c("g1","g1b","g2a","g2b")
#   # colnames(design) <- design_colnames
#   # contrasts <- makeContrasts("g1b - g1",
#   #                            "g2a - g1",
#   #                            "g2a - g1b",
#   #                            "g2b - g1",
#   #                            "g2b - g1b",
#   #                            "g2b - g2a",
#   #                            levels=design)
#   design_colnames <- c("g1","g2")
#   colnames(design) <- design_colnames
#   contrasts <- makeContrasts("g2 - g1",
#                              levels=design)
# }
design

eset <- exprs(gse)

##############################################
## Prepare the expression data if necessary ##
##############################################
if (seriesName == "GSE28735" || seriesName == "GSE62452"
    || seriesName == "GSE62165" #|| seriesName == "GSE57495"
    || seriesName == "GSE56560") {
  library(tidyr)
  
  # Split rows by Gene in GB_LIST duplicating the information

  eset_df <- as.data.frame(exprs(gse))
  # Add a new column with the GB_LIST  
  eset_df <- cbind(features$GB_LIST,eset_df)
  # Rename the column
  eset_df <- eset_df %>% 
    dplyr::rename(GB_LIST = "features$GB_LIST")
  # eset_df$GB_LIST
  
  # Check if there are repetitions among the probe IDs
  length(rownames(eset_df)) 
  length(unique(rownames(eset_df))) 
  
  delimiter <- ","
  if (seriesName == "GSE57495") {
    delimiter <- " "
  }
  
  # Duplicate each row as many times as elements are present in GB_LIST
  eset_df2 <- 
    eset_df %>%
    mutate(GB_LIST = strsplit(as.character(GB_LIST), delimiter)) %>%
    unnest(cols = c(GB_LIST)) %>%
    filter(GB_LIST != "")%>% 
    dplyr::rename(GB_ACC = GB_LIST) #%>%
    # select(V1, 1:-1)
  
  # Check if there is any missing value in this column  
  which(is.na(eset_df2$GB_ACC))
  
  eset2 <- limma::avereps(eset_df2[,2:ncol(eset_df2)], eset_df2$GB_ACC)
  ncol(eset2)
  nrow(eset2)
  
  eset <- eset2
  # eset["AB001736","GSM1527105"]
  # eset["AK123548","GSM1527105"]
}

##################################
## Discard samples if necessary ##
##################################

if (seriesName == "GSE45757") {
  # Remove columns of samples without the relevant metadata
  eset <- subset(eset, select=-c(GSM1113671,GSM1113672,GSM1113673,
                                 GSM1113674,GSM1113675,GSM1113676,
                                 GSM1113809,GSM1113810,GSM1113811))
}
if (seriesName == "GSE14426" || seriesName == "GSE21501" 
    || seriesName == "GSE62165" || seriesName == "GSE71729") {
  eset <- eset[,samples_to_keep]
}

###############################################################################
##                          Run DEG analysis with limma                      ##
###############################################################################

fit <- lmFit(eset, design)
head(fit$coefficients)

contrasts

fit2 <- contrasts.fit(fit, contrasts)

# Get differential expression statistics and p-values with empirical Bayes
fit2 <- eBayes(fit2)

# Results by contrast
topTable(fit2)
# topTable(fit2, coef=1)
# topTable(fit2, coef=2)
# topTable(fit2, coef=3)
# topTable(fit2, coef=4)

# Find differentially-expressed genes

if (seriesName == "GSE62165") {
  results <- decideTests(fit2, p.value = 0.3)  
} else {
  results <- decideTests(fit2)
}
# How many genes are differentially-expressed
table(results)

####################
## Add gene names ##
####################

# Print top 10 DEG with meaningful name
if (seriesName != "GSE28735" 
    && seriesName != "GSE62452"
    && seriesName != "GSE62165"
    && seriesName != "GSE56560") {
  
  gene_accession <- "GB_ACC"
  if (seriesName == "GSE14426") {
    gene_accession <- "Accession"
  }
  
  anno <- fData(gse)
  # anno
  anno <- dplyr::select(anno,all_of(gene_accession))
  fit2$genes <- anno
  
} else if (seriesName == "GSE71729") {
  fit2$genes <- rownames(eset)
  # fit2$genes
} else {
  eset2_df <- as.data.frame(eset2)
  eset2_df$GB_ACC <- rownames(eset2_df)
  # eset2_df$GB_ACC
  fit2$genes <- eset2_df$GB_ACC
}

topTable(fit2)

###############################################################################
##                              Venn Diagram                                 ##
###############################################################################

vennDiagram(results) 

######################
## Get common genes ##
######################
# results[,0]
# results[,1]
# results[,2]
# results[,3]
# results[,4]

# make a boolean index vector based on criteria
if (seriesName == "GSE14426" || seriesName == "GSE56560"
    || seriesName == "GSE77858") {
  iv <- results[,1] != 0 & results[,2] != 0
}
if (seriesName == "GSE112282" || seriesName == "GSE71729") {
  iv <- results[,1] != 0 & results[,2] != 0 & results[,3] != 0
}
if (seriesName == "GSE45757" || seriesName == "GSE28735"
    || seriesName == "GSE62452"|| seriesName == "GSE21501"
    || seriesName == "GSE15471" || seriesName == "GSE55643"
    || seriesName == "GSE11838" || seriesName == "GSE16515") {
  iv <- results[,1] != 0
}
# if (seriesName == "GSE71729") {
#   iv <- results[,1] != 0 & results[,2] != 0 & results[,3] != 0
# }

# use it to extract gene names
deg <- fit2$genes[iv]
length(deg)

#######################
# Convert GenBank IDs #
#######################
if (seriesName == "GSE71729") {
  length(unique(deg))
  mapped_to_hgnc <- deg
}
if (seriesName != "GSE71729") {
  library(org.Hs.eg.db)
  library(biomaRt)
  
  # gene_list <- unlist(fit2$genes)
  gene_list <- deg
  length(unique(gene_list))
  # 2842
  
  mapped <- select(org.Hs.eg.db, gene_list, c("ENTREZID","SYMBOL"), "ACCNUM")
  length(unique(mapped$ACCNUM))
  # 2842
  length(unique(mapped$ENTREZID))
  # 2288
  
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  # filters <- listFilters(mart)
  
  mapped_to_hgnc <- getBM(attributes=c('hgnc_symbol'),
                          filters = 'entrezgene_id',
                          mart = mart,
                          values = mapped$ENTREZID)
  length(unique(mapped_to_hgnc$hgnc_symbol))
  #mapped_to_hgnc
}

##################
# Save GENE list #
##################

write.table(unique(mapped_to_hgnc), 
            file=paste(seriesName,"common_genes.csv", sep="_"), 
            col.names = "HGNC", 
            row.names = FALSE,
            quote = FALSE) 

###############################################################################
##                              Volcano plot                                  #
###############################################################################

if (do_volcano_plots) {

  full_results <- topTable(fit2, coef=1, number=Inf)
  # full_results <- tibble::rownames_to_column(full_results,"ID")
  ggplot(full_results,aes(x = logFC, y=B)) + geom_point()
  
  ## change according to your needs
  p_cutoff <- 0.05
  fc_cutoff <- 1
  
  full_results %>% 
    mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
    ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()
  
  # Label some genes
  # library(ggrepel)
  # 
  # p_cutoff <- 0.05
  # fc_cutoff <- 1
  # topN <- 20
  # 
  # full_results %>% 
  #   mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  #   mutate(Rank = 1:n(), Label = ifelse(Rank < topN, GB_ACC,"")) %>% 
  #   ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")

}

# Save as CSV
# write_csv(harmonized_exp_data,file=paste(seriesName,"exp_data_harmonized.csv", sep="_")) 
# read.table("test.txt",header=TRUE,row.names=1) # says first column are rownames
