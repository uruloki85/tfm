library(dplyr)
library(readr)
require(GEOquery)
library(limma)
library(ggplot2)

########################
# Treatment vs Outcome #
########################
# GSE112282, GSE45757, GSE14426

##############################
# Gene expression vs Outcome #
##############################
# GSE21501, GSE28735, GSE62165, GSE71729, GSE56560

###############################################################################
##                              Get the data                                 ##
###############################################################################

seriesName <- "GSE56560"
do_volcano_plots <- FALSE

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
if (seriesName == "GSE45757"
    || seriesName == "GSE112282"
    || seriesName ==  "GSE14426"
    || seriesName == "GSE56560") {
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

if(seriesName == "GSE112282") { #x
  sampleInfo <- dplyr::select(sampleInfo, 
                       "cell line:ch1", 
                       "replicate info:ch1", 
                       "treatment:ch1")
  sampleInfo <- dplyr::rename(sampleInfo, 
                       line="cell line:ch1", 
                       replicate="replicate info:ch1", 
                       treatment="treatment:ch1")

  } else if (seriesName == "GSE45757") { # x
  sampleInfo <- dplyr::select(sampleInfo, 
                       "treated with:ch1", 
                       "cell line:ch1")
  sampleInfo <- dplyr::rename(sampleInfo, 
                       treated="treated with:ch1", 
                       line="cell line:ch1")
  # Fix typo: replace Mpanc96 with MPanc96
  sampleInfo$line[sampleInfo$line == "Mpanc-96"] <- "MPanc-96"
  
} else if (seriesName == "GSE14426") { # x
  # Select only samples for 24h and 168h
  library("stringr")  
  sampleInfoSubset <- sampleInfo[str_detect(sampleInfo$source_name_ch1, 
                                            "24hr|168hr"), ] 
  sampleInfo <- dplyr::select(sampleInfoSubset, 
                       "source_name_ch1")
  sampleInfo <- dplyr::rename(sampleInfo, 
                       source="source_name_ch1")
  
  samples_to_keep <- row.names(sampleInfo)
} else if (seriesName == "GSE71729") { # x
  # sampleInfo <- dplyr::select(sampleInfo, 
  #                             "tissue type:ch2")
  # sampleInfo <- rename(sampleInfo,
  #                      tissue="tissue type:ch2")
  # sampleInfo <- na.omit(sampleInfo, "tissue")
  
  sampleInfo <- dplyr::select(sampleInfo, 
                              "survival_months:ch2")
  sampleInfo <- dplyr::rename(sampleInfo,
                       OS="survival_months:ch2")
  sampleInfo$OS <- as.numeric(sampleInfo$OS)
  sampleInfo <- na.omit(sampleInfo, "OS")
  
  samples_to_keep <- row.names(sampleInfo)
  
  # 3 -> 822/30= 27.4
  # 2 -> 599/30= 19.9
  # 1 -> 474÷30= 15.8
  # 4a -> 324/30= 10.8
  # 4b -> 162/30= 5.4

  # Classify into a few groups
  
  # sampleInfo$stage[sampleInfo$OS <= 5.4] <- 'High'
  # sampleInfo$stage[sampleInfo$OS <= 10.8 & sampleInfo$OS > 5.4] <- 'Medium'
  # sampleInfo$stage[sampleInfo$OS > 10.8] <- 'Low'
  
  sampleInfo$stage[sampleInfo$OS <= 10.8] <- 'Advanced'
  sampleInfo$stage[sampleInfo$OS > 10.8] <- 'Early'
  
} else if (seriesName == "GSE56560") { # x
  sampleInfo <- dplyr::select(sampleInfo,
                              "grading:ch1")
  sampleInfo <- dplyr::rename(sampleInfo,
                       grading="grading:ch1")
  sampleInfo$grading[sampleInfo$grading == "N/A"] <- NA
  sampleInfo <- na.omit(sampleInfo, "grading")
  
  samples_to_keep <- row.names(sampleInfo)

} else if (seriesName == "GSE28735") { # x
  sampleInfo <- dplyr::select(sampleInfo, 
                              "survival_month:ch1", 
                              "tissue:ch1")
  sampleInfo <- dplyr::rename(sampleInfo,
                              OS="survival_month:ch1",
                              tissue="tissue:ch1")
  sampleInfo$OS <- as.numeric(sampleInfo$OS)
  sampleInfo <- na.omit(sampleInfo, "OS")
  # Select only Tumor samples
  # sampleInfo <- sampleInfo[sampleInfo$tissue == "T", ]
  
  samples_to_keep <- row.names(sampleInfo)
  
  # sampleInfo$stage[sampleInfo$OS <= 5.4] <- 'High'
  # sampleInfo$stage[sampleInfo$OS <= 10.8 & sampleInfo$OS > 5.4] <- 'Medium'
  # sampleInfo$stage[sampleInfo$OS > 10.8] <- 'Low'

  sampleInfo$stage[sampleInfo$OS <= 10.8] <- 'Advanced'
  sampleInfo$stage[sampleInfo$OS > 10.8] <- 'Early'
  
} else if (seriesName == "GSE21501") { # x
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
  
} else if (seriesName == "GSE62165") { # x
  # sampleInfo <- select(sampleInfo,
  #                      "grouped stage:ch1")
  # sampleInfo <- rename(sampleInfo,
  #                      stage="grouped stage:ch1")
  sampleInfo <- dplyr::select(sampleInfo,
                              "Stage:ch1")
  sampleInfo <- dplyr::rename(sampleInfo,
                              stage="Stage:ch1")
  # Remove samples with NA value
  sampleInfo[sampleInfo == "NA"] <- NA
  sampleInfo <- na.omit(sampleInfo, "stage")
  
  # Convert to 4 groups: Early (1,2) & Advanced (3,4) stages
  sampleInfo$stage[sampleInfo$stage == "1a"] <- "Early"
  sampleInfo$stage[sampleInfo$stage == "1b"] <- "Early"
  sampleInfo$stage[sampleInfo$stage == "2a"] <- "Early"
  sampleInfo$stage[sampleInfo$stage == "2b"] <- "Early"
  sampleInfo$stage[sampleInfo$stage == "3"] <- "Advanced"
  sampleInfo$stage[sampleInfo$stage == "4"] <- "Advanced"
  
  # sampleInfo$stage[sampleInfo$stage == "1a"] <- "1"
  # sampleInfo$stage[sampleInfo$stage == "1b"] <- "1"
  # sampleInfo$stage[sampleInfo$stage == "2a"] <- "2"
  # sampleInfo$stage[sampleInfo$stage == "2b"] <- "2"
  # sampleInfo$stage[sampleInfo$stage == "3"] <- "3"
  # sampleInfo$stage[sampleInfo$stage == "4"] <- "4"
  
  # sampleInfo$stage[sampleInfo$stage == "1a"] <- "1"
  # sampleInfo$stage[sampleInfo$stage == "1b"] <- "1"
  # sampleInfo$stage[sampleInfo$stage == "2a"] <- "2"
  # sampleInfo$stage[sampleInfo$stage == "2b"] <- "2"
  # sampleInfo$stage[sampleInfo$stage == "3"] <- "3"
  # sampleInfo$stage[sampleInfo$stage == "4"] <- "3"
  
  samples_to_keep <- row.names(sampleInfo)
  # samples_to_keep
  
  length(unique(sampleInfo$stage)) # 2 unique values
}

#################################################
## Get the relevant metadata from the platform ##
#################################################

features <- fData(gse)
# View(features)

# We keep the probe ID and the gene accession
if(seriesName == "GSE112282" #x
   || seriesName == "GSE45757" # x
   || seriesName == "GSE21501") {  # x
  features <- dplyr::select(features, ID, GB_ACC)
} else if (seriesName == "GSE14426") { # x
  features <- dplyr::select(features, ID, Accession)
} else if (seriesName == "GSE28735"  # x
    || seriesName == "GSE62165" # x
    || seriesName == "GSE56560") { # x
  features <- dplyr::select(features, ID, GB_LIST)
} else if (seriesName == "GSE71729") { # x
  features <- dplyr::select(features, ID)
}

###############################################################################
##                                Design                                     ##
###############################################################################

############################################
## Define the design for the DEG analysis ##
############################################
if(seriesName == "GSE112282") { #x
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

} else if (seriesName == "GSE45757") { # x
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
  
} else if (seriesName == "GSE14426") { # x
  design <- model.matrix(~0+sampleInfo$source)
  design_colnames <- c("ATRA168h","ATRA24h","Vehicle168h","Vehicle24h")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(Vehicle168h - ATRA168h, 
                             Vehicle24h - ATRA24h,
                             levels=design)
  
} else if (seriesName == "GSE21501") { # x
  design <- model.matrix(~0+sampleInfo$risk)
  design_colnames <- c("HighRisk","LowRisk")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(LowRisk - HighRisk, 
                             levels=design)
  
} else if (seriesName == "GSE62165") { # x
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
  design_colnames <- c("Advanced","Early")
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(Early - Advanced,
                             levels=design)
  # design_colnames <- c("g1","g2","g3")
  # colnames(design) <- design_colnames
  # contrasts <- makeContrasts("g2 - g1",
  #                            "g3 - g1",
  #                            "g3 - g2",
  #                            levels=design)

} else if (seriesName == "GSE71729") { # x
  # design <- model.matrix(~0+sampleInfo$tissue)
  # design_colnames <- c("Metastasis", "Normal", "Primary")
  # colnames(design) <- design_colnames
  # contrasts <- makeContrasts(Normal - Metastasis,
  #                            Primary - Metastasis,
  #                            Primary - Normal,
  #                            levels=design)
  design <- model.matrix(~0+sampleInfo$stage)
  # design_colnames <- c("High", "Low","Medium")
  design_colnames <- c("Advanced", "Early")
  # design_colnames <- c("High", "Low","Medium","Tumor")
  colnames(design) <- design_colnames
  # contrasts <- makeContrasts(Low - High,
  #                            Medium - Low,
  #                            Medium - High,
  #                            levels=design)
  contrasts <- makeContrasts(Early - Advanced,
                             levels=design)
  
} else if(seriesName == "GSE28735") { # x
  # design <- model.matrix(~0+sampleInfo$tissue)
  # design_colnames <- c("Metastasis", "Normal", "Primary")
  # colnames(design) <- design_colnames
  # contrasts <- makeContrasts(Normal - Metastasis,
  #                            Primary - Metastasis,
  #                            Primary - Normal,
  #                            levels=design)
  design <- model.matrix(~0+sampleInfo$stage
                         +sampleInfo$tissue)
  design_colnames <- c("Advanced", "Early", "Tumor")
  # design_colnames <- c("High", "Low","Medium","Tumor")
  colnames(design) <- design_colnames
  # contrasts <- makeContrasts(Low - High,
  #                            Medium - Low,
  #                            Medium - High,
  #                            levels=design)
  contrasts <- makeContrasts(Early - Advanced,
                             levels=design)
  
} else if (seriesName == "GSE56560") { # x
  # design <- model.matrix(~0+sampleInfo$tissue)
  # design_colnames <- c("Adjacent", "Normal","PDAC")
  design <- model.matrix(~0+sampleInfo$grading)
  design_colnames <- c("G2", "G3")
  colnames(design) <- design_colnames
  # contrasts <- makeContrasts(PDAC - Adjacent,
  #                            PDAC - Normal,
  #                            levels=design)
  contrasts <- makeContrasts(G3 - G2,
                             levels=design)
}
design

eset <- exprs(gse)

##############################################
## Prepare the expression data if necessary ##
##############################################
if (seriesName == "GSE28735"  # x
    || seriesName == "GSE62165" # x
    || seriesName == "GSE56560") {  # x
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

if (seriesName == "GSE45757") { # x
  # Remove columns of samples without the relevant metadata
  eset <- subset(eset, select=-c(GSM1113671,GSM1113672,GSM1113673,
                                 GSM1113674,GSM1113675,GSM1113676,
                                 GSM1113809,GSM1113810,GSM1113811))
} else if (seriesName == "GSE14426"  # x
           || seriesName == "GSE21501"  # x
    || seriesName == "GSE62165"  # x
    || seriesName == "GSE71729" # x
    || seriesName == "GSE28735"  # x
    || seriesName == "GSE56560") { # x
  eset <- eset[,samples_to_keep]
}

###############################################################################
##                          Run DEG analysis with limma                      ##
###############################################################################

ncol(eset)
nrow(design)

fit <- lmFit(eset, design)
# head(fit$coefficients)

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

if (seriesName == "GSE62165") { # x
  results <- decideTests(fit2, p.value = 0.6)  
} else if (seriesName == "GSE28735") { # x
  # tried: 0.8 (177 DEG), 0.85 (219), 0.9 (2560), 1.0 (all are detected as DEG)
  results <- decideTests(fit2, p.value = 0.85)
} else if (seriesName == "GSE71729") { # x
  # tried: 0.5 (26), 0.65 (235)
  results <- decideTests(fit2, p.value = 0.65)
} else if (seriesName == "GSE56560") { # x
  # tried: 0.5 (40), 0.6 (43), 0.65 (151), 0.7 (242)
  results <- decideTests(fit2, p.value = 0.65)
} else if (seriesName == "GSE14426") { # x
  results <- decideTests(fit2, p.value = 0.3)
} else {
  results <- decideTests(fit2)
}
# results <- decideTests(fit2, p.value = 0.9)

# How many genes are differentially-expressed
table(results)

##################
## Venn Diagram ##
##################

vennDiagram(results) 

####################
## Add gene names ##
####################

# Print top 10 DEG with meaningful name
if (seriesName != "GSE28735"  # x
    && seriesName != "GSE62165" # x
    && seriesName != "GSE56560" # x
    && seriesName != "GSE71729") { # x
  
  gene_accession <- "GB_ACC"
  if (seriesName == "GSE14426") { # x
    gene_accession <- "Accession"
  }
  
  anno <- fData(gse)
  # anno
  anno <- dplyr::select(anno,all_of(gene_accession))
  fit2$genes <- anno
  
} else if (seriesName == "GSE71729") { # x
  fit2$genes <- rownames(eset)
  # fit2$genes
} else {
  eset2_df <- as.data.frame(eset2)
  eset2_df$GB_ACC <- rownames(eset2_df)
  # eset2_df$GB_ACC
  fit2$genes <- eset2_df$GB_ACC
}

topTable(fit2)

######################
## Get common genes ##
######################
# results[,0]
# results[,1]
# results[,2]
# results[,3]
# results[,4]

# make a boolean index vector based on criteria
if (seriesName == "GSE14426") { # x
  iv <- results[,1] != 0 & results[,2] != 0
} else if (seriesName == "GSE112282") { #x
  iv <- results[,1] != 0 & results[,2] != 0 & results[,3] != 0
} else if (seriesName == "GSE45757"  # x
           || seriesName == "GSE28735" # x
    || seriesName == "GSE21501" # x
    || seriesName == "GSE62165"  # x
    || seriesName == "GSE71729" # x
    || seriesName == "GSE56560") { # x
  iv <- results[,1] != 0
}

# use it to extract gene names
deg <- fit2$genes[iv]
length(deg)
# deg

###############
# Get top DEG #
###############
do_top_x = FALSE
if (seriesName == "") {
  do_top_x = TRUE
  topValue <- 800
}

if (do_top_x) {
  # Get genes and p.values
  res_df <- data.frame(fit2$genes[iv], fit2$F.p.value[iv])
  colnames(res_df) <- c("Gene","F.p.value")
  res_df$Gene[res_df$Gene == ""] <- NA
  res_df <- na.omit(res_df, "Gene")
  # Order by p.value and get top 100
  top_x <- res_df[order(res_df$F.p.value),][0:topValue,]
  deg <- top_x$Gene
  # topTable(fit2, number=30)
}

#######################
# Convert GenBank IDs #
#######################
if (seriesName == "GSE71729") { # x
  length(unique(deg))
  mapped_to_hgnc <- deg
} else {
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
  if (seriesName == "GSE62165") { # x
    p_cutoff <- 0.3
  }
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
