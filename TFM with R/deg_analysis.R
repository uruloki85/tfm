library(dplyr)
library(readr)
require(GEOquery)
library(limma)
library(ggplot2)

# series_treatment_vs_outcome <- c("GSE37645", "GSE45757", "GSE112282", "GSE14426")
# series_gene_exp_vs_outcome <- c("GSE28735", "GSE62452", "GSE21501", "GSE62165", "GSE57495")

# Done: GSE112282, GSE45757, GSE14426
# Weird result: GSE37645

seriesName <- "GSE62165"

gse <- getGEO(seriesName, GSEMatrix=TRUE, getGPL = TRUE)
gse <- gse[[1]]
# show(gse)

## exprs get the expression levels as a data frame and get the distribution
summary(exprs(gse))
# if data is log2, will be between 0 and 16
if (seriesName == "GSE37645" || seriesName == "GSE45757" 
    || seriesName == "GSE112282" || seriesName ==  "GSE14426") {
  exprs(gse) <- log2(exprs(gse))
}

# Verify data has been normalized
boxplot(exprs(gse),outline=FALSE)

# Get sample info
sampleInfo <- pData(gse)
sampleInfo

# Select some columns
if(seriesName == "GSE112282") {
  sampleInfo <- select(sampleInfo, 
                       "cell line:ch1", 
                       "replicate info:ch1", 
                       "treatment:ch1")
  sampleInfo <- rename(sampleInfo, 
                       line="cell line:ch1", 
                       replicate="replicate info:ch1", 
                       treatment="treatment:ch1")
}
if (seriesName == "GSE37645") {
  sampleInfo <- select(sampleInfo, "characteristics_ch1.3")
  sampleInfo <- rename(sampleInfo, sensitive="characteristics_ch1.3")
}
if (seriesName == "GSE45757") {
  sampleInfo <- select(sampleInfo, "treated with:ch1", "cell line:ch1")
  sampleInfo <- rename(sampleInfo, treated="treated with:ch1", line="cell line:ch1")
}
if (seriesName == "GSE14426") {
  # Select only samples for 24h and 168h
  library("stringr")  
  sampleInfoSubset <- sampleInfo[str_detect(sampleInfo$title, "24hr|168hr"), ] 
  sampleInfo <- select(sampleInfoSubset, "source_name_ch1")
  sampleInfo <- rename(sampleInfo, source="source_name_ch1")
  
  row_names <- row.names(sampleInfo)
}

features <- fData(gse)
# View(features)
# We keep the probe ID and the gene accession
if(seriesName == "GSE112282" || seriesName == "GSE37645" || seriesName == "GSE45757") {
  features <- select(features,ID,GB_ACC)
}
if (seriesName == "GSE14426") {
  features <- select(features,ID,Accession)
}
full_output <- cbind(features,exprs(gse))
# write_csv(full_output,file=paste(seriesName,"exp_data_harmonized.csv", sep="_")) 

# design <- model.matrix(~0+sampleInfo$treatment)
# colnames(design) <- c("BET","BETMEK","MEK","VEHICLE")

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
if (seriesName == "GSE37645") {
  design <- model.matrix(~0+sampleInfo$sensitive)
  colnames(design) <- c("NonSensitive","Sensitive")
  contrasts <- makeContrasts(Sensitive - NonSensitive, levels=design)
}
if (seriesName == "GSE45757") {
  # TODO Fix typo Mpanc96=MPanc96
  design_colnames <- c("Treated","Untreated","Capan2","CFPAC1","COLO357",
                       "HPAFII","Hs766T","L33","L36pl","L36sl","MIAPaCa2",
                       "Mpanc96","MPanc96","Panc1","Panc1005","Panc203",
                       "Panc213","Panc327","Panc504","Panc603","Panc813",
                       "PL45","SU8686","SW1990")
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
design

eset <- exprs(gse)

if (seriesName == "GSE45757") {
  # Remove columns of samples without the relevant metadata
  eset <- subset(eset, select=-c(GSM1113671,GSM1113672,GSM1113673,
                                 GSM1113674,GSM1113675,GSM1113676,
                                 GSM1113809,GSM1113810,GSM1113811))
}
if (seriesName == "GSE14426") {
  eset <- eset[,row_names]
}

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

# How many genes are differentially-expressesd
results <- decideTests(fit2)
table(results)

vennDiagram(results)


# Print top 10 DEG with meaningful name
gene_accession <- "GB_ACC"
if (seriesName == "GSE14426") {
  gene_accession <- "Accession"
}

anno <- fData(gse)
# anno
anno <- select(anno,all_of(gene_accession))
fit2$genes <- anno
topTable(fit2)

################
# Volcano plot #
################

full_results <- topTable(fit2, coef=1, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")
ggplot(full_results,aes(x = logFC, y=B)) + geom_point()

## change according to your needs
p_cutoff <- 0.05
fc_cutoff <- 1

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()

# Label some genes
library(ggrepel)

p_cutoff <- 0.05
fc_cutoff <- 1
topN <- 20

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, GB_ACC,"")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")


# Save as CSV
# write_csv(harmonized_exp_data,file=paste(seriesName,"exp_data_harmonized.csv", sep="_")) 
# read.table("test.txt",header=TRUE,row.names=1) # says first column are rownames