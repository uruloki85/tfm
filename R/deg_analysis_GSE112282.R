library(dplyr)
library(readr)
require(GEOquery)
library(limma)
library(ggplot2)

seriesName <- "GSE112282"

gse <- getGEO(seriesName, GSEMatrix=TRUE, getGPL = TRUE)
gse <- gse[[1]]
# show(gse)

## exprs get the expression levels as a data frame and get the distribution
summary(exprs(gse))
# if data is log2, will be between 0 and 16

# Verify data has been normalized
exprs(gse) <- log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)

# Get sample info
sampleInfo <- pData(gse)
sampleInfo

# Select some columns
if(seriesName == "GSE112282") {
  sampleInfo <- select(sampleInfo, "cell line:ch1", "replicate info:ch1", "treatment:ch1")
  sampleInfo <- rename(sampleInfo, line="cell line:ch1", replicate="replicate info:ch1", treatment="treatment:ch1")
}

features <- fData(gse)
# View(features)
# We keep the probe ID and the gene accession
if(seriesName == "GSE112282") {
  features <- select(features,ID,GB_ACC)
}
full_output <- cbind(features,exprs(gse))
# write_csv(full_output,file=paste(seriesName,"exp_data_harmonized.csv", sep="_")) 

# design <- model.matrix(~0+sampleInfo$treatment)
# colnames(design) <- c("BET","BETMEK","MEK","VEHICLE")

if(seriesName == "GSE112282") {
  design_colnames <- c("BET","BETMEK","MEK","VEHICLE","COLO201","HPAFII","NCIH510","RKO","Replicate2")
  design <- model.matrix(~0+sampleInfo$treatment+sampleInfo$line+sampleInfo$replicate)
  colnames(design) <- design_colnames
  contrasts <- makeContrasts(BET - VEHICLE, BETMEK - VEHICLE, MEK - VEHICLE, levels=design)
}
design
design_colnames


fit <- lmFit(exprs(gse), design)
head(fit$coefficients)

contrasts

fit2 <- contrasts.fit(fit, contrasts)

# Get differential expression statistics and p-values with empirical Bayes
fit2 <- eBayes(fit2)

# Results by contrast
topTable(fit2)
topTable(fit2, coef=1)
topTable(fit2, coef=2)
topTable(fit2, coef=3)
 # topTable(fit2, coef=4)

# How many genes are differentially-expressesd
results <- decideTests(fit2)
table(results)

vennDiagram(results)


# Print top 10 DEG with meaningful name
gene_accession <- "GB_ACC"

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