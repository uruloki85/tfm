library(cBioPortalData)
library(AnVIL)

cbio <- cBioPortal()

study <- cBioDataPack(
  "paad_qcmg_uq_2016",
  use_cache = TRUE,
  names.field = c("Hugo_Symbol", "Entrez_Gene_Id", "Gene"),
  ask = TRUE
)

# Metadata
names(study@colData@listData)
# [1] "PATIENT_ID"                  "SAMPLE_ID"                   "SAMPLE_TYPE"                 "SAMPLE_CLASS"               
# [5] "TUMOR_GRADE"                 "ONCOTREE_CODE"               "CANCER_TYPE"                 "CANCER_TYPE_DETAILED"       
# [9] "SOMATIC_STATUS"              "ETHNICITY"                   "SEX"                         "AGE"                        
# [13] "COUNTRY"                     "SMOKER"                      "HISTOLOGICAL_SUBTYPE"        "LOCATION"                   
# [17] "AJCC_PATHOLOGIC_TUMOR_STAGE" "STATUS"                      "DAYS_TO_LAST_FOLLOWUP"      

# Get expression data
eset <- study@ExperimentList@listData[["RNA_Seq_v2_mRNA_median_all_sample_Zscores"]]@assays@data@listData[[1]]

