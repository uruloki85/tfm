library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

filters = listFilters(ensembl)

probeids=c('1007_s_at','1053_at','1053_at','117_at','117_at','121_at','121_at','1255_g_at','1255_g_at')

# probeids=c('200007_at', '200011_s_at', '200012_x_at')
bm <- getBM(attributes=c('affy_hg_u133_plus_2', 'hgnc_symbol'), 
      filters = 'affy_hg_u133_plus_2', 
      values = probeids, 
      mart = ensembl)
