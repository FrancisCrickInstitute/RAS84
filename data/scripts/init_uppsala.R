
library( DESeq2 )
library( biomaRt )
library( GEOquery )

DOWNLOAD_PATH <- file.path( "data", "downloads" )
GSE81089_DATA_FILE <- file.path( DOWNLOAD_PATH, "GSE81089_readcounts_featurecounts.tsv" )
UPPSALA_SE_FILE <- file.path( "data", "objects", "se_GSE81089.rds" )

ensembl <- useMart( "ensembl", dataset="hsapiens_gene_ensembl" )
gse <- getGEO( 'GSE81089',
               destdir = DOWNLOAD_PATH,
               GSEMatrix = TRUE )
eset <- gse[[ 1 ]]
sampleNames( eset ) <- as.character( eset$title )

## remove normals
eset <- eset[, pData( eset )$source_name_ch1 == "Human NSCLC tissue" ]
exprs_dat <- read.delim( file = GSE81089_DATA_FILE, row.names = 1 )

## two sample ids do not match between expression matrix and pData(eset)$title
## We remove those here.
common_sample_ids <- intersect( sampleNames( eset ), colnames( exprs_dat ) )

eset <- eset[, common_sample_ids ]
exprs_mat <- as.matrix( exprs_dat[, common_sample_ids ] )
ensembl_dat <- getBM( attributes = c( "ensembl_gene_id", "external_gene_name" ),
                      filters = "ensembl_gene_id",
                      values = rownames( exprs_mat ),
                      mart = ensembl,
		      useCache = FALSE )
rowdat <- data.frame( ensembl_gene_id = rownames( exprs_mat ) ) %>%
        left_join( ensembl_dat, by = "ensembl_gene_id" ) %>%
        mutate( rownames = ensembl_gene_id ) %>%
        column_to_rownames( var = "rownames" )
se <- SummarizedExperiment( assay = exprs_mat,
                               colData = pData( eset ),
                               rowData = rowdat )
se <- se[ !is.na( rowData( se )$external_gene_name ), ]
assays( se )$vst <- DESeq2::vst( assay( se ) )
saveRDS( se, UPPSALA_SE_FILE )
