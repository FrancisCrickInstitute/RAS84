uppsala_file <- file.path( adat$analysis.path, "se_GSE81089.rds" )
if( file.exists( uppsala_file ) ) {
    se_uppsala <- readRDS( uppsala_file )
} else {
    library("biomaRt")
    ensembl <- useMart( "ensembl", dataset="hsapiens_gene_ensembl" )
    gse <- getGEO( 'GSE81089',
                  destdir = adat$analysis.path,
                  GSEMatrix = TRUE )
    eset <- gse[[ 1 ]]
    ## remove normals
    eset <- eset[, pData( eset )$source_name_ch1 == "Human NSCLC tissue" ]
    exprs_dat <- read.delim( file = file.path( "validation_data", "uppsala", "GSE81089_readcounts_featurecou
nts.tsv" ),
                             row.names = 1 )
    ## two sample ids do not match between expression matrix and pData(eset)$title
    ## We remove those here.
    common_sample_ids <- intersect( sampleNames( eset ), colnames( exprs_dat ) )

    eset <- eset[, common_sample_ids ]
    exprs_mat <- as.matrix( exprs_dat[, common_sample_ids ] )
    ensembl_dat <- getBM( attributes = c( "ensembl_gene_id", "external_gene_name" ),
                         filters = "ensembl_gene_id",
                         values = rownames( exprs_mat ),
                         mart = ensembl )
    rowdat <- data.frame( ensembl_gene_id = rownames( exprs_mat ) ) %>%
        left_join( ensembl_dat, by = "ensembl_gene_id" ) %>%
        mutate( rownames = ensembl_gene_id ) %>%
        column_to_rownames( var = "rownames" )
    se <- SummarizedExperiment( assay = exprs_mat,
                               colData = pData( eset ),
                               rowData = rowdat )
    se <- se[ !is.na( rowData( se )$external_gene_name ), ]
    assays( se )$vst <- DESeq2::vst( assay( se ) )
    assays( se )$z <- t( scale( t( scale( assays( se )$vst ) ) ) )
    saveRDS( se, file = file.path( adat$analysis.path, "se_GSE81089.rds" ) )
}


