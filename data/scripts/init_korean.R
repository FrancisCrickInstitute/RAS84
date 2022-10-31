library( biomaRt )
library( affy )
library( DESeq2 )
library( tidyverse )
library( GEOquery )
library( DT )
library( ggrepel )
library( openxlsx )

KOREAN_EXPRESSION_FILE <- file.path( "..", "objects", "eSet_korean.rds"
KOREAN_EXPRESSION_DATA_FILE <- file.path( "..", "downloads", "GSE40419_LC-87_RPKM_expression.txt" )
KOREAN_MUTATION_DATA_FILE <- file.path( "..", "downloads", "SuppTable3.xlsx" )

eset_gse40419_l <- getGEO( 'GSE40419', GSEMatrix = TRUE )
eset_geo <- eset_gse40419_l[[ 1 ]]
eset_geo$title <- as.character( eset_geo$title )
mut_dat <- read.xlsx( KOREAN_MUTATION_DATA_FILE,
                     sheet = 3,
                     startRow = 2 ) %>%
    filter( Gene %in% c( "TP53", "EGFR", "KRAS" ) ) %>%
    select( -2 ) %>%
    gather( sample_id, value, -Gene ) %>%
    mutate( value = ifelse( value == "0|0", "0", "1" ) ) %>%
    spread( Gene, value )
pdat_geo <- pData( eset_geo ) %>%
    mutate( rownames = title ) %>%
    rename( sample_id = title ) %>%
    left_join( mut_dat, by = "sample_id" ) %>%
    column_to_rownames( var = "rownames" )
rpkm_dat_geo <- read.delim( file = KOREAN_EXPRESSION_DATA_FILE,
                           stringsAsFactors = FALSE )
exprs_mat_geo <- rpkm_dat_geo %>%
    select( -accession, -chrom, -start, -end, -strand, -gene ) %>%
    as.matrix( )
fdat_geo <- rpkm_dat_geo %>%
    rownames_to_column( var = "ID" ) %>%
    select( ID, gene, accession, chrom, start, end, strand ) 
eset_gse40419 <- ExpressionSet( assayData = log2( exprs_mat_geo[, pdat_geo$sample_id] + 0.1 ),
                               phenoData = new( "AnnotatedDataFrame", data = pdat_geo ),
                               featureData = new( "AnnotatedDataFrame", data = fdat_geo ) )
eset_gse40419 <- eset_gse40419[, eset_gse40419$source_name_ch1 == "Lung cancer cells" ]
saveRDS( eset_gse40419, KOREAN_EXPRESSION_FILE )
