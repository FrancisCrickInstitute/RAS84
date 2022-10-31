
library( SummarizeExperiment )

## set tumours vector to available tumour data
tumours <- c( "LUAD" )

## Paths
OBJECT_PATH <- file.path( "..", "objects" )
RAS84_FEATURE_MAP_FILE <- file.path( "..", "resources" )
PANCANCER_SE_FILE <- "se_pancancer.rds"
PANCANCER_SE_PATH <- file.path( OBJECT_PATH, PANCANCER_SE_FILE )

## read in feature id map
ras84_map <- read.xlsx( RAS84_FEATURE_MAP_FILE )

## read in individual tumour expression data. The SE objects here will have been downloaded previously using data/scripts/init_TCGA.R
se_l <- sapply( tumours, function( tumour ) {
    se_t_file <- paste0( "se_t_", tumour, ".RNA-Seq.legacy.biolinks.rda" )
    se_t_file <- file.path( OBJECT_PATH, se_t_file )
    load( file = se_t_file )
    assays( se_t )$zscore <- apply( assays( se_t )$vst, 2, function( x ) ( x - median( x ) )/ mad( x ) )
}, simplify = FALSE )

## identify common colData header names across tumour type SE objects
common_colnames <- table( unlist( lapply( se_l, function( se ) {
    colnames( colData( se ) )
} ) ) )
common_colnames <- names( common_colnames[ common_colnames == length( se_l ) ] )
se_pc <- do.call( 'cbind', lapply( se_l, function( se ) {
    colData( se ) <- colData( se )[, common_colnames ]
    se
} ) )
se_pc$tumour <- factor( se_pc$tumour )

## Calculate pancancer RAS index value
sig_f <- rowData( se_pc )$gene_id %in% ras84_map$RAS84
se_pc$ras_index_pc <- colMeans( assays( se_pc )$zscore[ sig_f, ] )

## save pancancer SE object
saveRDS( se_pc, file = PANCANCER_SE_PATH )

