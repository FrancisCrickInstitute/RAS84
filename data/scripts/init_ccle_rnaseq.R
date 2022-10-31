library( DESeq2 )
library( tidyverse )
library( biomaRt )
library( openxlsx )

CCLE_RNASEQ_DAT_FILE <- file.path( "..", "downloads", "CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz" )
CCLE_MUTATIONS_FILE <- file.path( "..", "downloads", "CCLE_DepMap_18q3_maf_20180718.txt" )
CCLE_RNASEQ_SE_FILE <- file.path( "..", "objects", "ccle_rnaseq_se.rda" ) 
RAS_PATHWAY_DEF_FILE <- file.path( "..", "downloads", "1-s2.0-S0092867418303593-mmc3.xlsx" )
CELL_DRIVER_MUTATION_FILE <- file.path( "..", "downloads", "1-s2.0-S009286741830237X-mmc1.xlsx" )

rasPathwayTab <- read.xlsx( RAS_PATHWAY_DEF_FILE, sheet = 8 )
rasGeneSets <- list( og = subset( rasPathwayTab, OG.TSG == "OG" )$Gene,
                    tsg = subset( rasPathwayTab, OG.TSG == "TSG" )$Gene )
rasGeneSets <- lapply( rasGeneSets, as.character )
ras_pathway_genes <- unlist( rasGeneSets )

## driver mutations from Cell paper
driver_muts_dat <- read.xlsx( CELL_DRIVER_MUTATION_FILE, sheet = 5, startRow = 2 ) %>%
    unite( driver_gene_index, Gene, Mutation, sep = "_", remove = FALSE ) %>%
    filter( Gene %in% ras_pathway_genes )

## Readin RNA-seq data
rnaseq_dat <- read.delim( CCLE_RNASEQ_DAT_FILE ) %>%
	separate( gene_id, c( "ensembl_gene_id", "ensembl_id_version" ), "\\.", remove = FALSE ) %>%
	dplyr::select( -transcript_ids )

## Gene annotation
row_dat <- rnaseq_dat %>%
        dplyr::select( gene_id, ensembl_gene_id, ensembl_id_version )
ensembl <- useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
ensembl_dat <- getBM( attributes = c( "ensembl_gene_id", "external_gene_name", "entrezgene_id" ),
                      filters = "ensembl_gene_id",
                      values = row_dat$ensembl_gene_id,
                      mart = ensembl )
ensembl_dat <- ensembl_dat %>%
        filter( !duplicated( ensembl_gene_id ) )
row_dat <- row_dat %>%
        left_join( ensembl_dat, by = "ensembl_gene_id" ) %>%
        mutate( entrezgene_id = as.character( entrezgene_id ) ) %>%
        mutate( rownames = gene_id ) %>%
        column_to_rownames( var = "rownames" )

## Count data
assay_dat <- rnaseq_dat %>%
        dplyr::select( -ensembl_gene_id, -ensembl_id_version ) %>%
        column_to_rownames( var = "gene_id" ) %>%
        as.matrix( )

## Mutation data
mutation_RAS_pway_dat <- read.delim( CCLE_MUTATIONS_FILE ) %>%
        filter( Variant_Classification != "Silent" ) %>%
        filter( Protein_Change != "" ) %>%
        dplyr::select( Hugo_Symbol, Tumor_Sample_Barcode, Protein_Change ) %>%
        unite( cc_gene_index, Tumor_Sample_Barcode, Hugo_Symbol, sep = "_", remove = FALSE ) %>%
        unite( driver_gene_index, Hugo_Symbol, Protein_Change, sep = "_", remove = FALSE ) %>%
        filter( Hugo_Symbol %in% ras_pathway_genes )
mutation_RAS_pway_dat_l <- split( mutation_RAS_pway_dat, mutation_RAS_pway_dat$cc_gene_index )
aa_changes_pergene <- map( mutation_RAS_pway_dat_l, function( df ) {
        data.frame( merged_protein_change = paste( df$Protein_Change, collapse = ";" ) )
    } ) %>% bind_rows( .id = "cc_gene_index" )
mutation_RAS_pway_drivers_df <- mutation_RAS_pway_dat %>%
        filter( driver_gene_index %in% driver_muts_dat$driver_gene_index ) %>%
        left_join( aa_changes_pergene, by = "cc_gene_index" ) %>%
        dplyr::select( -Protein_Change, -cc_gene_index, -driver_gene_index ) %>%
        distinct( ) %>%
        group_by( Tumor_Sample_Barcode ) %>%
        mutate( RAS_pathway = any( !is.na( merged_protein_change ) ) ) %>%
        ungroup( ) %>%
        spread( Hugo_Symbol, merged_protein_change )
KRAS_mutants_dat <- mutation_RAS_pway_dat %>%
        filter( Hugo_Symbol == "KRAS" )

## Cell line data data
col_dat <- data.frame( Tumor_Sample_Barcode = colnames( assay_dat ) ) %>%
        separate( Tumor_Sample_Barcode, c( "cell_line_id", "cancer_type" ), "_",
                 extra = "merge", remove = FALSE ) %>%
        mutate( CELL_LINE_ID = cell_line_id ) %>%
        mutate( CELL_LINE_NAME = cell_line_id ) %>%
        left_join( mutation_RAS_pway_drivers_df, by = "Tumor_Sample_Barcode" ) %>%
        mutate( RAS_pathway = ifelse( is.na( RAS_pathway ), FALSE, TRUE ) ) %>%
        mutate( KRASmut_driver = ifelse( !is.na( KRAS ), "mutant", "wt" ) ) %>%
        mutate( KRASmut_driver = factor( KRASmut_driver, levels = c( "wt", "mutant" ) ) ) %>%
        mutate( KRASmut = ifelse( Tumor_Sample_Barcode %in% KRAS_mutants_dat$Tumor_Sample_Barcode,
                                 "mutant", "wt" ) ) %>%
        mutate( KRASmut = factor( KRASmut, levels = c( "wt", "mutant" ) ) ) %>%
        mutate( rownames = Tumor_Sample_Barcode  ) %>%
        column_to_rownames( var = "rownames" )   

## SummarizedExperiment object
rnaseq_se <- SummarizedExperiment( assays = list( counts = round( assay_dat ) ),
                                   colData = col_dat[ colnames( assay_dat ), ],
                                   rowData = row_dat[ rownames( assay_dat ), ] )
## VST normalise
assays( rnaseq_se )$vst <- vst( assay( rnaseq_se ) )
saveRDS( rnaseq_se, CCLE_RNASEQ_SE_FILE )
