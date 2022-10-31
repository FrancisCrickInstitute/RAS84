
library( tidyverse )
library( SummarizedExperiment )
library( openxlsx )

GDSC_SE_L_FILE <- file.path( "data", "objects", "gdsc_se_l.rds" )
GDSC1_DATA_FILE <- file.path( "data", "downloads", "GDSC1_fitted_dose_response_25Feb20.xlsx" )
GDSC2_DATA_FILE <- file.path( "data", "downloads", "GDSC2_fitted_dose_response_25Feb20.xlsx" )

GDSC1_dat <- read.xlsx( GDSC1_DATA_FILE )
GDSC2_dat <- read.xlsx( GDSC2_DATA_FILE )

gdsc_dat_l <- list( GDSC1 = GDSC1_dat, GDSC2 = GDSC2_dat ) %>%
    map( group_by, DRUG_NAME, CELL_LINE_NAME ) %>%
    map( filter, LN_IC50 == min( LN_IC50 ) ) %>%
    map( ungroup )
gdsc_ic50_l <- gdsc_dat_l %>%
    map( dplyr::select, DRUG_NAME, CELL_LINE_NAME, LN_IC50 ) %>%
    map( spread, CELL_LINE_NAME, LN_IC50 ) %>%
    map( column_to_rownames, var = "DRUG_NAME" ) %>%
    map( as.matrix )
gdsc_auc_l <- gdsc_dat_l %>%
    map( dplyr::select, DRUG_NAME, CELL_LINE_NAME, AUC ) %>%
    map( spread, CELL_LINE_NAME, AUC ) %>%
    map( column_to_rownames, var = "DRUG_NAME" ) %>%
    map( as.matrix )
gdsc_coldat_l <- gdsc_dat_l %>%
    map( dplyr::select, CELL_LINE_NAME, TCGA_DESC ) %>%
    map( distinct ) %>%
    map( mutate, CELL_LINE_ID = gsub( "-", "", CELL_LINE_NAME ) ) %>%
    map( mutate, rownames = CELL_LINE_NAME ) %>%
    map( column_to_rownames, var = "rownames" )
gdsc_rowdat_l <- gdsc_dat_l %>%
    map( dplyr::select, DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME ) %>%
    map( mutate, TARGET_CATEGORY = PATHWAY_NAME ) %>%
    map( distinct ) %>%
    ## remove drug annotation errors
    map( function( df ) {
        df %>%
            filter( !grepl( "FGRF1", PUTATIVE_TARGET ) ) %>%
            filter( !DRUG_NAME %in% c( "Dactinomycin", "Fulvestrant" ) & !is.na( PUTATIVE_TARGET ) )
    } ) %>%
    map( mutate, rownames = DRUG_NAME ) %>%
    map( column_to_rownames, var = "rownames" )
gdsc_se_l <- map2( gdsc_ic50_l, names( gdsc_ic50_l ), function( ic50_mat, n ) {
        SummarizedExperiment( assays = list(
                              response = ic50_mat,
                              ln_IC50 = ic50_mat,
                              AUC = gdsc_auc_l[[ n ]] ),
                              colData = gdsc_coldat_l[[ n ]][ colnames( ic50_mat ), ],
                              rowData = gdsc_rowdat_l[[ n ]][ rownames( ic50_mat ), ] )
} )

saveRDS( gdsc_se_l, GDSC_SE_L_FILE )
