
library( affy )
library( limma )
library( org.Hs.eg.db )
library( ggplot2 )
library( reshape2 )
library( scales )
library( RColorBrewer )
library( ComplexHeatmap )
library( circlize )
library( colorspace )
library( GetoptLong )
library( dendextend )
library( matrixStats )
library( DESeq2 )
library( gplots )
library( tsne  )
library( intPredict )
library( pamr )
library( randomForest )
library( tibble )
library( dplyr )
library( plyr )
library( caret )
library( pROC )
library( ROCR )
library( grid )
library( gridExtra )
library( fastAdaboost )
library( tidyverse )

source( "/camp/stp/babs/working/eastp/code/R/R.packages.east01.functions.R")
source( "scripts/lib.Ras_HighLow.335.R" )
source( "scripts/lib.CCLE.335.R" )

my_grey <- "#707173"
my_red <- "#e3001a"
my_orange <- "#f6ad6e"
my_green <- "#7ab51d"
my_lightgreen <- "#adcf82"
my_purple <- "#bb90bd"
my_blue <- "#4066aa"

plot_formatter <- function() {
    theme_bw( ) +
        theme( panel.grid.major = element_blank( ),
               panel.grid.minor = element_blank( ),
               panel.border = element_blank( ),
               panel.background = element_blank( ),
               axis.line = element_line(color = "black"),
               axis.line.x = element_line(color="black", size = 0.1 ),
               axis.line.y = element_line(color="black", size = 0.1 ),
               text = element_text( size = 12 ) )
}

tumour <- "LUAD"
adat <- project.init( file.path( "RAS_short_signature", "TCGA", tumour ) )

## - load classified tumour samples
se_t_file <- file.path( "data", "LUAD_se_t_rasclass.rds" )
se_t <- readRDS( se_t_file )

ras84_map <- read.xlsx( file.path( "data", "RAS84_feature_map.xlsx" ) )
sig_f <- rowData( se_t )$gene_id %in% ras84_map$TCGA_feature_id
se_t_ras84 <- se_t[ sig_f, ]

##----------------------------
## modeling data
model_dat <-  data.frame(
    condition = se_t_ras84$RAS84,
    t( assays( se_t_ras84 )$z ) )

##-----------------------------------
## - QC check for near zero variance
nzv <- nearZeroVar( model_dat, saveMetrics= TRUE)

##--------------------------------------
## - identifying correlated predictors
cor_mat <- cor( model_dat[, -1 ] )

pw_ras84_gene_cor_gg <- data.frame( pw_cor = as.numeric( cor_mat[ upper.tri( cor_mat ) ] ) ) %>%
    ggplot( aes( x = pw_cor ) ) +
    geom_density( ) +
    labs( x = "RAS84 gene Pearson coefficients" ) +
    plot_formatter( )
cairo_pdf( file = file.path( adat$plot.path, "pw_ras84_gene_cor_gg.pdf" ),
          width = 3, height = 3)
print( pw_ras84_gene_cor_gg )
dev.off()

## Max abs correlation coefficient per gene
cor_mat[ cor_mat == 1 ] <- 0
apply( cor_mat, 1, function( x ) any( abs( x ) > 0.5 ) )

max_abs_cor_coef <- apply( cor_mat, 1, function( x ) x[ abs( x ) == max( abs( x ) ) ] )
max_ras84_gene_cor_gg <- data.frame( max_abs_cor_coef = max_abs_cor_coef ) %>%
    ggplot( aes( x = max_abs_cor_coef ) ) +
    geom_density( ) +
    lims( x = c( -1, 1 ) ) +
    labs( x = "Maximum RAS84 gene Pearson coefficient per-gene" ) +
    plot_formatter( )

cairo_pdf( file = file.path( adat$plot.path, "max_ras84_gene_cor_gg.pdf" ),
          width = 3, height = 3)
print( max_ras84_gene_cor_gg )
dev.off()

## - identify highly correlated genes
high_cor <- findCorrelation( cor_mat, cutoff = 0.75, verbose = TRUE )

##---------------------------------
## - Check linear dependencies
## - necessary for multi-factorial designs
findLinearCombos()

##----------------------------------
## - class distance
centroids <- classDist( model_dat[, 2:60 ], model_dat$condition )
distances <- predict(centroids, testBC)
distances <- as.data.frame(distances)
head(distances)

##----------------------------------
## split data train/test
train_test_dat_l_file <- file.path( "data", "test_train_dat_l.rda" )
if( file.exists( train_test_dat_l_file ) ) {
    load( file = train_test_dat_l_file )
} else{
    set.seed( 100 )
    train_index <- createDataPartition( model_dat$condition, p = 0.8, list = FALSE, times = 1 )
    train_test_dat_l <- list( train_dat = model_dat[ train_index, ],
                             test_dat = model_dat[ -train_index, ] )
    save( train_test_dat_l, file = train_test_dat_l_file )
}

##---------------------------------
## - feature selection using rft

## - method: resampling methods
## repeated cross-validation (repeatedcv) is used here. Other options
## include leave-one-out cross-validation and bootstrap
## (simple estimation or the 632 rule)

## - number: Either the number of folds or number of resampling iterations

## - repeats: For repeated k-fold cross-validation only: the number of
## complete sets of folds to compute

rfe_file <- file.path( adat$analysis.path, "rfe_res_l_84.rda" )
if( file.exists( rfe_file ) ) {
    load( file = rfe_file )
} else {
    rfe_funcs_l <- list( rf = rfFuncs, nb = nbFuncs, treebag = treebagFuncs )
    rfe_res_l <- lapply( names( rfe_funcs_l )[ 1 ], function( func_n ) {
        print( func_n )
        set.seed( 100 )
        rfe_control <- rfeControl( functions = rfe_funcs_l[[ func_n ]],
                                  method = "repeatedcv",
                                  number = 10,
                                  repeats = 5 )
        rfe( model_dat[ , -1 ], model_dat$condition,
            sizes = 1:84,
            rfeControl = rfe_control )
    } )
    names( rfe_res_l ) <- names( rfe_funcs_l )[ 1 ]
    save( rfe_res_l, file = rfe_file )
}

res <- rfe_res_l$rf

## genes / predictors
rfe_genes <- predictors( res )

## - plotting
p <- lapply( names( rfe_res_l ), function( l_n ) {
    p <- plot( rfe_res_l[[ l_n ]],
              type = c( "g", "o" ),
              main = l_n )
} )
grid.arrange( p[[1]], p[[2]], p[[3]], ncol = 1 )

##---------------------------------
## - build model
models <- list()

##----------------------------
## define training parameters
tr_control <- trainControl( method = "cv", number = 10  )

svm_rfgenes_file <- file.path( adat$analysis.path, "svm_rfgenes.rda" )
svm_models_l <- lapply( 2:81, function( i ) {
    print( i )
    training_genes <- predictors( res )[ 1:i ]
    train_dat <- train_test_dat_l$train[, c( "condition", training_genes ) ]
    test_dat <- train_test_dat_l$test[, c( "condition", training_genes ) ]
    model_name <- paste( "svm", i, sep = "." )
    modelLookup( "svmRadial" )
    ## - build model
    set.seed( 100 )
    model_svm <- train( condition ~ .,
                       data = train_dat,
                       method = "svmRadial",
                       tuneLength = 15,
                       trControl = tr_control )
    ## - variable importance
    var_imp_svm <- varImp( model_svm, scale = FALSE )
    ## - test
    test_res_svm <- predict( model_svm, test_dat[, -1 ] )
    ## - confusion matrix
    cm_svm <- confusionMatrix( reference = test_dat$condition,
                              data = test_res_svm,,
                              mode = 'everything',
                              positive = 'MM' )    
    list( model = model_svm,
         cm_test = cm_svm,
         varimp = var_imp_svm )
} )

names( svm_models_l ) <- as.character( 2:81 )
saveRDS( svm_models_l, file = svm_rfgenes_file )

svm_models_l <- readRDS( file = svm_rfgenes_file )
svm_models_l <- svm_models_l[ !names( svm_models_l ) %in% c( "1", "2", "3", "4", "5" ) ]

accuracy_df <- map( svm_models_l, function( l ) {
    l$cm_test$overall %>%
        t( ) %>%
        as.data.frame( ) %>%
        dplyr::select( Accuracy, AccuracyLower, AccuracyUpper )
} ) %>%
    bind_rows( .id = "rfgene_rank" ) %>%
    gather( metric, value, -rfgene_rank ) %>%
    mutate( greater_than_0.9 = value > 0.9 )

data <- filter( accuracy_df, metric == "Accuracy" ) %>%
    mutate( rfgene_rank = as.numeric( rfgene_rank ) )
loess_fit <- loess( value ~ rfgene_rank, data = data )
predict_perform <- map( 40:70, function( x ) predict( loess_fit, x ) ) %>%
    unlist( ) %>%
    setNames( 40:70 )
predict_perform[ predict_perform == max( predict_perform ) ]

accuracy_gg <- accuracy_df %>%
    ggplot( aes( x = as.numeric( rfgene_rank ), y = value, group = metric,
                linetype = metric, color = metric ) ) +
    geom_vline( xintercept = 55 ) +
    geom_hline( yintercept = 0.9, linetype = "dotted" ) +
    geom_point( ) +
    geom_line( ) +
    geom_smooth( data = filter( accuracy_df, metric == "Accuracy" ), method = "loess", se = FALSE ) +
    scale_color_manual( values = c( "black", "grey", "grey" ) ) +
    scale_linetype_manual( values = c( "solid", "dotted", "dotted" ) ) +
    labs( x = "RF gene rank count", y = "model performance" ) +
    plot_formatter( )

cairo_pdf( file = file.path( adat$plot.path, "model_performance.pdf" ),
          width = 5.5, height = 3 )
print( accuracy_gg )
dev.off( )

byClass_df <- map( svm_models_l, function( l ) {
    l$cm_test$byClass %>%
        as.data.frame( ) %>%
        dplyr::select( Sensitivity, Specificity ) %>%
        rownames_to_column( var = "RAG" )
} ) %>%
    bind_rows( .id = "rfgene_rank" ) %>%
    gather( metric, value, -rfgene_rank, -RAG ) %>%
    mutate( rfgene_rank = as.numeric( as.character( rfgene_rank ) ) )

byClass_gg <- byClass_df %>%
    ggplot( aes( x = rfgene_rank, y = value, color = RAG, group = RAG ) ) +
    #geom_point( ) +
    geom_line( ) +
    geom_smooth( method = "loess", se = FALSE ) +
    facet_wrap( ~metric, scales = "free_y" ) +
    lims( y = c( 0, 1 ) ) +
    labs( x = "RF gene rank count", y = "model performance" ) +
    plot_formatter( )

data <- filter( byClass_df, metric == "Sensitivity" & RAG == "Class: RAG-0" ) %>%
    mutate( rfgene_rank = as.numeric( rfgene_rank ) )
loess_fit <- loess( value ~ rfgene_rank, data = data )
predict( loess_fit, 43 )

cairo_pdf( file = file.path( adat$plot.path, "model_performance_byClass.pdf" ),
          width = 7, height = 3 )
print( byClass_gg )
dev.off( )
