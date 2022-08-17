
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
se_t_file <- file.path( "analysis", "sample_classification", "TCGA", tumour, "se_t_RAS84_rasclass.rda" )
load( file = se_t_file )
assay( se_t ) <- round( assay( se_t ) )
assays( se_t )$z <- t( scale( t( scale( assays( se_t )$vst ) ) ) )

## load signatures
load( file = file.path( "analysis", "validate_signatures", "TCGA", "LUAD", "signatures.rda" ) )

sig_f <- rowData( se_t )$gene_id %in% signatures$RAS84

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
train_test_dat_l_file <- file.path( adat$analysis.path, "test_train_dat_l.rda" )
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
RF_rank_file <- file.path( adat$analysis.path, "RF_ranked_RAS84_genes.txt" )
if( !file.exists( RF_rank_file ) ) {
    data.frame( genes = predictors( res ),
               RF_rank = 1:length( predictors( res ) ) ) %>%
        mutate( genes = sub( "\\..*", "", genes ) ) %>%
        write.table( file = ,
                    sep = "\t", quote = FALSE, row.names = FALSE )
}

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

## lowest gene rank > 0.9 model performance
## 34 genes
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

## Genes to achieve >90% specificity
byClass_df %>%
    filter( metric == "Specificity" ) %>%
    group_by( rfgene_rank ) %>%
    dplyr::summarize( threshold = all( value > 0.9 ) )
## 11 genes

byClass_df %>%
    filter( metric == "Sensitivity" ) %>%
    group_by( rfgene_rank ) %>%
    dplyr::summarize( threshold = all( value > 0.8 ) )
## 14 genes

## Top 34 genes
RAS34 <- predictors( res )[ 1:34 ]
RAS34 <- sub( "\\..*", "", RAS34 )
saveRDS( svm_models_l[[ "34" ]], file.path( adat$analysis.path, "RAS34_svm_model.rds" ) )
write.table( data.frame( C1 = RAS34 ), file = file.path( "data", "signatures", "ras", "active", "RAS34.csv" ),
             col.names = FALSE, row.names = FALSE, quote = FALSE )

models_file <- file.path( adat$analysis.path, "models.rda" )
if( file.exists( models_file ) ) {
    load( file = models_file )
} else {
    ##-------------------------
    ## - filter datasets
    for( gset in names( rfe_genes ) ) {
        
        train_dat <- train_test_dat_l$train[, c( "condition", rfe_genes[[ gset ]] ) ]
        test_dat <- train_test_dat_l$test[, c( "condition", rfe_genes[[ gset ]] ) ]
        
        ##--------------------------------
        ## Random Forest
        ##model_name <- paste( "rf", gset, sep = "." )
        ##modelLookup( "rf" )
        ##set.seed( 100 )
        ## - build model
        ##model_rf <- train( condition ~ .,
        ##                  data = train_dat,
        ##                  method = "rf",
        ##                  trControl = tr_control )
    
        ## - variable importance
        ##var_imp_rf <- varImp( model_rf, scale = FALSE )
        
        ## - test
        ##test_res_rf <- predict( model_rf, test_dat[, -1 ] )
        
        ## - confusion matrix
        ##cm_rf <- confusionMatrix( reference = test_dat$condition,
        ##                         data = test_res_rf,
        ##                         mode = 'everything',
        ##                         positive = 'MM' )
        ##
        ##models[[ model_name ]] <- list( model = model_rf,
        ##                               cm_test = cm_rf,
        ##                               varimp = var_imp_rf )
    
        ##---------------------------------------------------------
        ## svmRadial
        model_name <- paste( "svm", gset, sep = "." )
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
        
        models[[ model_name ]] <- list( model = model_svm,
                                       cm_test = cm_svm,
                                       varimp = var_imp_svm )
    }
    save( models, file = models_file )
}

##------------------------------------------------
## Create RAS33 SVM model

##-------------------------------------------------
## Compare models
    
models_compare <- resamples( lapply( models, function( l ) l$model ) )
summary( models_compare )
bwplot( models_compare )

##-------------------------------
## plot metrics
plotdat <- melt( lapply( models, function( l ) l$cm$byClass[, c( "Sensitivity", "Specificity" ) ] ) )

gg <- ggplot( plotdat, aes( x = Var1, y = value, color = L1, group = L1 ) ) +
    ##geom_bar( stat = "identity", position = "dodge" ) +
    geom_point( ) +
    geom_line( ) +
    facet_wrap( ~Var2 ) +
    theme( axis.text.x = element_text( angle = 90, hjust = 1 ) )
gg

## - overall accuracy
lapply( models, function( l ) l$cm_test$overall[ "Accuracy" ] )

## - write short signature
RAS47_genes <- sapply( strsplit( rfe_genes$rf, "\\." ), '[[', 1 )

write.table( RAS47_genes, file = file.path( adat$analysis.path, "RAS47.csv" ),
            row.names = FALSE, quote = FALSE, col.names = FALSE )

##-----------------------------------------------------
## Adaboost

modelLookup( "adaboost" )
set.seed( 100 )

model_ada_file <- file.path( adat$analysis.path, "model_ada.rda" )
if( file.exists( model_ada_file ) ) {
    load( file = model_ada_file )
} else {
    model_ada <- train( condition ~ .,
                       data = train_dat,
                       method = "adaboost",
                       tuneLength = 2,
                       trControl = tr_control )
    save( model_ada, file = model_ada_file )
}
plot( model_ada )

fitted_ada <- predict( model_ada )
## - accuracy validation
sum( fitted_ada == train_dat$condition ) / length( fitted_ada )

## - variable importance
var_imp_ada <- varImp( model_ada, scale = FALSE )

## - test
test_res_ada <- predict( model_ada, test_dat[, -1 ] )

## - accuracy 0.85
sum( test_res_ada == test_dat$condition ) / length( test_res_ada )

## - confusion matrix
cm <- confusionMatrix( reference = test_dat$condition,
                      data = test_res_ada,
                      mode = 'everything',
                      positive = 'MM' )

models$ada <- list( model = model_ada,
                   cm = cm,
                   varimp = var_imp_ada )

##---------------------------------------------------------
## xgBoost Dart
modelLookup( "xgbDART" )
set.seed( 100 )

model_xgb_file <- file.path( adat$analysis.path, "model_xgb.rda" )
if( file.exists( model_xgb_file ) ) {
    load( file = model_xgb_file )
} else {
    model_xgb <- train( condition ~ .,
                       data = train_dat,
                       method = "xgbDART",
                       tuneLength = 5,
                       trControl = tr_control )
    save( model_xgb, file = model_xgb_file )
}

plot( model_xgb )

fitted_xgb <- predict( model_xgb )
## - accuracy validation
sum( fitted_xgb == train_dat$condition ) / length( fitted_xgb )

## - variable importance
var_imp_xgb <- varImp( model_xgb, scale = FALSE )

## - test
test_res_xgb <- predict( model_xgb, test_dat[, -1 ] )

## - accuracy 0.85
sum( test_res_xgb == test_dat$condition ) / length( test_res_xgb )

## - confusion matrix
cm <- confusionMatrix( reference = test_dat$condition,
                      data = test_res_xgb,,
                      mode = 'everything',
                      positive = 'MM' )

models$xgb <- list( model = model_xgb,
                   cm = cm,
                   varimp = var_imp_xgb )

##-----------------------------------------
## classification using RAS47 via SVM

## - filter se_t for RAS47
sig_n <- "RAS47"
sig_f <- rowData( se_t )$gene_id %in% signatures[[ sig_n ]]
se_t_ras47 <- se_t[ sig_f, ]

## build SVM model
model_svm_ras47_file <- file.path( adat$analysis.path, "model_svm_ras47.rda" )
if( file.exists( model_svm_ras47_file ) ) {
    load( file = model_svm_ras47_file )
} else {
    ## - filter se
    sig_n <- "RAS47"
    sig_f <- rowData( se_t )$gene_id %in% signatures[[ sig_n ]]
    se_t_ras47 <- se_t[ sig_f, ]
    ## - create input data for RAS47
    model_dat <-  data.frame(
        condition = se_t_ras47$rasact_5,
        t( assays( se_t_ras47 )$vst ) )
    ## - build model
    set.seed( 100 )
    tr_control <- trainControl( method = "cv", number = 10,
                               classProbs = TRUE,
                               savePredictions = TRUE )
    model_svm_ras47 <- train( condition ~ .,
                             data = model_dat,
                             method = "svmRadial",
                             tuneLength = 15,
                             preProc = c( "center", "scale"), 
                             trControl = tr_control )
    save( model_svm, file = model_svm_ras47_file )
}

library(pROC)
f <- 
png( file = file.path( adat$plot.path, "ROC.png" ) )
plot.roc( model_svm_ras47$pred$obs, model_svm_ras47$pred$RASsig_max )
dev.off()

preds_ras47 <- predict( model_svm_ras47 )

## confusion matrix
cm_ras74 <- confusionMatrix( preds_ras47, se_t_ras47$rasact )

## - save RAS47 classifications to se_t
se_t$rasact_ras47 <- preds_ras47

save( se_t, file = file.path( adat$analysis.path, "se_t_ras47.rda" ) )

##------------------------------------------------
## - Short signature classification visualisation

hmcols <- colorpanel( 100, my_blue, "white", my_red )

## - All the samples
hm_anno_cols <- list( rasact_ras47 = c( RASsig_max = my_red, RASsig_3 = my_orange,
                                        RASsig_2 = my_green, RASsig_1 = "lightblue", RASsig_0 = my_blue ) )
## - create heatmap annotation object - need to add mutational load to this
hm_anno <- rowAnnotation( df = as.data.frame( colData( se_t ) )[, c( "rasact_ras47" ), drop = FALSE ],
                         col = hm_anno_cols,
                         annotation_width = 0.3 )

## - RAS47
sig_f <- rowData( se_t )$gene_id %in% signatures$RAS84

hmdat <- assays( se_t )$vst[ sig_f, ]
hmdat <- hmdat - rowMedians( hmdat )

hm_RAS47 <- Heatmap( t( hmdat ),
                    column_title = "RAS84 with RAS47 SVM classifications",
                    clustering_method_rows = "ward.D2",
                    clustering_method_columns = "ward.D2",
                    show_column_names = TRUE,
                    show_row_names = FALSE,
                    col = colorRamp2(
                        seq( from = -3, to = 3, length.out = length( hmcols ) ),
                        hmcols ),
                    split = 5,
                    column_names_gp = gpar( fontsize = 8 ) ) + hm_anno

pdf( file = file.path( adat$plot.path, "RAS84_RAS47_class_fullset_heatmaps.pdf" ),
     width = 10 )
hm_RAS47
dev.off( )

## - test samples
test_dat <- train_test_dat_l$test[, c( "condition", rfe_genes$rf ) ]

se_t_test <- se_t[, rownames( test_dat ) ]

preds_ras47_svm_test <- predict( models$svm.rf$model, test_dat[, -1 ] )

## - confusion matrix
cm_test <- confusionMatrix( test_dat$condition, preds_ras47_svm_test )

hmdat <- assays( se_t_test )$vst[ sig_f, ]
hmdat <- hmdat - rowMedians( hmdat )

hm_anno <- rowAnnotation( df = data.frame( rasact_ras47 = preds_ras47_svm_test ),
                          col = hm_anno_cols,
                          annotation_width = 0.3 )

hm_RAS47_testset <- Heatmap( t( hmdat ),
                            column_title = paste( "RAS84 with RAS47 SVM classes (20% test data)",
                                                   round( cm_test$overall[ 1 ], 2 ), "accuracy" ),
                            clustering_method_rows = "ward.D2",
                            clustering_method_columns = "ward.D2",
                            show_column_names = TRUE,
                            show_row_names = FALSE,
                            col = colorRamp2(
                                seq( from = -3, to = 3, length.out = length( hmcols ) ),
                                hmcols ),
                            split = se_t_test$rasact,
                            column_names_gp = gpar( fontsize = 8 ) ) + hm_anno

pdf( file = file.path( adat$plot.path, "RAS84_RAS47_class_testset_heatmaps.pdf" ),
    width = 10 )
hm_RAS47_testset
dev.off( )
