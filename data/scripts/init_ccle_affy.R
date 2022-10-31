
library( affy )

CCLE_EXPRESSION_FILE <- file.path( "data", "downloads", "CCLE_Expression_Entrez_2012-09-29.gct" )
CCLE_ONCOMAP_MAF_FILE <- file.path( "data", "downloads", "CCLE_Oncomap3_2012-04-09.maf" )
ESET_CCLE_AFFY_LUNG_FILE <- file.path( "data", "objects", "eset_ccle_affy_lung.rds" )

##---------------------
## read expression data
exprsDat <-  read.table( file = CCLE_EXPRESSION_FILE,                          
	sep = "\t", header = TRUE, skip = 2 )
fdat <- exprsDat[, 1:2 ]
exprsDat <- exprsDat[, -1:-2 ]
  
##-----------------------
## read mutation maf file
mutdat <- read.table( file = CCLE_ONCOMAP_MAF_FILE,
                      sep = "\t", header = TRUE, quote = "" )
  
##------------------------
## create pdat data.frame
pdat <- data.frame(row.names = colnames( exprsDat ),
                   stringsAsFactors = FALSE )
## Add tumour type column to pdat
pdat$tumour <- factor( sub( "^[A-Z0-9]+_", "", rownames( pdat ), perl = TRUE ) )
  
## Add mutation calls to mutations slot in pdat
mutmat <- data.frame( row.names = row.names( pdat ) )
for( i in levels( mutdat$Hugo_Symbol ) ) {
    ## extract KRAS mutations
    mutdatsub <- subset( mutdat, Hugo_Symbol == i )
    mutmat[[ i ]] <- "no"
    mutmat[ unique( as.character( mutdatsub$Tumor_Sample_Barcode ) ), i ] <- "yes" 
}
mutmat[ is.na( mutmat ) ] <- "no"
mutmat <- mutmat[ rownames( pdat ), ]
mutL <- as.list( as.data.frame( t( mutmat ), stringsAsFactors = FALSE ) )
pdat$mutations <- lapply( mutL[ rownames( pdat ) ], function( x ) {
    names( x ) <- levels( mutdat$Hugo_Symbol )
    x
} )
  
##---------------
## ExpressionSet
eSet <- ExpressionSet( 	as.matrix( exprsDat ),
                       	phenoData = new( 'AnnotatedDataFrame', 
			data = pdat[ colnames( exprsDat ),, drop = FALSE ] ),
                        featureData = new( 'AnnotatedDataFrame', data = fdat ) )
    
eSet <- eSet[, pData( eSet )$tumour %in% "LUNG" ]
pData( eSet ) <- droplevels( pData( eSet ) )
saveRDS( eSet, ESET_CCLE_AFFY_LUNG_FILE )
