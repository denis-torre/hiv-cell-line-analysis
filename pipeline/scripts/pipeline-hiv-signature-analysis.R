#################################################################
#################################################################
############### HIV Signature Analysis
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####
source('/Users/denis/Documents/Projects/scripts/Support.R')

##### 2. Other libraries #####


#######################################################
#######################################################
########## S1. Data Normalization
#######################################################
#######################################################

#############################################
########## 1. Variance Stabilizing Transform
#############################################


runVarianceStabilizingTransform <- function(rawcountDataframe, outfile) {

	# Load library
	library(DESeq2)

	# Fix rownames
	rownames(rawcountDataframe) <- rawcountDataframe$gene_symbol

	# Remove gene symbol column
	rawcountDataframe$gene_symbol <- NULL

	# Convet to integer matrix
	integerRawcountMatrix <- round(as.matrix(rawcountDataframe))

	# Perform VST normalization
	vstMatrix <- varianceStabilizingTransformation(integerRawcountMatrix)

	# Re-add gene symbol
	vstMatrix <- cbind(gene_symbol=rownames(vstMatrix), vstMatrix)

	# Save
	write.table(vstMatrix, file=outfile, row.names=FALSE, quote=FALSE, sep='\t')
}

#############################################
########## 2. Size Factor Normalization
#############################################

runSizeFactorNormalization <- function(rawcountDataframe, outfile) {

	# Load library
	library(DESeq2)

	# Fix rownames
	rownames(rawcountDataframe) <- rawcountDataframe$gene_symbol

	# Remove gene symbol column
	rawcountDataframe$gene_symbol <- NULL

	# Convet to integer matrix
	integerRawcountMatrix <- round(as.matrix(rawcountDataframe))

	# Calculate size factors
	sizeFactors <- estimateSizeFactorsForMatrix(integerRawcountMatrix)

	# Correct gene expression values
	normalizedMatrix <- t(t(integerRawcountMatrix)/sizeFactors)

	# Log-transform
	logNormalizedMatrix <- log10(normalizedMatrix + 1)

	# Re-add gene symbol
	logNormalizedMatrix <- cbind(gene_symbol=rownames(logNormalizedMatrix), logNormalizedMatrix)

	# Save
	write.table(logNormalizedMatrix, file=outfile, row.names=FALSE, quote=FALSE, sep='\t')
}

#######################################################
#######################################################
########## S3. Remove Batch Effects
#######################################################
#######################################################

#############################################
########## 1. Remove batch effects
#############################################

remove_batch_effects <- function(vstDataframe, annotationDataframe, outfile) {
	# Load libraries
	library(sva)
	
	# Create design matrix
	timepointDesign <- model.matrix(~timepoint, data=annotationDataframe)

	# Get batch
	batch <- annotationDataframe$batch

	# Correct data
	vstDataframeCorrected <- ComBat(dat=vstDataframe, batch=batch, mod=timepointDesign, par.prior=TRUE, prior.plots=FALSE)
	vstDataframeCorrected <- cbind(gene_symbol=rownames(vstDataframeCorrected), vstDataframeCorrected)

	# Save file
	write.table(vstDataframeCorrected, file=outfile, sep='\t', quote=FALSE, row.names=FALSE)
}

#######################################################
#######################################################
########## S4. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. Characteristic Direction
#############################################

run_characteristic_direction <- function(expressionDataframe, outfile) {

	# Read data
	source('pipeline/scripts/support/nipals.R')
	source('pipeline/scripts/support/chdir.R')

	# Get timepoints
	timepoints <- unique(sapply(colnames(expressionDataframe), function(x) strsplit(x, '.', fixed=TRUE)[[1]][1]))
	    
	# Get column names grouped by timepoint
	timepointColnames <- sapply(timepoints, function(x) colnames(expressionDataframe)[grepl(x, colnames(expressionDataframe))])
	    
	# Split and filter data
	timepointData <- sapply(timepointColnames, function(x) {
	    resultDataframe <- round(expressionDataframe[,x], digits=2);
	    geneVar <- apply(resultDataframe, 1, var);
	    resultDataframe <- resultDataframe[names(geneVar)[geneVar > 0],];
	    return(resultDataframe);
	}, simplify=FALSE)
	    
	# Get common genes
	commonGenes <- Reduce(intersect, sapply(timepointData, function(x) rownames(x)))[1:500]
	    
	# Filter data
	timepointData <- sapply(timepointData, function(x) x[commonGenes,])
	    
	# Run CD
	cdList <- sapply(timepoints[-1], function(x) chdir(timepointData[['h0']], timepointData[[x]], commonGenes)[commonGenes,])
	
	# Fix dimension names
	colnames(cdList) <- paste0('h0v', colnames(cdList))
	cdList <- cbind(gene_symbol=rownames(cdList), cdList)
	    
	# Save
	write.table(cdList, file=outfile, sep='\t', quote=FALSE, row.names=FALSE)
}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################