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

remove_batch_effects <- function(vstDataframe, annotationDataframe, outfile)
{
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
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################