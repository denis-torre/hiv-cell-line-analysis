#################################################################
#################################################################
############### Datasets2Tools Database Pipeline ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys, rpy2
import pandas as pd
import rpy2.robjects as robjects
import pandas.rpy.common as com

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
import Support as S 

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
rawDataFile = 'rawdata/new_norm_geneData.txt'
annotationFile = 'rawdata/hiv_sample_annotation.txt'

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-hiv-signature-analysis.R'
r = robjects.r
r.source(rSource)

#######################################################
#######################################################
########## S1. Process Dataset
#######################################################
#######################################################

#############################################
########## 1. Get rawcount table
#############################################

@follows(mkdir('f1-expression_data.dir'))

@files(rawDataFile,
	   'f1-expression_data.dir/hiv-rawcounts.txt')

def getRawcountTable(infile, outfile):

	# Read table
	expressionDataframe = pd.read_table(infile)

	# Define columns to keep
	columnsToKeep = ['Symbol', 'B10C', 'B11C', 'B11G', 'B11N', 'B12G', 'B12N', 'B1C', 'B1G', 'B1N', 'B2C', 'B2G', 'B2N', 'B3C', 'B3G', 'B3N', 'B8N', 'B9C', 'B9G', 'B9N', 'h0.2', 'h0.3', 'h0.4', 'h12.3', 'h12.4', 'h24.3', 'h24.4', 'h48.2', 'h48.3', 'h48.4', 'h6.2', 'h6.3', 'h6.4', 'NK1', 'NK2']

	# Get subset and rename column
	rawcountDataframe = expressionDataframe[columnsToKeep].rename(columns={'Symbol': 'gene_symbol'})

	# Write file
	rawcountDataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Get primary podocyte rawcounts
#############################################

@merge([getRawcountTable,
		annotationFile],
	   'f1-expression_data.dir/primary_podocyte-rawcounts.txt')

def getPrimaryPodocyteRawcounts(infiles, outfile):

	# Split infiles
	rawcountFile, annotationFile = infiles

	# Get rawcount dataframe
	rawcountDataframe = pd.read_table(rawcountFile).set_index('gene_symbol')

	# Get annotation dataframe
	annotationDataframe = pd.read_table(annotationFile).set_index('sample_name')

	# Get primary samples
	primarySamples = annotationDataframe.index[annotationDataframe['cell_type'] == 'primary_podocyte']

	# Get dataframe subset
	rawcountDataframeSubset = rawcountDataframe.loc[:, primarySamples].dropna()

	# Write file
	rawcountDataframeSubset.to_csv(outfile, sep='\t', index=True)

#############################################
########## 3. Get cell-line podocyte rawcounts
#############################################

@merge([getRawcountTable,
		annotationFile],
	   'f1-expression_data.dir/podocyte_cell_line-rawcounts.txt')

def getPodocyteCellLineRawcounts(infiles, outfile):

	# Split infiles
	rawcountFile, annotationFile = infiles

	# Get rawcount dataframe
	rawcountDataframe = pd.read_table(rawcountFile).set_index('gene_symbol')

	# Get annotation dataframe
	annotationDataframe = pd.read_table(annotationFile).set_index('sample_name')

	# Get primary samples
	cellLineSamples = annotationDataframe.index[annotationDataframe['cell_type'] == 'cell_line']

	# Get dataframe subset
	rawcountDataframeSubset = rawcountDataframe.loc[:, cellLineSamples].dropna()
	
	# Write file
	rawcountDataframeSubset.to_csv(outfile, sep='\t', index=True)


#######################################################
#######################################################
########## S2. Normalize Data
#######################################################
#######################################################

#############################################
########## 1. Variance Stabilizing Transform
#############################################

@follows(mkdir('f2-normalized_expression_data.dir'))

@transform([getPrimaryPodocyteRawcounts,
		    getPodocyteCellLineRawcounts],
		   regex(r'.*/(.*)-rawcounts.txt'),
		   r'f2-normalized_expression_data.dir/\1-vst.txt')

def runVstNormalization(infile, outfile):

	# Read infile
	rawcountDataframe = pd.read_table(infile)

	# Run function
	r.runVarianceStabilizingTransform(com.convert_to_r_dataframe(rawcountDataframe), outfile)

#############################################
########## 2. Log10 Size Factors
#############################################

@transform([getPrimaryPodocyteRawcounts,
		    getPodocyteCellLineRawcounts],
		   regex(r'.*/(.*)-rawcounts.txt'),
		   r'f2-normalized_expression_data.dir/\1-normcounts.txt')

def runSizeFactorNormalization(infile, outfile):

	# Read infile
	rawcountDataframe = pd.read_table(infile)

	# Run function
	r.runSizeFactorNormalization(com.convert_to_r_dataframe(rawcountDataframe), outfile)

#######################################################
#######################################################
########## S3. Remove Batch Effects
#######################################################
#######################################################

#############################################
########## 1. Remove batch effects
#############################################

@follows(mkdir('f3-adjusted_expression_data.dir'))

@files('f2-normalized_expression_data.dir/podocyte_cell_line-vst.txt',
	   'f3-adjusted_expression_data.dir/podocyte_cell_line-vst_corrected.txt')

def removeBatchEffects(infile, outfile):

	# Read data
	vstDataframe = pd.read_table(infile).set_index('gene_symbol')

	# Create annotation dataframe
	annotationDataframe = pd.DataFrame([[x, x.split('.')[0], x.split('.')[1]] for x in vstDataframe.columns], columns=['sample_id','timepoint','batch'])

	# Run function
	r.remove_batch_effects(com.convert_to_r_dataframe(vstDataframe), com.convert_to_r_dataframe(annotationDataframe), outfile)

#######################################################
#######################################################
########## S4. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. Characteristic Direction
#############################################

@follows(mkdir('f4-differential_expression.dir'))

@transform(['f2-normalized_expression_data.dir/podocyte_cell_line-vst.txt',
		    'f3-adjusted_expression_data.dir/podocyte_cell_line-vst_corrected.txt'],
		    regex(r'.*/(.*).txt'),
		    r'f4-differential_expression.dir/\1-differential_expression.txt')

def runCharacteristicDirection(infile, outfile):

	# Read infile
	expressionDataframe = pd.read_table(infile).set_index('gene_symbol')

	# Split column names
	treatedColumnBool = [x.split('.')[0] != 'h0' for x in expressionDataframe.columns]
	controlColumnBool = [not x for x in treatedColumnBool]

	# Get treated dataframe
	treatedExpressionDataframe = expressionDataframe.loc[:, treatedColumnBool]
	controlExpressionDataframe = expressionDataframe.loc[:, controlColumnBool]

	# Run function
	r.run_characteristic_direction(com.convert_to_r_dataframe(treatedExpressionDataframe), com.convert_to_r_dataframe(controlExpressionDataframe), outfile)

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
