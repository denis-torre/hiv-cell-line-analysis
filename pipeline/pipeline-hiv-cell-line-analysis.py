#################################################################
#################################################################
############### HIV Signature Analysis Pipeline #################
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
import numpy as np
import rpy2.robjects as robjects
import pandas.rpy.common as com

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support as S
import PipelineHivCellLineAnalysis as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
rawDataFile = 'rawdata.dir/new_norm_geneData.txt'

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-hiv-cell-line-analysis.R'
r = robjects.r
r.source(rSource)

#######################################################
#######################################################
########## S1. Raw Expression
#######################################################
#######################################################

#############################################
########## 1. Get rawcount table
#############################################

@follows(mkdir('f1-expression.dir'))

@files(rawDataFile,
	   'f1-expression.dir/hiv_cell_line-rawcounts.txt')

def getRawcountTable(infile, outfile):

	# Read table
	expressionDataframe = pd.read_table(infile)

	# Define columns to keep
	columnsToKeep = ['Symbol', 'h0.2', 'h0.3', 'h0.4', 'h12.3', 'h12.4', 'h24.3', 'h24.4', 'h48.2', 'h48.3', 'h48.4', 'h6.2', 'h6.3', 'h6.4']

	# Get subset and rename column
	rawcountDataframe = expressionDataframe[columnsToKeep].set_index('Symbol').dropna()

	# Convert to integer
	rawcountDataframe = rawcountDataframe.astype('int')

	# Write file
	rawcountDataframe.to_csv(outfile, sep='\t', index_label='gene_symbol')

#######################################################
#######################################################
########## S2. Annotation
#######################################################
#######################################################

#############################################
########## 1. Get sample annotations
#############################################

@follows(mkdir('f2-annotation.dir'))

@transform(getRawcountTable,
		   regex(r'.*/(.*)-rawcounts.txt'),
		   r'f2-annotation.dir/\1-annotation.txt')

def getAnnotationTable(infile, outfile):

	# Read infile
	rawcountDataframe = pd.read_table(infile, index_col='gene_symbol')

	# Get sample names
	sampleNames = rawcountDataframe.columns.tolist()

	# Define attribute dict
	annotationDict = {x:{} for x in sampleNames}

	# Loop through sample names
	for sampleName in sampleNames:
	    
	    # Get attributes
	    annotationDict[sampleName]['timepoint'] = sampleName.split('.')[0].replace('h', '')
	    annotationDict[sampleName]['batch'] = sampleName.split('.')[1].replace('h', '')
	    annotationDict[sampleName]['treatment'] = 'control' if annotationDict[sampleName]['timepoint'] == '0' else 'hiv_infected'

	# Convert to dataframe
	annotationDataframe = pd.DataFrame(annotationDict).T

	# Save file
	annotationDataframe.to_csv(outfile, sep='\t', index_label='sample_name')

#######################################################
#######################################################
########## S3. Normalize Expression
#######################################################
#######################################################

#############################################
########## 1. Variance Stabilizing Transform
#############################################

@follows(mkdir('f3-normalized_expression.dir'))

@transform(getRawcountTable,
		   regex(r'.*/(.*)-rawcounts.txt'),
		   r'f3-normalized_expression.dir/\1-vst.txt')

def runVst(infile, outfile):

	# Read expression dataframe
	rawcountDataframe = pd.read_table(infile, index_col='gene_symbol').fillna(0)

	# Run function
	vstMatrix = r.runVST(com.convert_to_r_dataframe(rawcountDataframe))

	# Convert to dataframe
	vstDataframe = com.convert_robj(vstMatrix)

	# Write file
	vstDataframe.to_csv(outfile, sep='\t', index_label='gene_symbol')

#######################################################
#######################################################
########## S4. Batch Effect Removal
#######################################################
#######################################################

#############################################
########## 1. Run ComBat
#############################################

@follows(mkdir('f4-combat.dir'))

@transform(runVst,
		   regex(r'.*/(.*)-vst.txt'),
		   add_inputs(getAnnotationTable),
		   r'f4-combat.dir/\1-combat.txt')

def runComBat(infiles, outfile):

	# Split infiles
	vstFile, annotationFile = infiles

	# Read expression dataframe
	vstDataframe = pd.read_table(vstFile, index_col='gene_symbol')

	# Read annotation dataframe
	annotationDataframe = pd.read_table(annotationFile, index_col='sample_name')

	# Run function
	combatMatrix = r.runComBat(com.convert_to_r_dataframe(vstDataframe), com.convert_to_r_dataframe(annotationDataframe), covariateFormula='~timepoint')

	# Convert to dataframe
	combatDataframe = com.convert_robj(combatMatrix)

	# Write file
	combatDataframe.to_csv(outfile, sep='\t', index_label='gene_symbol')

#######################################################
#######################################################
########## S5. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. Characteristic Direction
#############################################

@follows(mkdir('f5-characteristic_direction.dir'))

@transform(runComBat,
		   regex(r'.*/(.*).txt'),
		   add_inputs(getAnnotationTable),
		   r'f5-characteristic_direction.dir/\1_cd.txt')

def runCharacteristicDirection(infiles, outfile):

	# Split infiles
	expressionFile, annotationFile = infiles

	# Read expression file
	expressionDataframe = pd.read_table(expressionFile, index_col='gene_symbol')

	# Read annotation file
	annotationDataframe = pd.read_table(annotationFile, index_col='sample_name')

	# Get timepoint dict
	timepointDict = {str(x)+'h':annotationDataframe.index[annotationDataframe['timepoint'] == x].tolist() for x in set(annotationDataframe['timepoint'])}

	# Group 12h and 24h
	timepointDict['12-24h'] = timepointDict['12h'] + timepointDict['24h']
	del timepointDict['12h']
	del timepointDict['24h']

	# Get controls
	controlColumns = timepointDict.pop('0h')

	# Initialize empty dataframe
	resultDataframe = pd.DataFrame()

	# Loop through timepoints
	for timepoint in timepointDict.keys():

		# Get experiment samples
		experimentColumns = timepointDict[timepoint]

		# Run characteristic direction
		cdResults = r.runCharacteristicDirection(com.convert_to_r_dataframe(expressionDataframe), experimentColumns, controlColumns)

		# Convert to dataframe
		cdDataframe = com.convert_robj(cdResults).reset_index()

		# Add timepoint column
		cdDataframe['timepoint'] = timepoint

		# Append
		resultDataframe = pd.concat([resultDataframe, cdDataframe])

	# Pivot
	resultDataframeCast = resultDataframe.pivot(index='index', columns='timepoint', values='CD')

	# Save
	resultDataframeCast.to_csv(outfile, sep='\t', index_label='gene_symbol')

#######################################################
#######################################################
########## S6. Enrichment Analysis
#######################################################
#######################################################

#############################################
########## 1. Submit Genesets 
#############################################

@follows(mkdir('f6-enrichr.dir'))

@transform(runCharacteristicDirection,
		   regex(r'.*/(.*)cd.txt'),
		   r'f6-enrichr.dir/\1enrichr_links.txt')

def submitEnrichrGenesets(infile, outfile):

	# Read infile
	cdDataframe = pd.read_table(infile, index_col='gene_symbol').fillna(0)

	# Initialize link dataframe
	resultDataframe = pd.DataFrame()

	# Loop through timepoints
	for timepoint in cdDataframe.columns:

	    # Get Enrichr links
	    enrichrLinkDataframe = S.uploadToEnrichr(cdDataframe, timepoint)

	    # Add timepoint label
	    enrichrLinkDataframe['timepoint'] = timepoint

	    # Concatenate
	    resultDataframe = pd.concat([resultDataframe, enrichrLinkDataframe])

	# Save data
	resultDataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Get Enrichment results
#############################################

@transform(submitEnrichrGenesets,
		   suffix('links.txt'),
		   'results.txt')

def getEnrichrResults(infile, outfile):

	# Read infile
	enrichrLinkDataframe = pd.read_table(infile, index_col=['geneset','timepoint'])

	# Initialize result dataframe
	resultDataframe = pd.DataFrame()

	# Set libraries
	libraries = ['ChEA_2016', 'KEGG_2016', 'GO_Biological_Process_2015', 'GO_Cellular_Component_2015', 'GO_Molecular_Function_2015', 'VirusMINT']

	# Loop through timepoints, genesets and libraries
	for geneset in enrichrLinkDataframe.index.levels[0]:
	    for timepoint in enrichrLinkDataframe.index.levels[1]:
	        for library in libraries:

	            # Get enrichment results
	            enrichmentResultDataframe = S.getEnrichmentResults(enrichrLinkDataframe.loc[(geneset, timepoint), 'userListId'], library)

	            # Add labels
	            enrichmentResultDataframe['timepoint'] = timepoint
	            enrichmentResultDataframe['geneset'] = geneset
	            enrichmentResultDataframe['library'] = library

	            # Concatenate
	            resultDataframe = pd.concat([resultDataframe, enrichmentResultDataframe])

    # Write file
	resultDataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S7. L1000CDS2 Analysis
#######################################################
#######################################################

#############################################
########## 1. Submit analysis
#############################################

@follows(mkdir('f7-l1000cds2.dir'), getEnrichrResults)

@transform(runCharacteristicDirection,
		   regex(r'.*/(.*)cd.txt'),
		   r'f7-l1000cds2.dir/\1l1000cds2_links.txt')

def runL1000CDS2(infile, outfile):

	# Read infile
	cdDataframe = pd.read_table(infile, index_col='gene_symbol').fillna(0)

	# Initialize dataframes
	linkDataframe = pd.DataFrame()
	signatureDataframe = pd.DataFrame()

	# Loop through timepoints
	for timepoint in cdDataframe.columns:

	    # Run L1000CDS2
	    resultDict = S.getL1000CDS2Results(cdDataframe, timepoint)

	    # Add timepoint labels
	    resultDict['links']['timepoint'] = timepoint
	    resultDict['signatures']['timepoint'] = timepoint

	    # Append dataframes
	    linkDataframe = pd.concat([linkDataframe, resultDict['links']])
	    signatureDataframe = pd.concat([signatureDataframe, resultDict['signatures']])

	# Write files
	linkDataframe.to_csv(outfile, sep='\t', index=False)
	signatureDataframe.to_csv(outfile.replace('links', 'signatures'), sep='\t', index=False)

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
