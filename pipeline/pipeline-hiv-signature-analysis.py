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
import numpy as np
import rpy2.robjects as robjects
import pandas.rpy.common as com

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support as S
import PipelineHivSignatureAnalysis as P

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

@transform(removeBatchEffects,
		    regex(r'.*/(.*).txt'),
		    r'f4-differential_expression.dir/\1-differential_expression.txt')

def runCharacteristicDirection(infile, outfile):

	# Read infile
	expressionDataframe = pd.read_table(infile).set_index('gene_symbol')

	# Run function
	r.run_characteristic_direction(com.convert_to_r_dataframe(expressionDataframe), outfile)

#######################################################
#######################################################
########## S5. Enrichment Analysis
#######################################################
#######################################################

#############################################
########## 1. Geneset Submission
#############################################

@follows(mkdir('f5-geneset_enrichment.dir'))

@transform(runCharacteristicDirection,
	   	   regex(r'.*/(.*).txt'),
	   	   r'f5-geneset_enrichment.dir/\1_geneset_ids.txt')

def submitGenesets(infile, outfile):

	# Read infile
	cdDataframe = pd.read_table(infile).set_index('gene_symbol')

	# Set number of genes
	nGenes = 500

	# Get genesets
	genesets = {}
	genesets['upregulated'] = {x: cdDataframe[x].sort_values(ascending=False).index[:nGenes].tolist() for x in cdDataframe.columns}
	genesets['downregulated'] = {x: cdDataframe[x].sort_values(ascending=True).index[:nGenes].tolist() for x in cdDataframe.columns}

	# Define empty dict
	listIds = {x:{} for x in genesets.keys()}

	# Loop through directions
	for direction in genesets.keys():
	    
	    # Loop through timepoints
	    for timepoint in genesets[direction].keys():
	        
	        # Get geneset
	        geneset = genesets[direction][timepoint]
	        
	        # Submit genelist and get ID
	        listIds[direction][timepoint] = P.addGeneLists(geneset)

	# Create list dataframe
	genelistDataframe = pd.DataFrame()

	# Convert to dataframe
	for direction in listIds.keys():
	    
	    # Convert
	    genelistDataframeSubset = pd.DataFrame(listIds[direction]).T
	    
	    # Add label
	    genelistDataframeSubset['direction'] = direction
	    
	    # Fix URL
	    genelistDataframeSubset['URL'] = ['http://amp.pharm.mssm.edu/Enrichr/enrich?dataset='+x for x in genelistDataframeSubset['shortId']]
	    
	    # Append
	    genelistDataframe = pd.concat([genelistDataframe, genelistDataframeSubset])
		    
		# Write file
	genelistDataframe.to_csv(outfile, sep='\t', index=True, index_label='timepoint')	

#############################################
########## 2. Geneset Enrichment
#############################################

@transform(submitGenesets,
		   suffix('ids.txt'),
		   'enrichment.txt')

def getGenesetEnrichment(infile, outfile):

	# Read infile
	genelistDataframe = pd.read_table(infile)

	# Set geneset libraries
	genesetLibraries = ['ChEA_2016', 'KEGG_2016', 'GO_Biological_Process_2015', 'GO_Cellular_Component_2015', 'GO_Molecular_Function_2015', 'VirusMINT']

	# Define empty dataframe
	enrichmentDataframe = pd.DataFrame()

	# Loop through directions
	for direction in set(genelistDataframe['direction']):
	    
	    # Loop through timepoints
	    for timepoint in set(genelistDataframe['timepoint']):
	        
	        # Get geneset
	        listId = genelistDataframe.loc[(genelistDataframe['direction']==direction) & (genelistDataframe['timepoint']==timepoint), 'userListId'].values[0]

	        # Loop through
	        for genesetLibrary in genesetLibraries:

	        	# Get data
				enrichmentResultsDataframe = P.getEnrichmentResults(listId, gene_set_library=genesetLibrary)
				enrichmentResultsDataframe['direction'] = direction
				enrichmentResultsDataframe['timepoint'] = timepoint
				enrichmentResultsDataframe['gene_set_library'] = genesetLibrary
				enrichmentDataframe = pd.concat([enrichmentDataframe, enrichmentResultsDataframe])

	# Save outfile
	enrichmentDataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S6. L1000CDS2 Analysis
#######################################################
#######################################################

#############################################
########## 1. Geneset Submission
#############################################

@follows(mkdir('f6-l1000cds2_analysis.dir'))

@transform(runCharacteristicDirection,
	   	   regex(r'.*/(.*).txt'),
	   	   r'f6-l1000cds2_analysis.dir/\1_l1000cds2_results.txt')

def runL1000cds2Analysis(infile, outfile):

	# Read infile
	cdDataframe = pd.read_table(infile).set_index('gene_symbol')

	# Get results
	analysisResultDict = {timepoint:{aggravate: P.runL1000CDS2(cdDataframe, timepoint, aggravate) for aggravate in [True, False]} for timepoint in cdDataframe.columns}

	# Get link dictionary
	resultLinkDict = {timepoint:{aggravate:'http://amp.pharm.mssm.edu/L1000CDS2/#/result/'+analysisResultDict[timepoint][aggravate]['shareId'] for aggravate in analysisResultDict[timepoint].keys()} for timepoint in analysisResultDict.keys()}

	# Convert to dataframe
	resultLinkDataframe = pd.DataFrame(resultLinkDict).T.rename(columns={True:'mimic', False:'reverse'})

	# Get perturbation dataframe
	perturbationDataframe = pd.DataFrame()

	# Loop through timepoints
	for timepoint in analysisResultDict.keys():
	    
	    # Loop through aggravate
	    for aggravate in analysisResultDict[timepoint].keys():

	        # Get result
	        timepointPerturbationDataframe = pd.DataFrame(analysisResultDict[timepoint][aggravate]['topMeta']).drop('overlap', axis=1)

	        # Add labels
	        timepointPerturbationDataframe['timepoint'] = timepoint
	        timepointPerturbationDataframe['aggravate'] = aggravate

	        # Append
	        perturbationDataframe = pd.concat([perturbationDataframe, timepointPerturbationDataframe])

    # Add label
	perturbationDataframe['perturbation'] = [str(pert) if pert != '-666' else str(pubchem_id) for pert, pubchem_id in perturbationDataframe[['pert_desc','pubchem_id']].as_matrix()]

	# Save data
	linkOutfile = outfile.replace('results', 'links')
	resultLinkDataframe.to_csv(linkOutfile, sep='\t', index=True, index_label='timepoint')
	perturbationDataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S7. Clustergrammer Visualization
#######################################################
#######################################################

#############################################
########## 1. Prepare matrix
#############################################

@follows(mkdir('f7-clustergrammer.dir'))

@transform([removeBatchEffects,
			'f2-normalized_expression_data.dir/podocyte_cell_line-vst.txt'],
		   regex(r'.*/(.*).txt'),
		   r'f7-clustergrammer.dir/\1-clustergrammer_input.txt')

def createClustergrammerInputMatrix(infile, outfile):

	# Read infile
	expressionDataframe = pd.read_table(infile).set_index('gene_symbol')

	# Rank by variance
	rankedGenes = expressionDataframe.apply(np.var, 1).sort_values(ascending=False).index.tolist()

	# Get number of genes
	nGenes = 1000

	# Define empty list
	dataList = []

	# Get annotation
	dataList.append(['']+['Sample: '+x for x in expressionDataframe.columns])
	dataList.append(['']+['Timepoint: '+x.split('.')[0].replace('h','')+'h' for x in expressionDataframe.columns])
	dataList.append(['']+['Batch: '+x.split('.')[1] for x in expressionDataframe.columns])

	# Add genes
	for geneSymbol in rankedGenes[:nGenes]:
	    dataList.append(['Gene Symbol: '+geneSymbol]+[str(round(x, ndigits=2)) for x in expressionDataframe.loc[geneSymbol,:]])

	# Create string
	outfileString = '\n'.join(['\t'.join(x) for x in dataList])

	# Write
	with open(outfile, 'w') as openfile:
	    openfile.write(outfileString)

#############################################
########## 2. Upload matrix
#############################################

@transform([createClustergrammerInputMatrix,
			'f2-normalized_expression_data.dir/podocyte_cell_line-vst.txt'],
		   suffix('input.txt'),
		   'links.txt')

def uploadClustergrammerMatrix(infile, outfile):

	# Get Upload URL
	uploadUrl = 'http://amp.pharm.mssm.edu/clustergrammer/matrix_upload/'

	# Loop through
	r = requests.post(uploadUrl, files={'file': open(infile, 'rb')})

	# Get link
	link = r.text

	# Write
	with open(outfile, 'w') as openfile:
		openfile.write(link)

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
