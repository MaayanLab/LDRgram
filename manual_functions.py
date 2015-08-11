# This module contains functions that will be run manually 
def main():

	# # generate json of tf and target genes 
	# save_tf_targets_union_json('union')
	# save_tf_targets_union_json('int')

	# # generate tf clustergrams 
	# make_all_tf_cor_clustergrams()

	# # make clustergram of gene expression correlated with TF expression 
	# make_gene_exp_corr_with_tf_exp()

	# # make expression clustergrams 
	# make_expression_clustergrams()

	# chik specific code 
	#########################

	# make chik_log2 clustergram
	################################
	# cutoffs = [0.0,0.25,0.5, 0.75,1.0]
	# # cutoffs = [1.0]
	# for inst_cutoff in cutoffs:
	# 	make_chik_log2_clust(inst_cutoff)

	# # get chik up/dn lists 
	# get_chik_updn_lists()

	# LDR make clust 
	#########################
	make_ldr_clust()

def make_ldr_clust():
	import json_scripts
	import numpy as np
	import d3_clustergram 

	# load LDR data
	ldr = json_scripts.load_to_dict('LDR/ldr_mat.json')

	print(ldr.keys())

	ldr['mat'] = np.asarray(ldr['mat'])
	ldr['rl']['t'] = np.asarray(ldr['rl']['t'])
	ldr['rl']['f'] = np.asarray(ldr['rl']['f'])

	print( 'sum all \t' + str(np.sum(ldr['mat'])) )
	print( 'sum yes \t' + str(np.sum(ldr['rl']['t'])) )
	print( 'sum no  \t' + str(np.sum(ldr['rl']['f'])) )

	print(len(ldr['nodes']['as']))
	print(len(ldr['nodes']['cl']))
	print(ldr['mat'].shape)

	# define nodes: unfiltered
	nodes_uf = {}
	nodes_uf['row'] = ldr['nodes']['as']
	nodes_uf['col'] = ldr['nodes']['cl']

	# define parameters
	compare_cutoff = 0.05
	min_num_compare = 2

	# filter to remove nodes with no values 
	ldr['mat'], nodes = d3_clustergram.filter_sim_mat( ldr['mat'], nodes_uf, 1, 1 )
	# cherrypick using hte nodes 
	ldr['rl']['t'] = d3_clustergram.cherrypick_mat_from_nodes(nodes_uf, nodes, ldr['rl']['t'])
	ldr['rl']['f'] = d3_clustergram.cherrypick_mat_from_nodes(nodes_uf, nodes, ldr['rl']['f'])

	print( 'size all \t' + str(ldr['mat'].shape) )
	print( 'size yes \t' + str(ldr['rl']['t'].shape) )
	print( 'size no  \t' + str(ldr['rl']['f'].shape) )	
	print('\n')

	print( 'sum all \t' + str(np.sum(ldr['mat'])) )
	print( 'sum yes \t' + str(np.sum(ldr['rl']['t'])) )
	print( 'sum no  \t' + str(np.sum(ldr['rl']['f'])) )	
	print( 'total yes/no:\t' + str( np.sum(ldr['rl']['t']) + np.sum(ldr['rl']['f']) ) )

	print('\n\n\n')
	# print out nodes 
	for inst_row in nodes['row']:
		print(inst_row)

		
	print('\n\n\n')
	# print out nodes 
	for inst_row in nodes['row']:
		print(inst_row)

	print('\n\n\n')

	# cluster rows and columns 
	print('calculating clustering')
	clust_order = d3_clustergram.cluster_row_and_column( nodes, ldr['mat'], 'cosine', compare_cutoff, min_num_compare )

	print('finished calculating clustering')

	# write the d3_clustergram 
	base_path = 'static/networks/'
	full_path = base_path + 'LDR_as_cl.json'

	# add class information 
	row_class = {}
	col_class = {}

	print(len(nodes['row']))
	print(len(nodes['col']))

	# last minute cleaning up of row/col names 
	for i in range(len(nodes['col'])):
		nodes['col'][i] = nodes['col'][i].replace('/ single drugs','')


	for i in range(len(nodes['row'])):
		nodes['row'][i] = nodes['row'][i].replace('cell lines','')

	# write the clustergram 
	d3_clustergram.write_json_single_value( nodes, clust_order, ldr, full_path, row_class, col_class)




def get_chik_updn_lists():
	import json_scripts
	import numpy as np
	import d3_clustergram

	# load chik_log2 data and convert to array
	chik = json_scripts.load_to_dict('chik_log2.json')
	chik['mat'] = np.asarray(chik['mat'])

	# apply log2 cutoff of 1 
	inst_cutoff = 1

	# temporarily replace nans with zeros
	print('replace nans with zeros')
	chik['mat'] = np.nan_to_num(chik['mat'])

	# filter the matrix to only include rows/cols that have
	# values above a certain thresh n times 
	chik['mat'], chik['nodes'] = d3_clustergram.filter_sim_mat(chik['mat'],chik['nodes'],inst_cutoff,1)

	print(chik['nodes']['row'])
	print(len(chik['nodes']['row']))

	print('\n\n')
	print(chik['mat'].shape)
	print('\n\n')

	# find keratins
	chick_lists = {}
	chick_lists['keratin'] = []
	chick_lists['collagen'] = []
	chick_lists['up'] = []
	chick_lists['dn'] = []
	for i in range(len(chik['nodes']['row'])):

		inst_gene = chik['nodes']['row'][i]

		# get the protein level for each gene
		inst_value = chik['mat'][i,:]
		inst_mean = np.nanmean(inst_value)

		# collect keratins: dn 
		if 'KRT' in inst_gene:
			chick_lists['keratin'].append(inst_gene)
			# print( inst_gene + '\t' + str(inst_mean) )

		# collect collagens: up 
		if 'COL' in inst_gene:
			chick_lists['collagen'].append(inst_gene)
			# print( inst_gene + '\t' + str(inst_mean) )

		# collect up regulated genes 
		if inst_mean > 0 :
			chick_lists['up'].append(inst_gene)
			# print( 'up\t' + inst_gene + '\t' + str(inst_mean) )

		# collect dn regulated genes 
		if inst_mean < 0 :
			chick_lists['dn'].append(inst_gene)
			# print( 'dn\t' + inst_gene + '\t' + str(inst_mean) )

	# export lists to document 
	filename = 'chikungunya_lists.txt'
	fw = open(filename, 'w')

	# export all lists 
	for inst_key in chick_lists:

		# get list 
		inst_list = chick_lists[inst_key]

		# write list name 
		fw.write(inst_key +'\n')

		# write list components
		for inst_gene in inst_list:
			fw.write(inst_gene+'\n')

		fw.write('\n\n')

	fw.close()


def check_chik_prot_class(chik):
	import json_scripts

	# load gene class information from harmonogram 
	gc = json_scripts.load_to_dict('gene_classes_harmonogram.json')

	# iniitalize row_class dictionary
	row_class = {}

	# loop through the gene rows
	for inst_gene in chik['nodes']['row']:

		# check all protein classes 
		for inst_class in gc:

			# check specific protein class
			if inst_gene in gc[inst_class]:
				print('gene '+ inst_gene + ' is a ' + inst_class)

				# add gene class to dictionary 
				row_class[inst_gene] = inst_class

		# if there is no protein class then set to other
		if inst_gene not in row_class:
			row_class[inst_gene] = 'other'

	return row_class


def check_col_class(cols):

	col_class = {}
	for inst_col in cols:

		if 'lysate' in inst_col:

			col_class[inst_col] = 'lysate'
		else:
			col_class[inst_col] = 'supnt'

	return col_class

# make chik_log2 clustergram
def make_chik_log2_clust(inst_cutoff):
	import json_scripts
	import numpy as np
	import d3_clustergram

	# # only run once 
	# ################
	# # process the data 
	# import process_chikungunya_kea
	# process_chikungunya_kea.load_chikungunya_data()

	# quick attempt to fix recursion error
	import sys
	sys.setrecursionlimit(10000)

	# load chik_log2 data and convert to array
	chik = json_scripts.load_to_dict('chik_log2.json')
	chik['mat'] = np.asarray(chik['mat'])

	print(len(chik['nodes']['row']))

	# temporarily replace nans with zeros
	print('replace nans with zeros')
	chik['mat'] = np.nan_to_num(chik['mat'])

	# filter the matrix to only include rows/cols that have
	# values above a certain thresh n times 
	chik['mat'], chik['nodes'] = d3_clustergram.filter_sim_mat(chik['mat'],chik['nodes'],inst_cutoff,1)

	print(chik['mat'].shape)

	# define parameters
	compare_cutoff = 0.05
	min_num_compare = 2

	# cluster rows and columns 
	print('calculating clustering')
	clust_order = d3_clustergram.cluster_row_and_column( chik['nodes'], chik['mat'], 'cosine', compare_cutoff, min_num_compare )

	print('finished calculating clustering')

	# write the d3_clustergram 
	base_path = 'static/networks/'
	full_path = base_path + 'chik_cutoff_' + str(inst_cutoff) + '.json'

	# add class information 
	row_class = {}
	col_class = {}

	# add class information to the proteins 
	row_class = check_chik_prot_class( chik )	
	col_class = check_col_class( chik['nodes']['col'] )

	# write the clustergram 
	d3_clustergram.write_json_single_value( chik['nodes'], clust_order, chik['mat'], full_path, row_class, col_class)


# make multiple CCLE NSCLC expression clustergrams 
# at different zscore cutoffs 
def make_expression_clustergrams():
	import json_scripts
	import os

	print(os.getcwd())	

	# load gene class information from harmonogram 
	gc = json_scripts.load_to_dict('gene_classes_harmonogram.json')

	# define cutoffs
	# z_cutoffs = [ 3.5 ]
	z_cutoffs = [  0, 2, 2.5, 3, 3.5, 4 ]

	# define protein types 
	# prot_types = ['TF']
	prot_types = ['KIN','PP','ACT','DACT','MET','DMET','TF', 'GPCR', 'IC' 'all']
	# prot_types = ['GPCR','IC']
	# prot_types = gc.keys()

	# minimum number intersect
	min_num_compare = 3
	# compare cutoff - the minimum absolute value that will be compared
	compare_cutoff = 0.1

	# loop through zscore cutoffs
	for inst_z in z_cutoffs:
		# loop through protein types 
		for inst_type in prot_types:
			print( '\n\n' + inst_type + '\t' + str(inst_z) + '\n')
			# cluster and output json
			cluster_zscore(inst_type, inst_z, compare_cutoff, min_num_compare)

def cluster_zscore( inst_type, z_cutoff, compare_cutoff, min_num_compare ):
	import json_scripts
	import d3_clustergram
	import numpy as np

	# load ccle data with zscore normalization
	ccle = json_scripts.load_to_dict('CCLE/nsclc_allzc.json')
	# convert zscored data to nparray 
	ccle['data_z'] = np.asarray(ccle['data_z'], dtype = float)

	# load gs_list to filter for kinase, tf etc
	# gs_list = json_scripts.load_to_dict('enz_and_tf_lists_gmts/categories_gs_list.json')
	gs_list = json_scripts.load_to_dict('gene_classes_harmonogram.json')

	# generate node lists 
	nodes = {}
	# get all genes 
	nodes['row'] = filter_genes(ccle, gs_list, inst_type, z_cutoff)
	# get all cell lines from CCLE 
	nodes['col'] = ccle['cell_lines']

	# only make clustergram if there are genes remaining after filtering 
	# the minimum needed to produce a clustergram 
	if len(nodes['row']) > 1:
		# print(len(nodes['row']))

		print('there are ' + str(len(nodes['row'])) + ' genes being clustered' )

		# Generate data_mat: used to filter data from the original ccle for a subset of genes 
		# takes inputs: node lists for rows and columns, and primary data that will be used to make the matrix
		# the last two arguments are the names of the rows and columns in the original data
		data_mat = d3_clustergram.generate_data_mat_array( nodes, ccle, 'gene', 'cell_lines', 'data_z' )

		# cluster rows and columns 
		clust_order = d3_clustergram.cluster_row_and_column( nodes, data_mat, 'cosine', compare_cutoff, min_num_compare )

		# write the d3_clustergram 
		# base_path = '/Applications/XAMPP/xamppfiles/htdocs/cst_gram/networks/'
		base_path = 'static/networks/'
		full_path = base_path + inst_type + '_exp_std_' + str(z_cutoff) + '.json'

		# write the clustergram 
		d3_clustergram.write_json_single_value(nodes, clust_order, data_mat, full_path)

# gather all genes and cell lines 
def filter_genes(ccle, gs_list, inst_type, z_cutoff):
	import numpy as np 

	# get all gene names 
	all_genes = ccle['gene']

	# define the list of all genes as all the genes in ccle
	gs_list['all'] = all_genes 

	# filter for zscore here 
	#########################
	filtered_genes = []

	# loop through the genes and check if they are part of the gene class of interest 
	for inst_gene in all_genes:

		# check if gene of the type of interest 
		if inst_gene in gs_list[inst_type]:
			
			# get index of gene 
			inst_index =  all_genes.index(inst_gene)

			# get all zscores 
			inst_zscores = ccle['data_z'][inst_index,:]

			# get maximum zscore
			inst_max_z = np.amax(np.absolute(inst_zscores)) 

			# check if gene has zscore above cutoff 
			if inst_max_z > z_cutoff:
				filtered_genes.append(inst_gene)

	return filtered_genes


# make clustergram of gene expression correlated with TF expression 
def make_gene_exp_corr_with_tf_exp():
	import json_scripts
	import d3_clustergram
	import numpy as np
	import make_exp_clustergram

	print('make_gene_exp_corr_with_tf_exp')

	print('loading gene classes')
	# load gene class information from harmonogram
	gs_list = json_scripts.load_to_dict('gene_classes_harmonogram.json')

	print('loading tf targets - union')
	# load tf target information 
	tf_union = json_scripts.load_to_dict('enz_and_tf_lists_gmts/TF/tf_union.json')

	# load ccle data 
	##################
	# load ccle data with zscore normalization
	print('loading ccle data')
	ccle = json_scripts.load_to_dict('CCLE/nsclc_allzc.json')
	# convert zscored data to nparray 
	ccle['data_z'] = np.asarray(ccle['data_z'], dtype = float)

	# define cutoffs - this is redundant with later cutoffs 
	z_cutoffs = [2]

	# define protein types 
	prot_types = ['KIN']

	# the minimum number of comparisons of values above the threshold
	num_comp = 5
	# zscore cutoff - only compare data that is absolute value above zscore_cutoff
	# do not calculate the correlation of values with very low zscores since this is not interesting
	zscore_cutoff = 1.0
	# each row/col in the sim_mat must have a similarity value of absolute value greater than sim_cutoff
	sim_cutoff = 0.5
	# the minimum number of similarity values in a row/col that meet the similarity cutoff
	# the logic is that you do not want to have a lot of stray rows and cols with only one value since
	# these will not cluster well
	min_num_cutoff = 30
	min_meet_thresh = 7

	# generate the similarity clustergrams 
	for inst_z in z_cutoffs:

		# loop through the protein types
		for inst_type in prot_types:

			# make clustergrams and save json 

			# generate node lists 
			nodes = {}
			# get all genes 
			nodes['row'] = make_exp_clustergram.filter_genes(ccle, gs_list, inst_type, inst_z)
			# get all tfs 
			# only consider tfs that have a zscore greater than tf_zscore_cutoff in at least one cell line
			tf_zscore_cutoff = 2
			nodes['col'] = make_exp_clustergram.filter_genes(ccle, gs_list, 'TF', tf_zscore_cutoff) 

			# calculate the similarity matrix 
			sim_mat, nodes = d3_clustergram.generate_sim_mat_array( nodes, ccle, 'gene', 'gene', 'data_z', num_comp, zscore_cutoff, sim_cutoff, min_meet_thresh )

			# improve cosine distance clustering 

			print('clustering')

			# cluster rows and columns of sim_mat 
			# sim cutoff = 0
			# num_comp = 3
			clust_order = d3_clustergram.cluster_row_and_column( nodes, sim_mat, 'cosine', 0, 3)

			# write the d3_clustergram 
			base_path = 'static/networks/'
			full_path = base_path + 'TF_sim_' + inst_type + '_exp_std_' + str(inst_z) + '.json'

			print('writing network to json')

			# write the clustergram 
			# the last two arguments are targets and evidence 
			# targets are used for a clustergram of tfs vs genes and the rows are labeled as targets 
			# evidence is used to label individual tiles if the gene is known to be downstream of a tf 
			d3_clustergram.write_json_single_value(nodes, clust_order, sim_mat, full_path, [], tf_union)


# make a clustergram for each transcription factor 
def make_all_tf_cor_clustergrams():
	import json_scripts

	# load gene symbol lists 
	gs_list = json_scripts.load_to_dict('gene_classes_harmonogram.json')

	# load ccle data with zscore normalization
	ccle = json_scripts.load_to_dict('CCLE/nsclc_allzc.json')

	all_tfs = []

	# loop through all tfs
	for inst_tf in gs_list['TF']:

		# check that TF is a measured gene in the ccle data 
		if inst_tf in ccle['gene']:

			all_tfs.append(inst_tf)

			# generate the tf cor clustergram 
			tf_clust(inst_tf)

			# print('\n')

	print(len(all_tfs))

# save the gmt to a json 
def save_tf_targets_union_json(inst_type):
	import calc_enrichment_gl_gmt
	import json_scripts

	# load gmt into json 
	tf_union = calc_enrichment_gl_gmt.load_gmt('enz_and_tf_lists_gmts/TF/tf_' + inst_type + '.gmt')

	print( 'num tfs ' + str(len(tf_union.keys())) )

	# save json 
	json_scripts.save_to_json(tf_union,'enz_and_tf_lists_gmts/TF/tf_'+ inst_type +'.json','no_indent')

# make tf-sub clustergram 
def tf_clust(tf_name):
	import json_scripts
	import d3_clustergram
	import numpy as np
	from operator import itemgetter 

	# load gs_list to filter genes by type 
	gs_list = json_scripts.load_to_dict('gene_classes_harmonogram.json')

	# load ccle data with zscore normalization
	ccle = json_scripts.load_to_dict('CCLE/nsclc_allzc.json')

	# convert zscored data to nparray 
	ccle['data_z'] = np.asarray(ccle['data_z'], dtype = float)

	# load tf downstream gene data 
	tf_union = json_scripts.load_to_dict('enz_and_tf_lists_gmts/TF/tf_union.json')

	# find genes that are known to be targeted by transcription factor 
	target_genes = []
	if tf_name in tf_union:
		target_genes = tf_union[tf_name]

	# add tf to list 
	target_genes.extend(tf_name)

	# find top 100 most similar genes 
	##################################

	# define the minimum number of intersecting measurements 
	min_num_int = 3

	# find index of tf_name in ccle 
	inst_tf_index = ccle['gene'].index(tf_name)
	# get vector 
	inst_tf_vector = ccle['data_z'][inst_tf_index,:]

	cl_sort = np.argsort(inst_tf_vector)

	# calculate the distance of each gene to the inst_tf 
	all_dist = []
	for i in range(len(ccle['gene'])):

		# get gene_vect 
		gene_vect = ccle['data_z'][i,:]

		# calculate distance 
		inst_dist = d3_clustergram.calc_dist_vectors(inst_tf_vector, gene_vect, 'cosine', min_num_int) 

		# save dictionary 
		inst_dict = {}
		inst_dict['dist'] = inst_dist
		inst_dict['name'] = ccle['gene'][i]

		# add dist to distance vector 
		all_dist.append(inst_dict)

	print('there are ' + str(len(target_genes)) + ' genes targeted by ' + str(tf_name))

	# sort 
	all_dist = sorted(all_dist, key=itemgetter('dist'))

	# keep similar and antisimilar 
	# the first gene is inst_tf 
	tf_sim_dict = all_dist[:101] + all_dist[-100:]

	# get the list of sim and anti-sim genes 
	tf_sim_list = [d['name'] for d in tf_sim_dict]

	# generate clustergram 
	#########################

	# generate node lists 
	nodes = {}
	# get all genes 
	nodes['row'] = tf_sim_list
	# get all cell lines from CCLE 
	nodes['col'] = ccle['cell_lines']

	# minimum number intersect
	min_num_int = 3

	# only make clustergram if there are genes remaining after filtering 
	# the minimum needed to produce a clustergram 
	if len(nodes['row']) > 1:

		print('there are ' + str(len(nodes['row'])) + ' genes being clustered' )

		# Generate data_mat: used to filter data from the original ccle for a subset of genes 
		# takes inputs: node lists for rows and columns, and primary data that will be used to make the matrix
		# the last two arguments are the names of the rows and columns in the original data
		data_mat = d3_clustergram.generate_data_mat_array( nodes, ccle, 'gene', 'cell_lines', 'data_z' )

		# cluster rows and columns 
		clust_order = d3_clustergram.cluster_row_and_column( nodes, data_mat, 'cosine', min_num_int )

		# reset cluster ordering 
		clust_order['clust']['row'] = list( np.arange( len(clust_order['clust']['row']) ) )
		# reverse array 
		clust_order['clust']['row'] = clust_order['clust']['row'][::-1]
		# sort cell lines based on tf expression level 
		clust_order['clust']['col'] = list(cl_sort)

		# write the d3_clustergram 
		# base_path = '/Applications/XAMPP/xamppfiles/htdocs/cst_gram/networks/'
		base_path = 'static/networks/'
		full_path = base_path + 'tf_' + tf_name + '.json'

		# write the clustergram 
		d3_clustergram.write_json_single_value(nodes, clust_order, data_mat, full_path, target_genes)


	# generate matrix of transcription factor expression and target gene expression 

	# cluster

	# return clustergram 

	return {'tmp':'something'}

# run main 
main()