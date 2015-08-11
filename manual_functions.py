# This module contains functions that will be run manually 
def main():

	# LDR make clust 
	#########################
	make_ldr_clust()

def make_ldr_clust():
	import json_scripts
	import numpy as np
	import d3_clustergram 

	# load LDR data
	ldr = json_scripts.load_to_dict('ldr_mat.json')

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

	# # last minute cleaning up of row/col names 
	# for i in range(len(nodes['col'])):
	# 	nodes['col'][i] = nodes['col'][i].replace('/ single drugs','')
	# for i in range(len(nodes['row'])):
	# 	nodes['row'][i] = nodes['row'][i].replace('cell lines','')

	# write the clustergram 
	d3_clustergram.write_json_single_value( nodes, clust_order, ldr, full_path, row_class, col_class)


# run main 
main()