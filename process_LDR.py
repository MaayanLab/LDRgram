# I will process the LDR data with this script 

def main():

	# extract data and save to flattened matrices 
	###############################################
	# extract the nodes
	nodes = extract_nodes()
	# construct actual array
	construct_array(nodes)

	# LDR make clust 
	#########################
	make_ldr_clust()

def construct_array(nodes):
	import json_scripts
	import scipy

	print('\nconstructing array')

	# load the LDR data is json format 
	ldr = json_scripts.load_to_dict('LDR/LDR_api.json')

	# initialize matrix
	# mat = scipy.zeros([ len(nodes['as']), len(nodes['cl']), len(nodes['pt'])  ])
	# make 2d matrix for now 
	mat = scipy.zeros([ len(nodes['as']), len(nodes['cl']) ])

	# generate two released matrices 
	rl = {}
	rl['t'] = scipy.zeros([ len(nodes['as']), len(nodes['cl']) ])
	rl['f'] = scipy.zeros([ len(nodes['as']), len(nodes['cl']) ])

	# print(mat.shape)
	# print(mat[:,:,0].shape)
	# print(mat[:,:,1].shape)

	# print(nodes['as'])
	# print('\n')
	# print(nodes['cl'])
	# print('\n')
	# print(nodes['pt'])

	total = 0


	# loop through the ldf datasets 
	for inst_ldr in ldr:

		# get the inst_assay 
		inst_as = inst_ldr['datasetName'].strip()

		# print('inst_as: '+ inst_as)

		# get the cell line(s)
		inst_cls = [] 
		for inst_cl in inst_ldr['metadata']['cellLines']:
			if 'name' in inst_cl:
				inst_cls.append( inst_cl['name'].strip() )
				# print('\tcl\t'+inst_cl['name'].strip())


		# get the perturbations 
		inst_pts = []
		for inst_pt in inst_ldr['metadata']['perturbagens']:
			inst_pts.append( inst_pt['name'].strip() )
			# print('\tpt\t'+inst_pt['name'].strip())


		# add information to mat
		# get index of assay 
		index_as = nodes['as'].index(inst_as)


		# loop through cell lines
		for inst_cl in inst_cls:

			# get the index of the cell line
			index_cl = nodes['cl'].index(inst_cl)

			for inst_pt in inst_pts:

				# # get the index of the perturbation
				# index_pt = nodes['pt'].index(inst_pt)

				# check if the perturbation represents multiple perturbations 
				if 'compounds' in inst_pt and 'among' not in inst_pt:
					mult_pts = int(inst_pt.split(' ')[0])
				else:
					mult_pts = 0

				# # add value to matrix 
				# mat[ index_as, index_cl ] = mat[ index_as, index_cl ] + 1

				# print( 'released: ' + str(inst_ldr['released']) )

				# print('\n\n')

				if mult_pts == 0:
					mat[ index_as, index_cl ] = mat[ index_as, index_cl ] + 1

					# print('total: 1')

					# track number of released 
					if inst_ldr['released'] == True:
						# add to number of released 
						rl['t'][index_as, index_cl] = rl['t'][index_as, index_cl] + 1
						# print('released')
					else:
						rl['f'][index_as, index_cl] = rl['f'][index_as, index_cl] + 1
						# print('not released')

				else:
					mat[ index_as, index_cl ] = mat[ index_as, index_cl ] + mult_pts

					# track number of released 
					if inst_ldr['released'] == True:
						# add to number of released 
						rl['t'][index_as, index_cl] = rl['t'][index_as, index_cl] + mult_pts
					else:
						rl['f'][index_as, index_cl] = rl['f'][index_as, index_cl] + mult_pts

				# add to total 
				total = total + 1

				# if mat[ index_as, index_cl ] > 2:
				# 	print(mat[ index_as, index_cl ])				

				# print(index_assay)
				# print(index_cl)
				# print(index_pt)	

	# print('\n\n')
	# print(nodes['as'][0])
	# print(nodes['cl'][0])
	# print(mat[1,:])
	# print(mat[1,23])
	# print(nodes['cl'][23])

	# print('\n\n'+str(total))
	# save the matrix 
	mat = mat.tolist()
	rl['t'] = rl['t'].tolist() 
	rl['f'] = rl['f'].tolist() 

	# save the list 
	ldr_mat = {}
	ldr_mat['nodes'] = nodes
	ldr_mat['mat'] = mat
	ldr_mat['rl'] = rl

	json_scripts.save_to_json( ldr_mat, 'ldr_mat.json', 'no-indent' )

def extract_nodes():
	import json_scripts

	print('extracting nodes: as, cl, pt')

	# load the LDR data is json format 
	ldr = json_scripts.load_to_dict('LDR/LDR_api.json')

	# first generate lists of cell_lines, assays, and perturbagens 
	nodes = {}
	nodes['as'] = []
	nodes['cl'] = []
	nodes['pt'] = []
	nodes['ct'] = []

	# loop the ldr list 
	for inst_ldr in ldr:

		# add assay (datasetName)
		nodes['as'].append( inst_ldr['datasetName'].strip() )

		# get center name 
		nodes['ct'].append(inst_ldr['group']['name'])

		# # get release 
		# print( 'released: ' + str(inst_ldr['released']) )

		# add cell line(s)
		for inst_cl in inst_ldr['metadata']['cellLines']:
			if 'name' in inst_cl:
				nodes['cl'].append( inst_cl['name'].strip() )

		# add perturbation(s)
		for inst_pt in inst_ldr['metadata']['perturbagens']:
			nodes['pt'].append( inst_pt['name'].strip() )

	# get unique and sort 
	for inst_key in nodes:
		nodes[inst_key] = list(set(nodes[inst_key]))
		nodes[inst_key] = sorted(nodes[inst_key])

		print( 'there are ' + str(len(nodes[inst_key])) + ' ' + inst_key )

	return nodes

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