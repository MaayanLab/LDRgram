# I will process the LDR data with this script 

def main():

	# load Avi's dictionary of cell lines and assays 
	assay_cl_dict()

	# extract data and save to flattened matrices 
	construct_array()

	# LDR make clust 
	#########################
	make_ldr_clust()

def assay_cl_dict():
	import json_scripts

	f = open('LDR/assays_and_cl_lists_for_Avi-AM.txt', 'r')
	lines = f.readlines()
	f.close()

	# make names dictionary 
	names = {}
	names['as'] = {}
	names['cl'] = {}

	# will go through assays and cell lines 
	inst_data = ''

	# loop through the lines 
	for inst_line in lines:

		# strip the line 
		inst_line = inst_line.strip()

		if 'assays:' in inst_line:
			inst_data = 'as'
			# print(inst_data)

		if 'cell lines:' in inst_line:
			inst_data = 'cl'
			# print(inst_data)

		# load assays 
		##############
		if inst_data == 'as':
			# check if there is a short name 
			if '\t' in inst_line:
				inst_sn = inst_line.split('\t')[0]
				inst_ln =  inst_line.split('\t')[1]

				names[inst_data][inst_ln] = inst_sn

				# # add data to dictionary 
				# if inst_ln not in names[inst_data]:
				# # add short name to dictionary 
				# names[inst_data][inst_key].append(inst_ln)

			# if there is no short name add long name as key and value 
			elif len(inst_line) > 0:
				names[inst_data][inst_line] = inst_line
				# names[inst_data][inst_line].append(inst_line)

		# load cell lines 
		###################
		if inst_data == 'cl':
			# check if there is a short name 
			if '\t' in inst_line:
				inst_sn = inst_line.split('\t')[0]
				inst_ln  = inst_line.split('\t')[1]

				names[inst_data][inst_ln] = inst_sn

				# # add data to dictionary
				# if inst_key not in names[inst_data]:
				# 	names[inst_data][inst_key] = []
				# # add short name to dictionary 
				# names[inst_data][inst_key].append(inst_ln)

			# if tehre is no short name add long name as key and value 
			elif len(inst_line) > 0:
				names[inst_data][inst_line] = inst_line
				# names[inst_data][inst_line].append(inst_line)

	# print(len(names['as'].keys()))
	# print('\n')
	# print(len(names['cl'].keys()))
	# print('\n')
	# print( len(list(set(names['cl'].values()))) )
	# print('\n')
	# print( len(list(set(names['as'].values()))) )

	json_scripts.save_to_json(names,'as_cl_dict.json','indent')

def construct_array():
	import json_scripts
	import scipy

	print('\nconstructing array\n')

	# load the LDR data is json format 
	ldr = json_scripts.load_to_dict('LDR/LDR_api.json')

	# load cl and as dictionary 
	as_cl_dict = json_scripts.load_to_dict('as_cl_dict.json')

	# get nodes from 'short name' dictionary values 
	nodes = {}
	nodes['as'] = sorted(list(set(as_cl_dict['as'].values())))
	nodes['cl'] = list(set(as_cl_dict['cl'].values()))
	# add cell-free to list of cell lines 
	nodes['cl'].append('cell-free')
	nodes['cl'] = sorted(nodes['cl'])

	# print(nodes['as'])
	# print(nodes['cl'])
	# print('\n\n\n')



	# # run once - add back removed as and cl to Avi dictionary 
	# # find assays and cell lines that were removed from original list 
	# #####################################################################
	# all_nodes = extract_nodes()
	# for inst_data in as_cl_dict:
	# 	# get all nodes
	# 	tmp_dict = set( as_cl_dict[inst_data].keys() )
	# 	tmp_all = set( all_nodes[inst_data] )
	# 	not_found = list( tmp_all - tmp_dict )
	# 	print('\n')
	# 	print(inst_data)
	# 	for tmp in not_found:
	# 		print(tmp)
	# 	print('\n')



	# make 2d matrix for now 
	mat = scipy.zeros([ len(nodes['as']), len(nodes['cl']) ])

	# generate two released matrices 
	rl = {}
	rl['t'] = scipy.zeros([ len(nodes['as']), len(nodes['cl']) ])
	rl['f'] = scipy.zeros([ len(nodes['as']), len(nodes['cl']) ])

	total = 0

	# loop through the ldf datasets 
	for inst_ldr in ldr:

		# get the inst_assay: put name through dictionary 
		# print( inst_ldr['datasetName'].strip() )
		inst_as = as_cl_dict['as'][ inst_ldr['datasetName'].strip() ]
		print('inst_as: '+ inst_as)

		# get the cell line(s)
		inst_cls = [] 


		for inst_cl in inst_ldr['metadata']['cellLines']:
			if 'name' in inst_cl:
				#!! remove cell line 'TBD among cell ...'
				if 'TBD among' not in inst_cl['name'].strip():
					inst_cls.append( as_cl_dict['cl'][ inst_cl['name'].strip() ] )


		# get the perturbations 
		inst_pts = []
		for inst_pt in inst_ldr['metadata']['perturbagens']:
			inst_pts.append( inst_pt['name'].strip() )
			# print('\tpt\t'+inst_pt['name'].strip())


		# # check if there eare cases with perturbations and no cell lines 
		# if len(inst_ldr['metadata']['cellLines']) == 0 and len(inst_pts) > 0:
		# 	print('\n\nno cell lines found ')
		# 	print(inst_ldr['datasetName'])
		# 	print(inst_pts)
		# 	print('\n\n')

		# if the assay is kinomescan then set cell line to 'cell-free
		if inst_as == 'KINOMEscan':
			print('kinomescan')
			inst_cls.append( 'cell-free' )
			print(inst_cls)
			print('\n\n\n')


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