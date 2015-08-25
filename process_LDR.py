# I will process the LDR data with this script 

def main():

	# load Avi's dictionary of cell lines and assays 
	assay_cl_dict()

	# extract data and save to flattened matrices 
	perts = construct_array()

	# LDR make clust 
	make_ldr_clust(perts)

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

	# generate perturbation dictionary that will save perturbation 
	# information for assays and cell lines 
	perts = {}

	total = 0

	# loop through the ldf datasets 
	for inst_ldr in ldr:

		# get the inst_assay: put name through dictionary 
		# print( inst_ldr['datasetName'].strip() )
		inst_as = as_cl_dict['as'][ inst_ldr['datasetName'].strip() ]
		# print('inst_as: '+ inst_as)

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


		# if the assay is kinomescan then set cell line to 'cell-free
		if inst_as == 'KINOMEscan':
			# print('kinomescan')
			inst_cls.append( 'cell-free' )
			# print(inst_cls)
			# print('\n\n\n')


		# add information to mat
		# get index of assay 
		index_as = nodes['as'].index(inst_as)


		# loop through cell lines
		for inst_cl in inst_cls:

			# get the index of the cell line
			index_cl = nodes['cl'].index(inst_cl)

			for inst_pt in inst_pts:

				# check if the perturbation represents multiple perturbations 
				if 'compounds' in inst_pt and 'among' not in inst_pt:
					mult_pts = int(inst_pt.split(' ')[0])
				else:
					mult_pts = 0

				# track the number of perturbations and the released status 
				##############################################################
				if mult_pts == 0:
					mat[ index_as, index_cl ] = mat[ index_as, index_cl ] + 1

					# track number of released 
					if inst_ldr['released'] == True:
						rl['t'][index_as, index_cl] = rl['t'][index_as, index_cl] + 1
					else:
						rl['f'][index_as, index_cl] = rl['f'][index_as, index_cl] + 1

				else:
					mat[ index_as, index_cl ] = mat[ index_as, index_cl ] + mult_pts

					# track number of released 
					if inst_ldr['released'] == True:
						rl['t'][index_as, index_cl] = rl['t'][index_as, index_cl] + mult_pts
					else:
						rl['f'][index_as, index_cl] = rl['f'][index_as, index_cl] + mult_pts

				# keep track of perturbation information in the dictionary 
				##############################################################
				# genrate as cl tuple 
				inst_tuple = (inst_as, inst_cl) 
				# initailize list if necessary 
				if inst_tuple not in perts:
					perts[inst_tuple] = []
				# generate pert_dict
				pert_dict = {}
				pert_dict['name'] = inst_pt 
				pert_dict['release'] = inst_ldr['released']
				pert_dict['_id'] = inst_ldr['_id']
				# add dictionary to list 
				perts[inst_tuple].append(pert_dict)

				# add to total 
				total = total + 1

	# check perts dictionary 
	print('perts dictionary - the number of found as/cl combinations')
	print(len(perts.keys()))
	# print(perts)

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
	# ldr_mat['perts'] = perts

	json_scripts.save_to_json( ldr_mat, 'ldr_mat.json', 'no-indent' )

	# return perts since tuple dictionaries do not save as json easily
	return perts 

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

def make_ldr_clust(perts):

	import json_scripts
	import numpy as np
	import d3_clustergram 
	from d3_clustergram_class import Network 

	# load LDR data - stored as:
	# released status (rl)
	# nodes, and mat 
	ldr = json_scripts.load_to_dict('ldr_mat.json')
	print('\nclustering')
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
	print('\n')

	print( 'size all \t' + str(ldr['mat'].shape) )
	print( 'size yes \t' + str(ldr['rl']['t'].shape) )
	print( 'size no  \t' + str(ldr['rl']['f'].shape) )	
	print('\n')

	print( 'sum all \t' + str(np.sum(ldr['mat'])) )
	print( 'sum yes \t' + str(np.sum(ldr['rl']['t'])) )
	print( 'sum no  \t' + str(np.sum(ldr['rl']['f'])) )	
	print( 'total yes/no:\t' + str( np.sum(ldr['rl']['t']) + np.sum(ldr['rl']['f']) ) )

	# define nodes: unfiltered
	nodes_uf = {}
	nodes_uf['row'] = ldr['nodes']['as']
	nodes_uf['col'] = ldr['nodes']['cl']

	# initialize a new network class 
	##################################
	net = Network()

	net.dat['nodes']['row'] = nodes_uf['row']
	net.dat['nodes']['col'] = nodes_uf['col']
	net.dat['mat'] = ldr['mat']
	net.dat['mat_up'] = ldr['rl']['t']
	net.dat['mat_dn'] = -ldr['rl']['f']

	# filter the matrix using cutoff and min_num_meet
	###################################################
	# filtering matrix 
	cutoff_meet = 1
	min_num_meet = 1
	net.filter_network_thresh( cutoff_meet, min_num_meet )

	# cluster 
	#############
	cutoff_comp = 3
	min_num_comp = 4
	net.cluster_row_and_col('cos', cutoff_comp, min_num_comp, dendro=False)

	# export data visualization to file 
	######################################
	net.write_json_to_file('viz', 'static/networks/class_network.json','indent')

	
# run main
main()