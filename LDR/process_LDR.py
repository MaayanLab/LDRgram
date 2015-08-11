# I will process the LDR data with this script 

def main():

	# # initial 
	# ##############
	# # first look at the json 
	# first_look()
	# # construct 3d array
	# example_construct_array()

	# extract data and save to flattened matrices 
	###############################################
	# extract the nodes
	nodes = extract_nodes()
	# construct actual array
	construct_array(nodes)

	# # collapse into three matrices 
	# collapse_mat()

	# # export tsv for clustergrammer visualization 
	# ##############################################
	# export_tsv()




def construct_array(nodes):
	import json_scripts
	import scipy

	print('\nconstructing array')

	# load the LDR data is json format 
	ldr = json_scripts.load_to_dict('LDR_api.json')

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

				print('\n\n')

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


def example_construct_array():
	import scipy
	print('constructing example array')


	mat = scipy.zeros([3, 5	, 5])

	print(mat)

	# the ndim array is basically 
	# three 5x5 matrices stacked on each other 

	# [[[ 0.  0.  0.  0.  0.]
	#  [ 0.  0.  0.  0.  0.]
	#  [ 0.  0.  0.  0.  0.]
	#  [ 0.  0.  0.  0.  0.]
	#  [ 0.  0.  0.  0.  0.]]

	# [[ 0.  0.  0.  0.  0.]
	#  [ 0.  0.  0.  0.  0.]
	#  [ 0.  0.  0.  0.  0.]
	#  [ 0.  0.  0.  0.  0.]
	#  [ 0.  0.  0.  0.  0.]]

	# [[ 0.  0.  0.  0.  0.]
	#  [ 0.  0.  0.  0.  0.]
	#  [ 0.  0.  0.  0.  0.]
	#  [ 0.  0.  0.  0.  0.]
	#  [ 0.  0.  0.  0.  0.]]]

	print('\n\n')

	# printing one slice of the matrix 
	print(mat[0,:,:])
	
	# [[ 0.  0.  0.  0.  0.]
	# [ 0.  0.  0.  0.  0.]
	# [ 0.  0.  0.  0.  0.]
	# [ 0.  0.  0.  0.  0.]
	# [ 0.  0.  0.  0.  0.]]

def first_look():
	import json_scripts

	print('having a first look at the data')

	# load the LDR data is json format 
	ldr = json_scripts.load_to_dict('LDR_api.json')

	# the data is a list of jsons, each json has information on a study 

	print('Assay: datasetName')
	print(ldr[0]['datasetName'])
	print('\n')

	print('CellLines')
	print(ldr[0]['metadata']['cellLines'][0]['name'])
	print('\n')

	print('perturbagens')
	print(ldr[0]['metadata']['perturbagens'][0]['name'])
	print('\n')

def extract_nodes():
	import json_scripts

	print('extracting nodes: as, cl, pt')

	# load the LDR data is json format 
	ldr = json_scripts.load_to_dict('LDR_api.json')

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



def export_tsv():
	import json_scripts
	import numpy as np

	ldr_mat = json_scripts.load_to_dict('ldr_mat_collapse.json')

	ldr_mat['as_cl_mat'] = np.asarray(ldr_mat['as_cl_mat'])

	print(ldr_mat.keys())

	print('as\t'+str(objectlen(ldr_mat['nodes']['as'])))
	print('cl\t'+str(len(ldr_mat['nodes']['cl'])))
	print('pt\t'+str(len(ldr_mat['nodes']['pt'])))

	# write tab separated file 
	filename = 'tmp_export.txt'
	# write first line 
	fw = open(filename, 'w')

	fw.write('\t')

	# write row
	for i in range(len(ldr_mat['nodes']['cl'])):

		fw.write( ldr_mat['nodes']['cl'][i] + '\t')

	fw.write('\n')


	# get matrix dimensions 
	inst_dim = ldr_mat['as_cl_mat']

	x = ldr_mat['as_cl_mat'].shape[0]
	y = ldr_mat['as_cl_mat'].shape[1]

	# print(inst_dim)

	for i in range(x):

		# print rows 
		fw.write( ldr_mat['nodes']['as'][i] + '\t' )

		# print columns 
		for j in range(y):

			# write cols 
			fw.write( str(ldr_mat['as_cl_mat'][i,j]) + '\t' )

		fw.write('\n')

	fw.close()

def collapse_mat():
	import json_scripts
	import scipy
	import numpy as np

	# load ldr_mat
	ldr_mat = json_scripts.load_to_dict('ldr_mat.json')

	# save matrix as array 
	ldr_mat['mat'] = np.asarray(ldr_mat['mat'])

	# matrix
	# mat[ as, cl, pt ]

	# make three matrices 
	# as vs cl, sum pt
	# as vs pt, sum cl 
	# cl vs pt, sum as 

	# make list of lists
	comp = []
	comp.append(['as','cl'])
	comp.append(['as','pt'])
	comp.append(['cl','pt'])

	flat_index = [2,1,0]

	print(ldr_mat['mat'].shape)
	print('as\t'+str(len(ldr_mat['nodes']['as'])))
	print('cl\t'+str(len(ldr_mat['nodes']['cl'])))
	print('pt\t'+str(len(ldr_mat['nodes']['pt'])))
	print('\n')

	# initialize three matrices 
	for i in range(len(comp)):

		# get inst_comp
		inst_comp = comp[i]

		# get dimensions 
		inst_x = len(ldr_mat['nodes'][inst_comp[0]])
		inst_y = len(ldr_mat['nodes'][inst_comp[1]])

		print('\n')
		print(flat_index[i])
		print(inst_comp)
		print( np.sum(ldr_mat['mat'], axis=flat_index[i]).shape )

		# generate collapsed matrices 
		ldr_mat[inst_comp[0]+'_'+inst_comp[1]+'_mat'] = np.sum(ldr_mat['mat'], axis=flat_index[i])

		# convert matrix to list 
		ldr_mat[inst_comp[0]+'_'+inst_comp[1]+'_mat'] = ldr_mat[inst_comp[0]+'_'+inst_comp[1]+'_mat'].tolist()

	# convert ldr_mat['mat'] to list
	ldr_mat['mat'] = ldr_mat['mat'].tolist() 

	# save to json 
	json_scripts.save_to_json( ldr_mat, 'ldr_mat_collapse.json', 'no-indent' )




# run main
main()