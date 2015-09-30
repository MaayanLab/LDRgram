def main():
	import json
	from d3_clustergram_class import Network

	net = Network()

	mat_info = {}

	mat_info[str((1,1))] = 1

	print(mat_info[ str((1,1)) ])	


	print(mat_info)	
	print(type(mat_info))

	tmp = json.dumps(mat_info)

	print(tmp)
	print(type(tmp))

main()