
def read_continents():
	# check first whether the continents have already been pickled
	path = "../"
	# read in the continents file
	fh = open(path + "map_continents_lowres.txt")
	lines = fh.readlines()			
	fh.close()
	# get the info line by line - each continent takes up 3 lines
	continent_points = []
	c_cont = 0
	for i in range(0, len(lines), 3):
		if (lines[i][0] != "#"):
			print "error in file"
		lons = lines[i+1].split(",")
		lats = lines[i+2].split(",")
		# build a list of the continents lat & lon coords
		cont = [[float(lons[0]), float(lats[0])]]
		phase = 0
		for j in range(1, len(lons)):
			pt = [float(lons[j]), float(lats[j])]
			cont.append(pt)
		c_cont += 1
		continent_points.append(cont)
	return continent_points