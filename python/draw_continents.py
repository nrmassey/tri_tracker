
import sys,os, pickle
from matplotlib.path import Path
import matplotlib.patches as patches

###############################################################################

def read_continents():
	# check first whether the continents have already been pickled
	path = "/Users/massey/Coding/cpdn_analysis/map_plot/data/"
	pik_name = path + "map_continents.pik"
	if os.path.exists(pik_name):
		fh = open(pik_name)
		continents = pickle.load(fh)
		fh.close()
	else:
		# read in the continents file
		fh = open(path + "map_continents_lowres.txt")
		lines = fh.readlines()			
		fh.close()
		continents = []
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
				# check for wrap around date line
				if abs(pt[0] - cont[-1][0]) < 5.0:
					cont.append(pt)
				else:
					# create new continent - special case for europe / africa and antartica
					if c_cont == 0:
						# draw a box to the new poitn
						cont.append([-180,-90])
						cont.append([180,-90])
						cont.append(pt)
					elif c_cont == 186:
						phase += 1
						if phase == 1:
							continent_points.append(cont)
							cont = [pt]
						elif phase == 2:
							continent_points.append(cont)
							cont = continent_points[-2]
							cont.append(pt)

			c_cont += 1
			continent_points.append(cont)
		
		# create the path commands for the continent
		for cont in continent_points:
			cont_codes = [Path.MOVETO]
			for pt_n in range(1, len(cont)):
				cont_codes.append(Path.LINETO)
					
			# create a path from cont and cont_path
			cont_path = Path(cont, cont_codes)				
			continents.append(cont_path)
		# write the data out as a pickle so we only have to construct it once
		fo = open(pik_name, "w")
		pickle.dump(continents, fo)
		fo.close()
		
	return continents

###############################################################################
		
def draw_continents(sp, fillcolor='none', cn=0, lw=0.5, alpha=1.0):
	"""draw the continents on the map"""
	res=1
	if fillcolor == 'none':
		ec = 'k'
	else:
		ec = fillcolor
	continents = read_continents()
	cont = continents[cn]
	for cont in continents:
		cont_p = patches.PathPatch(cont, edgecolor=ec, facecolor=fillcolor, lw=lw, zorder=1,
								   alpha=alpha)
		p = sp.add_patch(cont_p)
