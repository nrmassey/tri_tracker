from plot_all_extrema import calc_bearing
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import math

def curvature_cost(P, Q, C):
	# calculate cost of appending C to track containing P and Q
	b1 = calc_bearing(P[0],P[1], Q[0],Q[1])
	b2 = calc_bearing(Q[0],Q[1], C[0],C[1])
	c = abs(b2-b1)
	if (c > 180.0):
		c -= 90.0
	return c / 180.0
	

track_points = [[30,30],[31,29]] # // longitude / latitude of points already in track
d=5
a=1.0/math.sqrt(2)*d
tp = track_points[1]
candidate_points = [[tp[0],tp[1]+d],[tp[0]+a,tp[1]+a],[tp[0]+d,tp[1]],[tp[0]+a,tp[1]-a],
                    [tp[0],tp[1]-d],[tp[0]-a,tp[1]-a],[tp[0]-d,tp[1]],[tp[0]-a,tp[1]+a]]

projection = ccrs.PlateCarree()
sp0 = plt.subplot(111, projection=projection)

for tp in track_points:
	sp0.plot(tp[0], tp[1], marker='o', mec='r', mfc='r', ms=3.5)
	
sp0.plot([track_points[0][0], track_points[1][0]], 
		 [track_points[0][1], track_points[1][1]], 'r')
		 
for tp in candidate_points:
	sp0.plot(tp[0], tp[1], marker='o', mec='b', mfc='b', ms=3.5)
	c = curvature_cost(track_points[0], track_points[1], tp)
	sp0.text(tp[0], tp[1]-0.5, str(c))
	
#sp0.set_xlim([track_points[0][0]-1, track_points[1][0]+d+1])
#sp0.set_ylim([track_points[0][1]-d-1, track_points[1][1]+d+1])

plt.show()