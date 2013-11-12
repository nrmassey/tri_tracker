/******************************************************************************
** Program : minima_background.cpp
** Author  : Neil Massey
** Date    : 07/08/13
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded, and removing the background field
******************************************************************************/

#include "minima_background.h"
#include "haversine.h"
#include "geo_convert.h"
#include <sstream>
#include <math.h>

/******************************************************************************/

minima_background::minima_background(void) : extrema_locator(), bck_field_ds(NULL)
{
}

/******************************************************************************/

minima_background::~minima_background(void)
{
	delete bck_field_ds;
}

/******************************************************************************/

void minima_background::locate(void)
{
	process_data();
	find_extrema();
	if (grid_level > 5)
	{
		merge_objects();
		refine_objects();
	}
	find_objects();
	if (grid_level > 5)
		trim_objects();
	merge_objects();
	ex_points_from_objects();	
}

/******************************************************************************/

void minima_background::parse_arg_string(std::string method_string)
{
	// arguments are:
	// arg[0] = file to read background
	// arg[1] = level to use as the background field
	// arg[2] = averaging period to take background field over
	// arg[3] = mesh level to detect large scale minima at
	// arg[4] = contour value
	// parameters for minima location with background removal
	// get the first bracket
	int c_pos = method_string.find_first_of("(")+1;
	int e_pos = method_string.find(")", c_pos);
	char dummy;
	if (method_string.substr(c_pos, e_pos-c_pos) == "help")
		throw(std::string("minima_back parameters = (file name to take background field from, level in mesh to use as background field, averaging period of background field, contour value)"));
	
	int b_pos = method_string.find(",", c_pos);
	bck_field_file = method_string.substr(c_pos, b_pos-c_pos);
	std::stringstream stream(method_string.substr(b_pos+1, e_pos-b_pos));
	stream >> bck_mesh_lvl >> dummy 
		   >> bck_avg_period >> dummy
		   >> min_mesh_lvl >> dummy
		   >> contour_value >> dummy
		   >> min_delta;
		   
	// add to the metadata
	std::stringstream ss;
	meta_data["method"] = "minima_back";
	meta_data["background_file"] = bck_field_file;
	ss << bck_mesh_lvl;
	meta_data["background_mesh_level"] = ss.str();
	ss.str(""); ss << bck_avg_period;
	meta_data["background_averaging_period"] = ss.str();
	ss.str(""); ss << min_mesh_lvl;
	meta_data["minima_mesh_level"] = ss.str();
	ss.str(""); ss << contour_value;
	meta_data["contour_value"] = ss.str();
	ss.str(""); ss << min_delta;
	meta_data["minimum_delta"] = ss.str();
}

/******************************************************************************/

void minima_background::get_child_values(QT_TRI_NODE* c_node, int max_depth, 
										 int t_step, 
										 std::vector<FP_TYPE>& child_values)
{
	// walk the tree to either a leaf node or the max_depth and add the
	// value to the vector
	if (c_node->is_leaf() or c_node->get_level()==max_depth)
		child_values.push_back(ds.get_data(t_step, c_node->get_data()->get_ds_index()));
	else
		for (int n=0; n<4; n++)
			get_child_values(c_node->get_child(n), max_depth, t_step, child_values);
}

/******************************************************************************/

bool minima_background::is_extrema(indexed_force_tri_3D* tri, int t_step)
{
	FP_TYPE tri_val = ds.get_data(t_step, tri->get_ds_index());
	// if it's the missing value then return false
	if (fabs(tri_val) >= fabs(0.99*ds.get_missing_value()))
		return false;
	if (tri_val > -min_delta)	// must be a negative anomaly
		return false;

	// find the minima in the triangle n grid levels up and its value
	QT_TRI_NODE* gf_tri_node = tg.get_triangle_node(tri->get_label());
	// go up the tree to the desired level
	int nups = grid_level - min_mesh_lvl;
	for (int l=0; l<nups; l++)
		gf_tri_node = gf_tri_node->get_parent();
	// get the triangle
	indexed_force_tri_3D* gf_tri = gf_tri_node->get_data();	
	// get the data
	FP_TYPE gf_tri_val = ds.get_data(t_step, gf_tri->get_ds_index());
	// if it's the missing value then return false
	if (fabs(gf_tri_val) >= fabs(0.99*ds.get_missing_value()))
		return false;
	// if it is not negative then return false
	if (gf_tri_val > -min_delta)
		return false;
		
	// record the number of surrounding triangles that are greater than current triangle
	FP_TYPE n_st = 0;	
	// loop through all the adjacent triangles
	const LABEL_STORE* gf_adj_tri_labels = gf_tri->get_adjacent_labels(adj_type);
	// nearest neighbour search of all adjacent triangles for a minimum
	for (LABEL_STORE::const_iterator gf_adj_it = gf_adj_tri_labels->begin();
		 gf_adj_it != gf_adj_tri_labels->end(); gf_adj_it++)
	{
		// get the triangle from the label
		indexed_force_tri_3D* gf_adj_tri = tg.get_triangle(*gf_adj_it);
		// get the value of the adjacent triangle		
		FP_TYPE gf_adj_val = ds.get_data(t_step, gf_adj_tri->get_ds_index());
		// if it's the missing value then continue onto next one
		if (fabs(gf_adj_val) >= fabs(0.99*ds.get_missing_value()))
			continue;
		// if the middle triangle is less than or equal to this surrounding triangle
		FP_TYPE V = (gf_tri_val - gf_adj_val);
		if (V <= 0.0)
			n_st += 1.0;
	}
	int min_sur = adj_type == POINT ? 11 : 4;
	return (n_st >= min_sur);
}

/******************************************************************************/

bool minima_background::is_in_object(indexed_force_tri_3D* O_TRI, 
				  					 indexed_force_tri_3D* C_TRI, int t_step)
{
	// O_TRI - original triangle
	// C_TRI - candidate triangle - triangle being tested for inclusion
	bool is_in = false;

	// get the candidate triangle value
	FP_TYPE cl_v = ds.get_data(t_step, C_TRI->get_ds_index());
	FP_TYPE ol_v = ds.get_data(t_step, O_TRI->get_ds_index());

	// quick check	
	is_in = cl_v <= ol_v;
	// not the mv
	is_in &= fabs(cl_v) <= fabs(0.99 * ds.get_missing_value());
	is_in &= (cl_v <= -min_delta);
	return is_in;
}

/******************************************************************************/

void minima_background::calculate_background_field(void)
{
	std::cout << "# Calculating background field" << std::endl;
	
	// load in the field first
	bck_field_ds = new data_store();
	bck_field_ds->load(bck_field_file);
	
	// get the size of the current datastore
	int n_ts = bck_field_ds->get_number_of_time_steps();
	int n_idx = bck_field_ds->get_number_of_indices();
	
	// do we need to take a mean?
	if (bck_avg_period > 1 && n_ts > bck_avg_period)
	{
		// create a new background field
		data_store* new_bck_field_ds = new data_store();
		// set the scaling for each averaging period
		FP_TYPE scale = 1.0 / bck_avg_period;
		// create the datastore
		new_bck_field_ds->set_size(n_ts/bck_avg_period, n_idx);
		new_bck_field_ds->set_missing_value(bck_field_ds->get_missing_value());
		// loop through the data producing the average
		for (int t=0; t<n_ts; t++)
		{
			int dest_pos = t / bck_avg_period;
			for (int i=0; i<n_idx;  i++)
			{
				// get the current data value, add the value from the original ds store
				// multiplied by the scaler
				FP_TYPE c_val = new_bck_field_ds->get_data(dest_pos, i);
				FP_TYPE t_val = scale * bck_field_ds->get_data(t, i) + c_val;
				// contour the data
				t_val = float(int(t_val/contour_value)+0.5) * contour_value;
				new_bck_field_ds->set_data(dest_pos, i, t_val);
			}
		}
		// assign the bck_field to be the new_field and delete the old one
		delete bck_field_ds;
		bck_field_ds = new_bck_field_ds;
	}
}

/******************************************************************************/

void minima_background::contour_data(void)
{
	// contour the data - i.e. round it to intervals of contour_value
	for (int t=0; t<ds.get_number_of_time_steps(); t++)
		for (int i=0; i<ds.get_number_of_indices(); i++)
		{
			FP_TYPE val = ds.get_data(t, i);
			val = float(int(val/contour_value)+0.5) * contour_value;
			ds.set_data(t, i, val);
		}
}

/******************************************************************************/

bool minima_background::process_data(void)
{
	// process the data
	// this removes the background field from the data
	// the background field may be averaged over a time period - i.e. remove
	// the monthly average / 5 day average or just daily average
	// the difference is taken between the triangle and an ancestor triangle
	// (remember that data is propagated up the tree to the parent and ancestors)
	// the background to be subtracted is taken to be a weighted sum of all the ancestor's
	// surrounding triangles (using the adjacency list) based on the distance between them
	// and the original triangle.  this is to avoid edge effects, where a number of minima
	// are identified along the edge of a triangle because the difference at that edge is
	// greater than elsewhere in the triangle
	
	// first calculate the background field
	calculate_background_field();
	
	std::cout << "# Processing data" << std::endl;
	int n_ts = ds.get_number_of_time_steps();
	FP_TYPE mv = ds.get_missing_value();

	// we want to process every level, not just the grid level
	for (int c_level = 0; c_level <= grid_level; c_level ++)
	{	
		// get a list of all the triangles at the required level
		std::list<QT_TRI_NODE*> tris = tg.get_triangles_at_level(c_level);
		// how many times do we have to go up the tree to get to the required parent?
		int n_ups = c_level - bck_mesh_lvl;
		for (std::list<QT_TRI_NODE*>::iterator it = tris.begin();
			 it != tris.end(); it++)		 
		{
			// get the index of the current triangle
			int o_idx = (*it)->get_data()->get_ds_index();
			// get the position of the triangle in the original grid
			vector_3D o_tri_centroid = (*it)->get_data()->centroid();
			FP_TYPE o_lon, o_lat;
			cart_to_model(o_tri_centroid, o_lon, o_lat);
		
			// get the index of the ancestor triangle
			QT_TRI_NODE* anc_tri = (*it);
			for (int p=0; p<n_ups; p++)
				// ascend the tree
				anc_tri = anc_tri->get_parent();
			
			// indices and weights for the ancestor and surrounding triangles
			std::vector<FP_TYPE> weights;
			std::vector<int> ds_idxs;
		
			// get the lat, lon of the ancestor triangle
			FP_TYPE s_lon, s_lat, dist;
			cart_to_model(anc_tri->get_data()->centroid(), s_lon, s_lat);
			dist = haversine(o_lon, o_lat, s_lon, s_lat, 1.0);
			weights.push_back(dist);
			ds_idxs.push_back(anc_tri->get_data()->get_ds_index());
			// want to find the maximum distance
			FP_TYPE max_dist = dist;
		
			// get the list of adjacent triangles
			const LABEL_STORE* anc_adj_list = anc_tri->get_data()->get_adjacent_labels(POINT);
			for (LABEL_STORE::const_iterator it_anc_tris = anc_adj_list->begin();
				 it_anc_tris != anc_adj_list->end(); it_anc_tris++)
			{
				// get a surrounding triangle
				indexed_force_tri_3D* anc_sur_tri = tg.get_triangle(*it_anc_tris);
				// get the centroid and convert to lat / lon
				cart_to_model(anc_sur_tri->centroid(), s_lon, s_lat);
				// calculate the distance between the original triangle and this triangle
				dist = haversine(o_lon, o_lat, s_lon, s_lat, 1.0);
				if (dist > max_dist)
					max_dist = dist;
				// add the index and distance to the above lists
				weights.push_back(dist);
				ds_idxs.push_back(anc_sur_tri->get_ds_index());
			}
			// the "weights" at this point contain just the distances - normalize by the max_dist and
			// set the nearest to have the highest weight
			for (unsigned int w=0; w<weights.size(); w++)
				weights[w] = 1.0 - weights[w]/max_dist;
			// loop over the time and transform to the background value
			int bck_field_nts = bck_field_ds->get_number_of_time_steps();
			for (int t=0; t<n_ts; t++)
			{
				FP_TYPE o_val = ds.get_data(t, o_idx);			// get the data of the original triangle
				// get the sum of the surrounding triangle values * weights, and the sum of the weights
				float sum_v = 0.0, sum_w = 0.0;
				for (unsigned int w=0; w<weights.size(); w++)
				{
					// check for stepping outside of the array
					int bck_t = t/bck_avg_period;
					if (bck_t >= bck_field_nts)
						bck_t = bck_field_nts-1;
					// get the background value
					FP_TYPE val = bck_field_ds->get_data(bck_t, ds_idxs[w]);
					if (fabs(val) < 0.99*fabs(mv))				// numerical inaccuracy fix!
					{
						sum_v += weights[w] * val;
						sum_w += weights[w];
					}
				}
				FP_TYPE t_val;
				if (fabs(o_val) > 0.99*fabs(mv) || sum_w == 0.0)
					t_val = mv;
				else
					t_val = o_val - sum_v / sum_w;				// do the difference
				ds.set_data(t, o_idx, t_val);					// write the data back
			}
		}
	}
	// smooth the data
	smooth_data();
	// contour the data
	if (contour_value != 1.0)
		contour_data();
//	ds.save("ds_background_removed.rgd");
	return true;
}
/******************************************************************************/

void minima_background::smooth_data(void)
{
	// smooth the data as the mean average of the surrounding triangles
	FP_TYPE mv = ds.get_missing_value();
	// get a list of all the triangles at the required level
	std::list<QT_TRI_NODE*> tris = tg.get_triangles_at_level(grid_level);
	// require a new location - create a datastore
	data_store smooth_ds;
	smooth_ds = ds;
	// loop over each triangle
	for (std::list<QT_TRI_NODE*>::iterator it = tris.begin(); it != tris.end(); it++)
	{
		// a list to keep the indices of the surrounding triangles of the current triangle
		std::vector<int> ds_idxs;
		// get the location in the data store of the current triangle
		int o_idx = (*it)->get_data()->get_ds_index();
		ds_idxs.push_back(o_idx);
		// build the rest of the vector of indices
		const LABEL_STORE* anc_adj_list = (*it)->get_data()->get_adjacent_labels(POINT);
		for (LABEL_STORE::const_iterator it_anc_tris = anc_adj_list->begin();
			 it_anc_tris != anc_adj_list->end(); it_anc_tris++)
		{
			// get a surrounding triangle
			indexed_force_tri_3D* anc_sur_tri = tg.get_triangle(*it_anc_tris);
			ds_idxs.push_back(anc_sur_tri->get_ds_index());
		}				
		// loop over the time info
		for (int t=0; t<ds.get_number_of_time_steps(); t++)
		{
			FP_TYPE val = 0.0;
			int n_vals = 0;
			for (int d=0; d<ds_idxs.size(); d++)
			{
				FP_TYPE c_val = ds.get_data(t, ds_idxs[d]);
				if (fabs(c_val) < 0.99*fabs(mv))
				{
					val += c_val;
					n_vals++;
				}
			}
			if (n_vals > 0)
				smooth_ds.set_data(t, o_idx, val/n_vals);
			else
				smooth_ds.set_data(t, o_idx, mv);			
		}
	}
	ds = smooth_ds;
}

/******************************************************************************/

FP_TYPE minima_background::calculate_point_weight(FP_TYPE V, FP_TYPE min_v, FP_TYPE max_v)
{
	FP_TYPE w = 1.0 - (V-min_v) / (max_v-min_v);
	w = w * w;	// square to penalise the values further away from the min_v
	return w;
}

/******************************************************************************/

void minima_background::trim_objects(void)
{
	std::cout << "# Trimming objects, timestep: ";
	FP_TYPE min_diameter = 50;
	FP_TYPE max_diameter = 1500;

	// calculate surface area of a triangle
	LABEL tri_lab_0 = ex_list.get(0,0)->object_labels[0];
	indexed_force_tri_3D* tri_0 = tg.get_triangle(tri_lab_0);
	FP_TYPE tri_surf_area = tri_0->surface_area(6371);

	for (int t=0; t<ex_list.size(); t++)
	{
		tstep_out_begin(t);
		int o_s = ex_list.number_of_extrema(t);
		for (int o1=0; o1<o_s; o1++)
		{
			LABEL_STORE* o1_labs = &(ex_list.get(t, o1)->object_labels);
			FP_TYPE object_area = o1_labs->size() * tri_surf_area;
			FP_TYPE object_diameter = sqrt(object_area / M_PI);
			if (object_diameter < min_diameter || object_diameter > max_diameter)
			{
				o1_labs->clear();		// delete!
			}
		}
		tstep_out_end(t);	
	}
	std::cout << std::endl;
}

/******************************************************************************/

void minima_background::refine_objects(void)
{
	// refine the object (after the initial merge) so that any triangle which is not
	// the minimum value is removed from this initial object
	std::cout << "# Refining objects, timestep: ";
	
	for (int t=0; t<ex_list.size(); t++)
	{
		tstep_out_begin(t);
		int o_s = ex_list.number_of_extrema(t);
		for (int o1=0; o1<o_s; o1++)
		{
			FP_TYPE min_v, max_v;
			get_min_max_values(min_v, max_v, o1, t);
			LABEL_STORE* object_labels = &(ex_list.get(t, o1)->object_labels);
			LABEL_STORE new_labels;
			for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
				 it_ll != object_labels->end(); it_ll++)	// tri indices
			{
				indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
				FP_TYPE val = ds.get_data(t, c_tri->get_ds_index());
				bool add_label = false;
				if (val <= min_v + contour_value)
					add_label = true;					
				if (add_label)
					new_labels.push_back(*it_ll);
			}
			ex_list.get(t, o1)->object_labels = new_labels;
		}
		tstep_out_end(t);		
	}
	std::cout << std::endl;
}

/******************************************************************************/

indexed_force_tri_3D* minima_background::get_original_triangle(int o, int t)
{
	// get the original triangle - i.e. the one at the centre of the object
	// for the minima this is defined as the triangle with the lowest value
	FP_TYPE c_val = 2e20f;
	indexed_force_tri_3D* o_tri = NULL;
	LABEL_STORE* object_labels = &(ex_list.get(t, o)->object_labels);
	// if there are no labels do not try to find the min/max
	if (object_labels->size() == 0)
		return NULL;
	// get the missing value
	FP_TYPE mv = ds.get_missing_value();
	// loop over all the triangles in the object
	for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
		 it_ll != object_labels->end(); it_ll++)	// tri indices
	{
		// get the triangle
		indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
		FP_TYPE val = ds.get_data(t, c_tri->get_ds_index());
		// find the min and max values
		if (fabs(val) < 0.99*fabs(mv) && val < c_val)
		{
			c_val = val;
			o_tri = c_tri;
		}
	}
	return o_tri;
}