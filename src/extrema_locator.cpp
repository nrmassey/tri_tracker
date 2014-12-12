/******************************************************************************
** Program : extrema_locator.cpp
** Author  : Neil Massey
** Date    : 10/07/13
** Purpose : class that searches for extrema in data regridded onto a 
**           regional hierarchical triangular mesh - provided by regrid
**           this class is abstract and should be inherited from.  See
**           minima_locator.h and maxima_locator.h for an example
******************************************************************************/

#include <iostream>
#include <math.h>
#include <sstream>
#include "extrema_locator.h"
#include "geo_convert.h"
#include "haversine.h"
#include "concentric_shell.h"
#include "read_from_string.h"

/*****************************************************************************/

void extrema_locator::tstep_out_begin(int t)
{
	std::cout << t;
	std::cout.flush();
}

/*****************************************************************************/

void extrema_locator::tstep_out_end(int t)
{
	int e = t;
	if (t == 0)
		std::cout << "\b";
	while (e > 0)
	{
		e = e / 10;
		std::cout << "\b";
	}
}

/*****************************************************************************/

extrema_locator::extrema_locator(void) : sv(NULL)
{
}

/*****************************************************************************/

extrema_locator::~extrema_locator(void)
{
}

/*****************************************************************************/

void extrema_locator::locate(void)
{
	find_extrema();
	find_objects();
	merge_objects();
	ex_points_from_objects();
}

/*****************************************************************************/

void extrema_locator::save(std::string output_fname, bool save_text)
{
	int n_ex = 0;
	for (int t=0; t<ex_list.size(); t++)
		n_ex += ex_list.number_of_extrema(t);
	std::cout << "# Number of extrema: " << n_ex << std::endl;
	ex_list.set_meta_data(&meta_data);

	if (sv != NULL)
		ex_list.set_meta_data(sv->get_meta_data());
	ex_list.save(output_fname, ds.get_missing_value());
	if (save_text)
	{
		std::cout << output_fname + ".txt" << std::endl;
		ex_list.save_text(output_fname+".txt", &tg);
	}
}

/*****************************************************************************/

void extrema_locator::set_steering_vector(steering_vector* isv)
{
	sv = isv;
}

/*****************************************************************************/

void extrema_locator::set_inputs(std::string input_fname, std::string mesh_fname,
								 int i_grid_level, ADJACENCY i_adj_type)
{
	// load the data
	ds_fname = input_fname;
	ds.load(input_fname);
	tg.load(mesh_fname);
	// set the other inputs
	grid_level = i_grid_level;
	adj_type = i_adj_type;
	
	// create the meta data
	std::stringstream ss;
	meta_data["mesh_file_name"] = mesh_fname;
	meta_data["input_file_name"] = input_fname;
	ss << i_grid_level;
	meta_data["extrema_grid_level"] = ss.str();
	meta_data["adjacency_type"] = i_adj_type == POINT ? "point" : "edge";
}

/*****************************************************************************/

void extrema_locator::find_extrema(void)
{
	// resize the extrema list to be the correct size for the number of timesteps
	ex_list.set_size(ds.get_number_of_time_steps());
	std::cout << "# Locating extrema, timestep: ";
	// get a list of all the triangles at the required level
	std::list<QT_TRI_NODE*> tris = tg.get_triangles_at_level(grid_level);
	// get the missing value
	FP_TYPE mv = ds.get_missing_value();
	// repeat over all timesteps
	for (int t=0; t<ds.get_number_of_time_steps(); t++)
	{
		tstep_out_begin(t);	
		// repeat over all the nodes in this level
		for (std::list<QT_TRI_NODE*>::iterator it = tris.begin();
			 it != tris.end(); it++)
		{

			indexed_force_tri_3D* c_tri = (*it)->get_data();
			if (is_extrema(c_tri, t))
			{
				// add to the extrema list, via the object list - the location and
				// value of the svex is currently just filled with the missing value
				steering_extremum svex(mv, mv, mv, mv, mv, mv);
				svex.object_labels.push_back(c_tri->get_label());
				// now add to the extrema list
				ex_list.add(t, svex);				
			}
		}
		tstep_out_end(t);
	}
	std::cout << std::endl;
}

/*****************************************************************************/

void extrema_locator::find_objects(void)
{
	// at this point the extrema list contains the labels of the triangles
	// which have extrema points.  We want to grow these points into objects
	// which encompass (for example) the entire low pressure system
	std::cout << "# Locating objects from extrema, timestep: ";
	concentric_shell c_shell;
	LABEL_STORE shell_in_object;
	for (int t=0; t<ds.get_number_of_time_steps(); t++)
	{
		tstep_out_begin(t);
		for (int e=0; e<ex_list.number_of_extrema(t); e++)
		{
			c_shell.clear();
			// get the extremum for this object
			steering_extremum* svex = ex_list.get(t, e);
			// check to see whether a label actually exists
			if (svex->object_labels.size() == 0)
				continue;
			// get the "original triangle" for this object 
			// - i.e. the one at the centre of the object
			indexed_force_tri_3D* O_TRI = tg.get_triangle(svex->object_labels[0]);
			bool keep_adding_to_obj = true;
			// continue until no more triangles are added to the object
			c_shell.calculate_inner_ring(&tg, svex);	// calculate the initial shape
			// calculate the 1st concentric shell
			c_shell.calculate(&tg, svex);
			while (keep_adding_to_obj)
			{
				keep_adding_to_obj = false;
				LABEL_STORE* c_shell_labs = c_shell.get_labels();
				shell_in_object.clear();
				// loop over all the labels
				for (LABEL_STORE::iterator it_c_shell_labs = c_shell_labs->begin();
					 it_c_shell_labs != c_shell_labs->end(); it_c_shell_labs++)
				{
					// get the candidate triangle
					indexed_force_tri_3D* C_TRI = tg.get_triangle(*it_c_shell_labs);
					if (is_in_object(O_TRI, C_TRI, t))
					{
						keep_adding_to_obj = true;
						// add to object but only if not already added
						if (std::find(svex->object_labels.begin(), svex->object_labels.end(), 
							*it_c_shell_labs) == svex->object_labels.end())
						{
							svex->object_labels.push_back(*it_c_shell_labs);
							shell_in_object.push_back(*it_c_shell_labs);
						}
					}
				}
				c_shell.recalculate(&tg, &shell_in_object);
			}
		}
		tstep_out_end(t);
	}
	std::cout << std::endl;
}

/*****************************************************************************/

bool extrema_locator::objects_share_nodes(const LABEL_STORE* o1, 
						 				  const LABEL_STORE* o2)
{
	// are two labels equal in the object lists
	// check that each object actually has a label!
	if (o1->size() == 0 || o2->size() == 0)
		return false;
	
	for (LABEL_STORE::const_iterator it_o1 = o1->begin();
		 it_o1 != o1->end(); it_o1++)
	{
		if (std::find(o2->begin(), o2->end(), *it_o1) != o2->end())
			return true;
	}
	return false;
}

/*****************************************************************************/

void extrema_locator::merge_objects(void)
{
	// merge objects together - two or more objects may have emerged from
	// extrema located close together
	std::cout << "# Merging objects, timestep: ";
	for (int t=0; t<ex_list.size(); t++)
	{
		tstep_out_begin(t);
		int o_s = ex_list.number_of_extrema(t);
		// calculate the concentric shells for each object
		std::vector<concentric_shell> obj_c_shells;
		for (int o1=0; o1<o_s; o1++)
		{
			concentric_shell o_c_shell;
			o_c_shell.calculate(&tg, ex_list.get(t, o1));
			obj_c_shells.push_back(o_c_shell);
		}
		for (int o1=0; o1<o_s; o1++)
		{			
			LABEL_STORE* o1_shell_labs = obj_c_shells[o1].get_labels();
			if (o1_shell_labs->size() == 0)	// deleted object as above
				continue;
			indexed_force_tri_3D* O1_TRI = tg.get_triangle((*o1_shell_labs)[0]);
			// get the labels in the object as well as the shell
			LABEL_STORE* o1_labs = &(ex_list.get(t, o1)->object_labels);
			for (int o2=0; o2<o_s; o2++)
			{
				if (o1 == o2)
					continue;
				LABEL_STORE* o2_shell_labs = obj_c_shells[o2].get_labels();
				if (o2_shell_labs->size() == 0)	// deleted object
					continue;
				// get the labels for both objects
				LABEL_STORE* o2_labs = &(ex_list.get(t, o2)->object_labels);
				// is a label in the shell found in the object?
				// test at different levels - 1st shells overlap?
				bool test = objects_share_nodes(o1_shell_labs, o2_shell_labs);
				// 2nd - 1st object shell overlaps with 2nd object
				if (!test) test = objects_share_nodes(o1_shell_labs, o2_labs);
				// 3rd - 2nd object shell overlaps with 1st object
				if (!test) test = objects_share_nodes(o1_labs, o2_shell_labs);
				// 4th - labels in 1st object overlaps 2nd object
				if (!test) test = objects_share_nodes(o1_labs, o2_labs);

				if (test)
				{
					// add to the first object
					for (LABEL_STORE::iterator it_o2 = o2_labs->begin(); 
						 it_o2 != o2_labs->end(); it_o2++)
					{
						// if it's not already in the object
						if (find(o1_labs->begin(), o1_labs->end(), *it_o2) == o1_labs->end())
							o1_labs->push_back(*it_o2);
					}
					obj_c_shells[o2].get_labels()->clear();	// delete the object
					o2_labs->clear();
				}
			}
			// sort mostly for debugging purposes
//			sort(o1_labs->begin(), o1_labs->end());
		}
		tstep_out_end(t);
	}
	// remove any deleted objects from the extrema list - not needed as having 
	// object_indices as size 0 marks out a deleted object / extrema
	std::cout << std::endl;
}

/*****************************************************************************/

void extrema_locator::get_min_max_values(FP_TYPE& min_v, FP_TYPE& max_v, 
										 int o, int t)
{
	// get the minimum and maximum values for an object containing a number
	// of labels
	min_v = 2e20f;
	max_v = -2e20f;
	LABEL_STORE* object_labels = &(ex_list.get(t, o)->object_labels);
	// if there are no labels do not try to find the min/max
	if (object_labels->size() == 0)
		return;
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
		if (fabs(val) < 0.99*fabs(mv) && val > max_v && fabs(val) < 1e10)
			max_v = val;
		if (fabs(val) < 0.99*fabs(mv) && val < min_v && fabs(val) < 1e10)
			min_v = val;
	}	
}

/*****************************************************************************/

void extrema_locator::calculate_steering_vector(int o, int t)
{
	// get the extremum first and the list of labels in the object
	steering_extremum* svex = ex_list.get(t, o);
	FP_TYPE mv = ds.get_missing_value();
	sv->calculate_steering_vector(&tg, svex, t, mv);
}

/*****************************************************************************/

void extrema_locator::ex_points_from_objects(void)
{
	// get the triangle from the tri_grid, via the label
	std::cout << "# Generating points from objects, timestep ";
	for (int t=0; t<ex_list.size(); t++)	// time
	{
		tstep_out_begin(t);
		for (int o=0; o<ex_list.number_of_extrema(t); o++)	// objects
		{
			calculate_object_position(o, t);
			if (sv != NULL)
				calculate_steering_vector(o, t);
			calculate_object_delta(o, t);
			calculate_object_intensity(o, t);
		}
		tstep_out_end(t);
	}
	ex_list.consolidate(ds.get_missing_value());	// remove any undefined objects	
	std::cout << std::endl;
}