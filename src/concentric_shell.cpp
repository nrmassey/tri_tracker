/******************************************************************************
** Program : concentric_shell.cpp
** Author  : Neil Massey
** Date    : 07/10/13
** Purpose : class that creates a concentric shell from an object list
**           the concentric shell is just the triangles that surround the object
******************************************************************************/

#include "concentric_shell.h"
#include <string>
#include <list>
#include <algorithm>

const ADJACENCY adjacency_type = EDGE;

/*****************************************************************************/

concentric_shell::concentric_shell(void)
{
	labels.clear();
}

/*****************************************************************************/

void concentric_shell::calculate(tri_grid* tg, extremum* ex)
{
	// calculate the concentric shell
	labels.clear();
	// get the labels from the object
	LABEL_STORE& obj_labels = ex->object_labels;
	// loop over them
	for (LABEL_STORE::iterator it_obj_lab = obj_labels.begin();
		 it_obj_lab != obj_labels.end(); it_obj_lab++)
	{
		// get the triangle from the label
		indexed_force_tri_3D* obj_tri = tg->get_triangle(*it_obj_lab);
		// get the labels of the adjacent triangles
		const LABEL_STORE* obj_adj_labels = obj_tri->get_adjacent_labels(adjacency_type);
		// loop through and check that the labels are either
		// 1. not in the object
		// 2. not already in the concentric shell
		for (LABEL_STORE::const_iterator it_adj_lab = obj_adj_labels->begin();
			 it_adj_lab != obj_adj_labels->end(); it_adj_lab++)
		{
			if (std::find(labels.begin(), labels.end(), *it_adj_lab) == labels.end() &&
				std::find(obj_labels.begin(), obj_labels.end(), *it_adj_lab) == obj_labels.end())
			{
				labels.push_back(*it_adj_lab);
			}
		}
	}
}

/*****************************************************************************/

void concentric_shell::calculate_inner_ring(tri_grid* tg, extremum *ex)
{
	// The inner ring (outer ring of an object) is defined as having less than 
	// four edge adjacent triangles that are within the object itself
	inner_ring.clear();
	// get the labels from the object
	LABEL_STORE& obj_labels = ex->object_labels;
	// loop over them
	for (LABEL_STORE::iterator it_obj_lab = obj_labels.begin();
		 it_obj_lab != obj_labels.end(); it_obj_lab++)
	{
		// get the triangle from the label
		indexed_force_tri_3D* obj_tri = tg->get_triangle(*it_obj_lab);
		// get the labels of the adjacent triangles
		const LABEL_STORE* obj_adj_labels = obj_tri->get_adjacent_labels(adjacency_type);
		int label_count = 0;
		for (LABEL_STORE::const_iterator it_adj_lab = obj_adj_labels->begin();
			 it_adj_lab != obj_adj_labels->end(); it_adj_lab++)
		{
			if(std::find(obj_labels.begin(), obj_labels.end(), *it_adj_lab) != obj_labels.end())
				label_count ++;			
		}
		if (label_count < 3)
			inner_ring.push_back(*it_obj_lab);
	}
}
/*****************************************************************************/

void concentric_shell::recalculate(tri_grid* tg, LABEL_STORE* shell_in_object)
{
	// shell_in_object contains all the labels which the feature detection
	// algorithm deems to be in the object
	// take the triangle labels (edge) adjacent to these labels and add to
	// the new shell if they are not already in the "inner ring"
	// then add the shell in object to the inner ring
	
	LABEL_STORE new_shell;
	// loop over every label in the shell_in_object
	for (LABEL_STORE::iterator it_sio = shell_in_object->begin();
		 it_sio != shell_in_object->end(); it_sio ++)
	{
		// get the adjacent labels
		indexed_force_tri_3D* sio_tri = tg->get_triangle(*it_sio);
		const LABEL_STORE* sio_adj_labels = sio_tri->get_adjacent_labels(adjacency_type);
		// loop over the adjacent labels
		for (LABEL_STORE::const_iterator it_sio_adj = sio_adj_labels->begin();
			 it_sio_adj != sio_adj_labels->end(); it_sio_adj++)
		{
			// if this label is not in the inner ring add to the new shell
			if (std::find(inner_ring.begin(), inner_ring.end(), *it_sio_adj) == inner_ring.end())
				new_shell.push_back(*it_sio_adj);
		}
		// add the label to the inner ring
		inner_ring.push_back(*it_sio);
	}
	// assign the labels to be the new shell
	labels = new_shell;
}

/*****************************************************************************/

LABEL_STORE* concentric_shell::get_inner_ring(void)
{
	return &inner_ring;
}

/*****************************************************************************/

LABEL_STORE* concentric_shell::get_labels(void)
{
	return &labels;
}

/*****************************************************************************/

void concentric_shell::clear(void)
{
	labels.clear();
}