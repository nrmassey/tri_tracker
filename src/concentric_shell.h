/******************************************************************************
** Program : concentric_shell.h
** Author  : Neil Massey
** Date    : 07/10/13
** Purpose : class that creates a concentric shell from an object list
**           the concentric shell is just the triangles that surround the object
******************************************************************************/

#ifndef CONCENTRIC_SHELL_H
#define CONCENTRIC_SHELL_H

#include "tri_grid.h"
#include "extremum.h"

class concentric_shell
{
	public:
		concentric_shell(void);
		void calculate(tri_grid* tg, extremum* ex);
		void recalculate(tri_grid* tg, LABEL_STORE* shell_in_object);
		LABEL_STORE* get_labels(void);
		void clear(void);
		void calculate_inner_ring(tri_grid* tg, extremum *ex);
		LABEL_STORE* get_inner_ring(void);
		
	private:
		LABEL_STORE labels;
		LABEL_STORE inner_ring;	// this is the outer ring of the object
};

#endif