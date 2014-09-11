/******************************************************************************
** Program : regridder.h
** Author  : Neil Massey
** Date    : 28/06/13
** Purpose : class that takes a regional triangular grid and a netCDF file
**           the netCDF file is regridded onto the mesh
******************************************************************************/

#ifndef REGRIDDER_H
#define REGRIDDER_H

#include <string>
#include "tri_grid.h"
#include "ncdata.h"
#include "data_store.h"

/*****************************************************************************/

class regridder
{
	public:
		regridder(const std::string mesh_fname, const std::string nc_fname, 
				  const std::string nc_vname, const int z_level,
				  FP_TYPE smooth_weight, const int p_method);
		~regridder(void);
		void regrid(void);
		void save(std::string out_fname);
		void save_text(std::string out_fname);

	private:
		
		void smooth_data(void);
	
		/*********************************************************************/
		
		tri_grid tg;
		ncdata nc_data;
		data_store* ds;
		int z_level;
		FP_TYPE smooth_weight;
		FP_TYPE p_method;
		
		/*********************************************************************/
};

#endif
