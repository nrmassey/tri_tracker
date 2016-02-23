/******************************************************************************
** Program : extremum.h
** Author  : Neil Massey
** Date    : 13/07/13
** Purpose : class to hold extremum data and function manipulations
******************************************************************************/

#include <fstream>
#include "extremum.h"
#include "bin_file_utils.h"

/*****************************************************************************/

void steering_extremum::set(FP_TYPE ilon, FP_TYPE ilat, FP_TYPE iintensity, FP_TYPE idelta,
	 				        FP_TYPE iu, FP_TYPE iv) 
{
	lon = ilon; lat = ilat;
	intensity = iintensity; 
	delta = idelta;
	sv_u = iu; sv_v = iv;
	deleted = false;
}

/*****************************************************************************/

void steering_extremum::set(FP_TYPE ilon, FP_TYPE ilat, FP_TYPE iintensity, FP_TYPE idelta,
     				        FP_TYPE iu, FP_TYPE iv, LABEL_STORE iobject_labels)
{
	set(ilon, ilat, iintensity, idelta, iu, iv);
	deleted = false;
	object_labels.clear();
	for (LABEL_STORE::iterator it_obj_lab = iobject_labels.begin();
		 it_obj_lab != iobject_labels.end(); it_obj_lab++)
    	object_labels.push_back(*it_obj_lab);
}

/*****************************************************************************/

void steering_extremum::save(std::ofstream& out)
{
	// write the extremum out as a binary chunk
	// write the timestep, longitude, latitude, intensity, delta, geostrophic wind u & v
	write_float(out, lon);
	write_float(out, lat);
	write_float(out, intensity);
	write_float(out, delta);
	write_float(out, sv_u);
	write_float(out, sv_v);
	// write number of object labels first
	write_int(out, object_labels.size());
	for (LABEL_STORE::iterator it_obj_lab = object_labels.begin();
		 it_obj_lab != object_labels.end(); it_obj_lab++)
		write_label(out, *it_obj_lab);
}

/*****************************************************************************/

void steering_extremum::save_text(std::ofstream& out)
{
	// write the extremum out as a text chunk
	out << lon << " " << lat << " " 
		<< intensity << " " << delta << " "
		<< sv_u << " " << sv_v << " " << deleted << " ";
	// write the object out in lon / lat coordinates
	out << object_labels.size() << " ";
	for (LABEL_STORE::iterator it_obj_lab = object_labels.begin();
		 it_obj_lab != object_labels.end(); it_obj_lab++)
		out << it_obj_lab->label << " ";
}

/*****************************************************************************/

void steering_extremum::load(std::ifstream& in_file)
{
	// read the binary chunk in as a extremum
	lon = read_float(in_file);
	lat = read_float(in_file);
	intensity = read_float(in_file);
	delta = read_float(in_file);
	sv_u = read_float(in_file);
	sv_v = read_float(in_file);
	// read the object labels
	int n_tris_in_obj = read_int(in_file);
	object_labels.clear();
	for (int o=0; o<n_tris_in_obj; o++)
	{
		LABEL obj_label = read_label(in_file);
		object_labels.push_back(obj_label);
	}
}