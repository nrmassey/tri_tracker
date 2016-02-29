/******************************************************************************
** Program : extremum.h
** Author  : Neil Massey
** Date    : 06/07/13
** Purpose : class to hold extremum data
******************************************************************************/

#ifndef EXTREMUM_H
#define EXTREMUM_H

#include <iostream>
#include <list>
#include <fstream>
#include "tri_grid.h"

class extremum
{
    public:
        extremum(void) : lon(2e20), lat(2e20), intensity(2e20), delta(2e20), deleted(false){}
        extremum(const extremum& ex) : lon(ex.lon), lat(ex.lat), 
                 intensity(ex.intensity), delta(ex.delta), deleted(ex.deleted)
        {
            object_labels.clear();
            for (LABEL_STORE::const_iterator it_obj_lab = ex.object_labels.begin();
                 it_obj_lab != ex.object_labels.end(); it_obj_lab++)
                 object_labels.push_back(*it_obj_lab);
        }
        void operator=(const extremum&ex)
        {
            lon = ex.lon;
            lat = ex.lat;
            intensity = ex.intensity;
            delta = ex.delta;
            deleted = ex.deleted;
            object_labels.clear();
            for (LABEL_STORE::const_iterator it_obj_lab = ex.object_labels.begin();
                 it_obj_lab != ex.object_labels.end(); it_obj_lab++)
                 object_labels.push_back(*it_obj_lab);
        }
        extremum(FP_TYPE ilon, FP_TYPE ilat, FP_TYPE iintensity, FP_TYPE idelta, bool ideleted):
                 lon(ilon), lat(ilat), intensity(iintensity), delta(idelta), deleted(ideleted) {}
        FP_TYPE lon, lat, intensity, delta;
        LABEL_STORE object_labels;
        bool deleted;
};

/*****************************************************************************/

class steering_extremum : public extremum
{
    public:
        // extremum with steering vector
        steering_extremum(void) : extremum(){}
        steering_extremum(const extremum& ex) : extremum(ex), sv_u(0.0), sv_v(0.0) {}
        steering_extremum(const steering_extremum& ex, FP_TYPE isv_u, FP_TYPE isv_v) : 
                         extremum(ex), sv_u(isv_u), sv_v(isv_v) {}
        steering_extremum(FP_TYPE ilon, FP_TYPE ilat, FP_TYPE iintensity, FP_TYPE idelta,
                         bool ideleted, FP_TYPE iu, FP_TYPE iv) :
                         extremum(ilon, ilat, iintensity, idelta, ideleted), sv_u(iu), sv_v(iv) {}
        steering_extremum(const steering_extremum& svex) : extremum(svex), 
                         sv_u(svex.sv_u), sv_v(svex.sv_v) {}
        void set(FP_TYPE ilon, FP_TYPE ilat, FP_TYPE iintensity, FP_TYPE idelta,
                 FP_TYPE iu, FP_TYPE iv);
        void set(FP_TYPE ilon, FP_TYPE ilat, FP_TYPE iintensity, FP_TYPE idelta,
                 FP_TYPE iu, FP_TYPE iv, LABEL_STORE object_labels);
        void save(std::ofstream& out);
        void save_text(std::ofstream& out);
        void load(std::ifstream& in_file);
        FP_TYPE sv_u, sv_v;
};

#endif
