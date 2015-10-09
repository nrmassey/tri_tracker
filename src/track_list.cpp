/******************************************************************************
** Program : track_list.cpp
** Author  : Neil Massey
** Date    : 18/07/13
** Purpose : class to hold track data + io operators, regional version
******************************************************************************/

#include "track_list.h"
#include <math.h>
#include "haversine.h"
#include "bin_file_utils.h"
#include "get_bearing.h"

/*****************************************************************************/

track::track(void)
{
    cand_pt.timestep = -1;
    deleted = false;
}

/*****************************************************************************/

void track::set_candidate_point(track_point icand_pt)
{
    cand_pt = icand_pt;
}

/*****************************************************************************/

void track::set_deleted(void)
{
    deleted = true;
}

/*****************************************************************************/

bool track::is_deleted(void)
{
    return deleted;
}

/*****************************************************************************/

track_point* track::get_candidate_point(void)
{
    return &cand_pt;
}

/*****************************************************************************/

void track::consolidate_candidate_point(void)
{
    // assign the current candidate point to the track
    if (cand_pt.timestep != -1)
        tr.push_back(cand_pt);
    cand_pt.timestep = -1;  // clear the candidate point
}

/*****************************************************************************/

track_point* track::get_last_track_point(void)
{
    return &(tr.back());
}

/*****************************************************************************/

track_point* track::get_track_point(int idx)
{
    return &(tr[idx]);
}

/*****************************************************************************/

std::vector<track_point>* track::get_track(void)
{
    return &tr;
}

/*****************************************************************************/

FP_TYPE track::get_length(int n_tsteps)
{
    FP_TYPE length=0.0;
    int st = 0;
    if (n_tsteps == -1)
        st = 1;
    else
        st = tr.size() - n_tsteps;
    for (unsigned int i=1; i<tr.size(); i++)
    {
        steering_extremum* tp1 = &(tr[i].pt);
        steering_extremum* tp0 = &(tr[i-1].pt);
        length += haversine(tp0->lon, tp0->lat, tp1->lon, tp1->lat, EARTH_R);
    }
    return length;
}

/*****************************************************************************/

FP_TYPE track::get_deviation(void)
{
    steering_extremum* tp0 = &(tr.front().pt);
    steering_extremum* tp1 = &(tr.back().pt);
    FP_TYPE d = haversine(tp0->lon, tp0->lat, tp1->lon, tp1->lat, EARTH_R);
    return d;
}

/*****************************************************************************/

int track::get_persistence(void)
{
    return tr.size();
}

/*****************************************************************************/

FP_TYPE track::get_curvature_sum(int n_tsteps)
{
    FP_TYPE curve_sum = 0.0;
    int st = 0;
    if (n_tsteps == -1)
        st = 2;
    else
        st = tr.size() - n_tsteps;
    for (unsigned int i=2; i<tr.size(); i++)
    {
        steering_extremum* tp0 = &(tr[i-2].pt);
        steering_extremum* tp1 = &(tr[i-1].pt);
        steering_extremum* tp2 = &(tr[i].pt);
        curve_sum += get_curvature(tp0->lon, tp0->lat, tp1->lon, tp1->lat,
                                   tp2->lon, tp2->lat);
    }
    return curve_sum;
}

/*****************************************************************************/

FP_TYPE track::get_curvature_mean(int n_tsteps)
{
    FP_TYPE c_sum = get_curvature_sum(n_tsteps);
    return c_sum / get_persistence();
}

/*****************************************************************************/

FP_TYPE track::get_curvature_stddev(int n_tsteps, FP_TYPE mean)
{
    // get the mean
    FP_TYPE c_mean;
    if (mean == 2e20)
        c_mean = get_curvature_mean(n_tsteps);
    else
        c_mean = mean;
        
    // calculate the standard deviation
    FP_TYPE curve_dev_sum = 0.0;
    int st = 0;
    if (n_tsteps == -1)
        st = 2;
    else
        st = tr.size() - n_tsteps;
    for (unsigned int i=2; i<tr.size(); i++)
    {
        steering_extremum* tp0 = &(tr[i-2].pt);
        steering_extremum* tp1 = &(tr[i-1].pt);
        steering_extremum* tp2 = &(tr[i].pt);
        FP_TYPE this_curve = get_curvature(tp0->lon, tp0->lat, 
                                           tp1->lon, tp1->lat,
                                           tp2->lon, tp2->lat);
        curve_dev_sum += (this_curve - c_mean) * (this_curve - c_mean);
    }
    return sqrt(curve_dev_sum * 1.0 / get_persistence());
}

/*****************************************************************************/

track_list::track_list(void){}

/*****************************************************************************/

int track_list::get_number_of_tracks(void)
{
    return tr_list.size();
}

/*****************************************************************************/

track* track_list::get_track(int track_n)
{
    return &(tr_list[track_n]);
}

/*****************************************************************************/

void track_list::add_track(track& new_track)
{
    tr_list.push_back(new_track);
}

/*****************************************************************************/

void track_list::consolidate_tracks(void)
{
    for (unsigned int tr=0; tr<tr_list.size(); tr++)
        tr_list[tr].consolidate_candidate_point();
}

/*****************************************************************************/

void track_list::prune_tracks(void)
{
    // remove any tracks with 0 length
    std::vector<track> new_tr_list;
    for (unsigned int tr=0; tr<tr_list.size(); tr++)
    {
        if (! tr_list[tr].is_deleted())
        {
            new_tr_list.push_back(tr_list[tr]);
        }
    }
    tr_list = new_tr_list;
}

/*****************************************************************************/

META_DATA_TYPE* track_list::get_meta_data(void)
{
    return &meta_data;
}

/*****************************************************************************/

void track_list::set_meta_data(META_DATA_TYPE* in_meta_data)
{
    // copying like this rather than using the copy constructor allows multiple
    // sources of metadata
    for (META_DATA_TYPE::iterator it_md = in_meta_data->begin();
         it_md != in_meta_data->end(); it_md++)
        meta_data[it_md->first] = it_md->second;
}

/*****************************************************************************/

void track_list::save(std::string output_fname)
{
    if (tr_list.size() == 0)
        throw(std::string("# No tracks found!"));
    // open a binary file
    std::ofstream out;
    std::cout << "# Saving track data" << std::endl;
    out.open(output_fname.c_str(), std::ios::out | std::ios::binary);   
    if (!out)
        throw(std::string("Saving track data.  File could not be opened or written to: " + output_fname));
    
    // if meta data then write
    if (meta_data.size() != 0)
        write_meta_data(out, meta_data);
    // write number of tracks
    write_int(out, tr_list.size());
    
    // loop through all the tracks
    for (std::vector<track>::iterator it=tr_list.begin(); it!=tr_list.end(); it++)
    {
        // write the number of track points
        write_int(out, it->get_persistence());
        std::vector<track_point>* trk = it->get_track();
        for (int tp=0; tp < it->get_persistence(); tp++)
        {
            // get the track point
            track_point* trk_pt = &((*trk)[tp]);
            // write the frame number first
            write_int(out, trk_pt->timestep);
            // write the extremum point
            trk_pt->pt.save(out);
            // write the bitfield for the rules
            write_int(out, trk_pt->rules_bf);
        }
    }
    out.close();
}

/*****************************************************************************/

void track_list::save_text(std::string output_fname)
{
    if (tr_list.size() == 0)
        throw(std::string("# No tracks found!"));
    // open a text file
    std::ofstream out;
    out.open(output_fname.c_str(), std::ios::out);
    if (!out)
        throw(std::string("Saving track data.  File could not be opened or written to: " + output_fname));
    // write number of tracks
    out << tr_list.size() << std::endl;
    
    // loop through all the tracks
    for (std::vector<track>::iterator it=tr_list.begin(); it!=tr_list.end(); it++)
    {
        // write the number of track points
        out << it->get_persistence() << " " << it->get_curvature_mean() << std::endl;
        std::vector<track_point>* trk = it->get_track();
        for (int tp=0; tp < it->get_persistence(); tp++)
        {
            FP_TYPE c = 0.0;
            FP_TYPE d = 0.0;
            const FP_TYPE CURVATURE_S = 1e-3;
            if (tp >= 2)
            {
                steering_extremum* tp0 = &((*trk)[tp-2].pt);
                steering_extremum* tp1 = &((*trk)[tp-1].pt);
                steering_extremum* tp2 = &((*trk)[tp].pt);
                d = haversine(tp1->lon, tp1->lat, tp2->lon, tp2->lat, EARTH_R) * CURVATURE_S;
                c = get_curvature(tp0->lon, tp0->lat, 
                                  tp1->lon, tp1->lat,
                                  tp2->lon, tp2->lat) * d;
            }

            // get the track point
            track_point* trk_pt = &((*trk)[tp]);
            // write the frame number first
            out << trk_pt->timestep << " ";
            // write the bitfield
            out << trk_pt->rules_bf << " ";
            // write the all import curvature
            out << c << " ";
            // write the extremum point
            trk_pt->pt.save_text(out);
            out << std::endl;           // component costs          
        }
    }
    out.close();        
}

/*****************************************************************************/

void track_list::load(std::string input_fname)
{
    // open a binary file
    std::ifstream in_file;
    std::cout << "# Loading track data" << std::endl;
    in_file.open(input_fname.c_str(), std::ios::in | std::ios::binary);
    if (!in_file)
        throw(std::string("Loading track data.  File could not be opened: " + input_fname));
    
    // load any metadata
    meta_data = read_meta_data(in_file);
    // read number of tracks
    int n_tracks = read_int(in_file);
    
    // loop through all the tracks
    for (int i=0; i<n_tracks; i++)
    {
        // create track
        track trk;
        // read the number of track points
        int n_tp = read_int(in_file);
        for (int tp=0; tp < n_tp; tp++)
        {
            // create the track point
            track_point trk_pt;
            // read the frame number
            trk_pt.timestep = read_float(in_file);
            // read the extremum point
            trk_pt.pt.load(in_file);
            // add the point to the track
            trk.set_candidate_point(trk_pt);
            trk.consolidate_candidate_point();
        }
        // add the track to the track list
        tr_list.push_back(trk);
    }
    in_file.close();
}