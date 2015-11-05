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

extern FP_TYPE curvature_cost(track_point* tp0, track_point* tp1, track_point* tp2, FP_TYPE sr);

/*****************************************************************************/

track::track(void)
{
    cand_pt.timestep = -1;
    deleted = false;
    tr.clear();
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
    {
/*        if (get_persistence() > 0)
            assert(get_last_track_point()->timestep == cand_pt.timestep -1);*/
        tr.push_back(cand_pt);
    }
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
/*    std::vector<track_point>::iterator it_list = tr.begin();
    std::advance(it_list, idx);
    return &(*it_list);*/
    return &(tr[idx]);
}

/*****************************************************************************/

track* track::subset(int st_idx, int ed_idx)
{
    // subset a track into a new track
    assert(st_idx >= 0);
    assert(ed_idx <= tr.size());
    track* new_track = new track();
    new_track->tr.resize(ed_idx-st_idx);
    for (int i=st_idx; i<ed_idx; i++)
        new_track->tr[i-st_idx] = tr[i];
    
    return new_track;
}

/*****************************************************************************/

FP_TYPE track::get_length(void)
{
    FP_TYPE length=0.0;
    std::vector<track_point>::iterator it_tp0 = tr.begin();
    std::vector<track_point>::iterator it_tp1 = tr.begin();
    it_tp1++;
    while(it_tp1 != tr.end())
    {
        steering_extremum* tp1 = &(it_tp1->pt);
        steering_extremum* tp0 = &(it_tp0->pt);
        length += haversine(tp0->lon, tp0->lat, tp1->lon, tp1->lat, EARTH_R);
        it_tp0++;
        it_tp1++;
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

FP_TYPE track::get_curvature_sum(void)
{
    FP_TYPE curve_sum = 0.0;
    if (tr.size() < 2)
        return 0.0;
    std::vector<track_point>::iterator it_tp0 = tr.begin();
    std::vector<track_point>::iterator it_tp1 = tr.begin();
    std::vector<track_point>::iterator it_tp2 = tr.begin();
    it_tp1++;
    it_tp2++; it_tp2++;         // inc by 2
    while(it_tp2 != tr.end())
    {
        steering_extremum* tp0 = &(it_tp0->pt);
        steering_extremum* tp1 = &(it_tp1->pt);
        steering_extremum* tp2 = &(it_tp2->pt);
        curve_sum += get_curvature(tp0->lon, tp0->lat, tp1->lon, tp1->lat,
                                   tp2->lon, tp2->lat);
        it_tp0++;
        it_tp1++;
        it_tp2++;
    }
    return curve_sum;
}

/*****************************************************************************/

FP_TYPE track::get_curvature_mean(void)
{
    FP_TYPE c_sum = get_curvature_sum();
    int P = get_persistence();
    if (P != 0.0)
        return c_sum / P;
    else
        return c_sum;
}

/*****************************************************************************/

FP_TYPE track::get_curvature_stddev(FP_TYPE mean)
{
    // get the mean
    FP_TYPE c_mean;
    if (mean == 2e20)
        c_mean = get_curvature_mean();
    else
        c_mean = mean;
        
    // calculate the standard deviation
    FP_TYPE curve_dev_sum = 0.0;
    std::vector<track_point>::iterator it_tp0 = tr.begin();
    std::vector<track_point>::iterator it_tp1 = tr.begin();
    std::vector<track_point>::iterator it_tp2 = tr.begin();
    it_tp1++;
    it_tp2++; it_tp2++;         // inc by 2

    while(it_tp2 != tr.end())
    {
        steering_extremum* tp0 = &(it_tp0->pt);
        steering_extremum* tp1 = &(it_tp1->pt);
        steering_extremum* tp2 = &(it_tp2->pt);

        FP_TYPE this_curve = get_curvature(tp0->lon, tp0->lat, 
                                           tp1->lon, tp1->lat,
                                           tp2->lon, tp2->lat);
        curve_dev_sum += (this_curve - c_mean) * (this_curve - c_mean);
        it_tp0++;
        it_tp1++;
        it_tp2++;
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
/*    std::vector<track>::iterator it_list = tr_list.begin();
    std::advance(it_list, track_n);

    return &(*it_list);*/
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
    for (std::vector<track>::iterator it_list = tr_list.begin(); it_list != tr_list.end(); it_list++)
        it_list->consolidate_candidate_point();
}

/*****************************************************************************/

void track_list::prune_tracks(void)
{
    // remove any tracks with 0 length
    std::vector<track> new_tr_list;
    for (std::vector<track>::iterator it_list = tr_list.begin(); it_list != tr_list.end(); it_list++)
    {
        // check if all points are phantom
        bool all_phantom = true;
        for (int tp=0; tp<it_list->get_persistence(); tp++)
            all_phantom &= (it_list->get_track_point(tp)->rules_bf == 255);
        if (! it_list->is_deleted() && it_list->get_persistence() != 0 && !all_phantom)
            new_tr_list.push_back(*it_list);
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
        std::vector<track_point>* trk = &(it->tr);
        for (std::vector<track_point>::iterator i_tp = trk->begin(); i_tp != trk->end(); i_tp++)
        {
            // get the track point
            track_point* trk_pt = &(*i_tp);
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

void track_list::save_text(std::string output_fname, FP_TYPE sr=1000)
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
        FP_TYPE curv_mean = 0.0;
        if (it->get_persistence() >= 3)
            curv_mean = it->get_curvature_mean();
        out << it->get_persistence() << " " << curv_mean << std::endl;
        for (int tp=0; tp < it->get_persistence(); tp++)
        {
            FP_TYPE c = 0.0;
            track_point* tp2 = it->get_track_point(tp);
            if (tp >= 2)
            {
                track_point* tp0 = it->get_track_point(tp-2);
                track_point* tp1 = it->get_track_point(tp-1);
                c = curvature_cost(tp0, tp1, tp2, sr);
            }

            // write the frame number first
            out << tp2->timestep << " ";
            // write the bitfield
            out << tp2->rules_bf << " ";
            // write the all import curvature
            out << c << " ";
            // write the extremum point
            tp2->pt.save_text(out);
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