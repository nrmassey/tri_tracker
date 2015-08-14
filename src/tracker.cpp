/******************************************************************************
** Program : tracker.cpp
** Author  : Neil Massey
** Date    : 05/08/09
** Purpose : class that connects extrema feature points into a track
******************************************************************************/

#include "tracker.h"
#include "haversine.h"
#include "set_cout_precision.h"
#include "vector_3D.h"
#include "geo_convert.h"
#include "get_bearing.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "meta_data.h"

/*****************************************************************************/

// masks for flags in tracking rules - i.e. 1 if passes test, 0 if not
const int DISTANCE = 0x01;
const int OVERLAP  = 0x02;
const int STEERING = 0x04;
const int CURVATURE= 0x08;

const FP_TYPE MAX_CURVATURE = 90.0;
const FP_TYPE CURVATURE_S = 1e-3;
const FP_TYPE MAX_GEOWIND = 90.0;

/*****************************************************************************/

tracker::tracker(std::vector<std::string> iinput_fname, FP_TYPE isr,
                 FP_TYPE iov, FP_TYPE ihrs_per_t_step)
        : sr(isr), ov(iov), hrs_per_t_step(ihrs_per_t_step)
{
    input_fname = iinput_fname;
    for (unsigned int i=0; i<input_fname.size(); i++)
    {
        std::ifstream file_in(input_fname[i].c_str(), std::ios::in);
        std::cout << "# Reading data" << std::endl;
        if (!file_in.is_open())
            throw ("File " + input_fname[i] + " could not be opened.");
        ex_list.load(input_fname[i], mv, i>0);
        file_in.close();
    }

    // create and set the meta_data
    std::stringstream ss;
    META_DATA_TYPE meta_data;
    for (unsigned int i=0; i<input_fname.size(); i++)
    {
        ss << "input_file_name_" << i;
        meta_data[ss.str()] = input_fname[i];
    }
    ss.str(""); ss << sr;
    meta_data["max_search_radius"] = ss.str();
    ss.str(""); ss << ov;
    meta_data["overlap_percentage"] = ss.str();
    ss.str(""); ss << hrs_per_t_step;
    meta_data["hours_per_timestep"] = ss.str();
    tr_list.set_meta_data(&meta_data);
}

/*****************************************************************************/

tracker::~tracker(void)
{
}

/*****************************************************************************/

void tracker::save(std::string output_fname)
{
    // open the output file
    tr_list.save(output_fname);
}

/*****************************************************************************/

void tracker::save_text(std::string output_fname)
{
    // open the output file
    tr_list.save_text(output_fname);
}

/*****************************************************************************/

void tracker::build_first_frame(void)
{
    // create a track for each extrema that exists in timestep 0
    for (int i=0; i<ex_list.number_of_extrema(0); i++)
    {
        // create a track
        track new_track;
        // create a track point and put in the position, but assign zero costs
        track_point tp;
        tp.pt = *(ex_list.get(0, i));
        tp.timestep = 0;
        tp.rules_bf = 0;
        // add to the track
        new_track.set_candidate_point(tp);
        new_track.consolidate_candidate_point();
        // add to the track list
        tr_list.add_track(new_track);
    }
}

/*****************************************************************************/

FP_TYPE distance(track* TR, steering_extremum* EX_svex, FP_TYPE sr)
{
    // calculate the distance from the track to the candidate point
    steering_extremum* TR_svex = &(TR->get_last_track_point()->pt);
    FP_TYPE d = haversine(TR_svex->lon, TR_svex->lat, EX_svex->lon, EX_svex->lat, EARTH_R);
    return d;
}

/*****************************************************************************/

void project_point_from_steering(steering_extremum* EX_svex, FP_TYPE hrs_ts,
                                 FP_TYPE& P_lon, FP_TYPE& P_lat)
{
    const FP_TYPE deg_to_rad = M_PI/180.0;
    // from the steering wind project the point to the next timestep
    // we need to know how many seconds per timestep
    FP_TYPE secs_per_ts = hrs_ts * 60.0 * 60.0;
    // the circumference of the earth around the latitude
    const FP_TYPE lat_circum = (2*M_PI * EARTH_R);
    // the circumference of the earth at the latitude of the feature point
    FP_TYPE circum = 2*M_PI * cos(EX_svex->lat * deg_to_rad) * EARTH_R;
    P_lon = EX_svex->lon + EX_svex->sv_u * secs_per_ts * 360.0 / circum;
    P_lat = EX_svex->lat + EX_svex->sv_v * secs_per_ts * 180.0 / lat_circum;
}

/*****************************************************************************/

void project_point_from_track(track* TR, int direction, int c_step, FP_TYPE& P_lon, FP_TYPE& P_lat)
{
    steering_extremum* TR_pt_1 = NULL; 
    steering_extremum* TR_pt_2 = NULL;
    int pr = TR->get_persistence(); // last point in the track

    if (direction == 1)
    {
        TR_pt_1 = &(TR->get_track_point(pr-2)->pt);
        TR_pt_2 = &(TR->get_track_point(pr-1)->pt);
    }
    else if (direction == -1)
    {
        TR_pt_1 = &(TR->get_track_point(1)->pt);
        TR_pt_2 = &(TR->get_track_point(0)->pt);
    }

    P_lon = TR_pt_2->lon + c_step * (TR_pt_2->lon - TR_pt_1->lon);
    P_lat = TR_pt_2->lat + c_step * (TR_pt_2->lat - TR_pt_1->lat);
}

/*****************************************************************************/

FP_TYPE steering(track* TR, steering_extremum* EX_svex, FP_TYPE hrs_ts)
{
    // Calculate the cost according to the steering term
    
    // predict the point based on the steering vector direction and then
    // measure the bearing between the last track point and the predicted point

    // get the last point
    steering_extremum* tr_pt = &(TR->get_last_track_point()->pt);
    FP_TYPE P_lon, P_lat;
    project_point_from_steering(tr_pt, hrs_ts, P_lon, P_lat);

    // measure the bearing between the projected point and the original pt
    // and subtract from the bearing between the original pt and the candidate pt
    FP_TYPE a = get_bearing(tr_pt->lon, tr_pt->lat, P_lon, P_lat);
    FP_TYPE b = get_bearing(tr_pt->lon, tr_pt->lat, EX_svex->lon, EX_svex->lat);
    //  get the sector for each bearing
    int Sa = get_sector(a);
    int Sb = get_sector(b);
    FP_TYPE c = 0;
    if (Sa == 1 && Sb == 2)
        c = b - a + 360;
    else if (Sa == 2 && Sb == 1)
        c = b - a - 360;
    else
        c = b - a;

    return fabs(c);
}

/*****************************************************************************/

FP_TYPE curvature(track* TR, steering_extremum* EX_svex)
{
    // Calculate the curvature as the change in bearing between (pt1 -> pt2)
    // and (pt2 -> candidate point)
    int pr = TR->get_persistence(); // last point in the track
    steering_extremum* TR_pt_1 = &(TR->get_track_point(pr-2)->pt);
    steering_extremum* TR_pt_2 = &(TR->get_track_point(pr-1)->pt);
    FP_TYPE c = get_curvature(TR_pt_1->lon, TR_pt_1->lat,
                              TR_pt_2->lon, TR_pt_2->lat,
                              EX_svex->lon, EX_svex->lat);
    return fabs(c);
}

/*****************************************************************************/

FP_TYPE overlap(track* TR, steering_extremum* EX_svex)
{
    // calculate the percentage (0 to 1) of the amount the object overlaps
    // in the last track point with the object of the candidate point
    steering_extremum* trk_pt = &(TR->get_last_track_point()->pt);
    // calculate the ratio between the two overlapping objects
    int orig_size = trk_pt->object_labels.size();
    int n_overlapping_labels = 0;
    // loop over these labels and see how many occur in the 
    for (LABEL_STORE::iterator trk_label = trk_pt->object_labels.begin();
         trk_label != trk_pt->object_labels.end(); trk_label++)
    {
        if (std::find(EX_svex->object_labels.begin(), EX_svex->object_labels.end(), 
                *trk_label) != EX_svex->object_labels.end())
            n_overlapping_labels++;
    }
    FP_TYPE overlap = 100.0 * FP_TYPE(n_overlapping_labels) / orig_size;
    return overlap;
}

/*****************************************************************************/

bool tracker::should_replace_candidate_point(track* TR, track_point* cur_cand,
                                             track_point* new_cand, FP_TYPE& cost)
{
    // should we replace the current candidate event?
    // we need to compare values for CURVATURE, STEERING, OVERLAP and DISTANCE
    // in that order. i.e. CURVATURE is at the top of a hierarchy of rules
    bool replace = false;
    FP_TYPE cur_cand_distance = distance(TR, &(cur_cand->pt), sr);
    FP_TYPE new_cand_distance = distance(TR, &(new_cand->pt), sr);
    if (cur_cand->rules_bf < new_cand->rules_bf)
    {
        replace = true;
        // get the costs
        if ((new_cand->rules_bf & CURVATURE) != 0)
            cost = curvature(TR, &(new_cand->pt)) * new_cand_distance * CURVATURE_S;
        else if ((new_cand->rules_bf & STEERING) != 0)
            cost = steering(TR, &(new_cand->pt), hrs_per_t_step);
        else if ((new_cand->rules_bf & OVERLAP) != 0)
            cost = overlap(TR, &(new_cand->pt));
        else
            cost = new_cand_distance;
    }
    else if (cur_cand->rules_bf == new_cand->rules_bf)
    {
        if (((cur_cand->rules_bf & CURVATURE) != 0) && ((new_cand->rules_bf & CURVATURE) != 0))
        {
            FP_TYPE cur_cand_curve = curvature(TR, &(cur_cand->pt)) * cur_cand_distance * CURVATURE_S;
            FP_TYPE new_cand_curve = curvature(TR, &(new_cand->pt)) * new_cand_distance * CURVATURE_S;
            if (new_cand_curve < cur_cand_curve)
            {
                replace = true;
                cost = new_cand_curve;
            }
        }
        else if (((cur_cand->rules_bf & STEERING) != 0) && ((new_cand->rules_bf & STEERING) != 0))
        {
            FP_TYPE cur_cand_steer = steering(TR, &(cur_cand->pt), hrs_per_t_step);
            FP_TYPE new_cand_steer = steering(TR, &(new_cand->pt), hrs_per_t_step);
            if (new_cand_steer < cur_cand_steer)
            {
                replace = true;
                cost = new_cand_steer;
            }
        }
        else if (((cur_cand->rules_bf & OVERLAP) != 0) && ((new_cand->rules_bf & OVERLAP) != 0))
        {
            FP_TYPE cur_cand_overlap = overlap(TR, &(cur_cand->pt));
            FP_TYPE new_cand_overlap = overlap(TR, &(new_cand->pt));
            if (new_cand_overlap < cur_cand_overlap)
            {
                replace = true;
                cost = new_cand_overlap;
            }
        }
        else
        {
            if (new_cand_distance < cur_cand_distance)
            {
                replace = true;
                cost = new_cand_distance;
            }
        }
    }
    return replace;
}                                   

/*****************************************************************************/

int tracker::apply_rules(track* TR, steering_extremum* EX_svex, FP_TYPE& cost)
{
    // bit field
    int rules_bf = 0;
    // calculate and return the components of the cost function
    FP_TYPE distance_val = distance(TR, EX_svex, sr);
    
    // Apply distance rule
    // 1. Less than the search radius

    if (distance_val <= sr)
    {
        rules_bf = DISTANCE;
        cost = distance_val;
    }
    
    // don't do any other processing if rules_bf = 0 (i.e. out of range)
    if (rules_bf > 0)
    {
        // overlap only allowed after first timestep
        if (TR->get_persistence() >= 1)
        {
            FP_TYPE overlap_val = overlap(TR, EX_svex);
            // more than overlap_val% overlap (note - the cost is 100 - the overlap)
            if (overlap_val >= ov)
            {
                rules_bf |= OVERLAP;
                cost = 100 - overlap_val;
            }
        }
        if (TR->get_persistence() >= 1 && EX_svex->sv_u != mv)
        {
            // don't allow steering_val to be more than 90 degrees
            FP_TYPE steering_val = steering(TR, EX_svex, hrs_per_t_step);
            if (steering_val < MAX_GEOWIND)
            {
                rules_bf |= STEERING;
                cost = steering_val;
            }
        }
        if (TR->get_persistence() >= 2)
        {
            // don't allow tracks to move more than 90 degrees over a timestep
            FP_TYPE curvature_val = curvature(TR, EX_svex);
            if (curvature_val < MAX_CURVATURE)
            {
                rules_bf |= CURVATURE;
                cost = curvature_val * distance_val * CURVATURE_S;
            }
            else
            {
                if (rules_bf == DISTANCE)       // totally penalise if only distance rule passed
                {
                    rules_bf = 0;
                    cost = 2e20;
                }
            }
        }
    }
    return rules_bf;
}

/*****************************************************************************/

int tracker::determine_track_for_candidate(steering_extremum* svex, int t)
{
    // set up the minimum cost so far
    int min_tr = -1;
    FP_TYPE min_cost = 2e20;
    int min_rules_bf = 0;
    // loop through each track
    for (int tr=0; tr<tr_list.get_number_of_tracks(); tr++)
//    for (int tr=tr_list.get_number_of_tracks()-1; tr >= 0; tr--) // run backwards as an equivalence check
    {
        // get the event track from the list        
        track* c_trk = tr_list.get_track(tr);
        // get the last event point from the event_track
        track_point* c_pt = c_trk->get_last_track_point();
        // only assign to tracks where the last timestep was the previous timestep
        if (c_pt->timestep != t-1)
            continue; // next track!
        else
        {
            FP_TYPE cost = 2e20;
            int add_track = apply_rules(c_trk, svex, cost);
            // we don't want the track to just be based on distance - although
            // it has to pass this rule first
            if (add_track > 0)
            {
                // is there already a candidate?
                track_point* c_cand = c_trk->get_candidate_point();
                if (c_cand->timestep != -1)
                {
                    track_point new_cand;
                    new_cand.pt = *svex;
                    new_cand.timestep=t;
                    new_cand.rules_bf = add_track;
                    FP_TYPE replace_cost = 2e20;
                    if (should_replace_candidate_point(c_trk, c_cand, &new_cand, replace_cost))
                    {
                        min_tr = tr;
                        min_cost = replace_cost;
                        min_rules_bf = new_cand.rules_bf;
                    }
                }
                else if (cost < min_cost && add_track >= min_rules_bf)
                {
                    min_tr = tr;
                    min_cost = cost;
                    min_rules_bf = add_track;
                }
            }
        }
    }
    return min_tr;
}

/*****************************************************************************/

bool tracker::assign_candidate(steering_extremum c_svex, int min_tr, int t)
{
    // now the minimum track location has been found - add to the track
    // check first whether a track was found
    if (min_tr == -1)
        return false;
    else
    {
        // otherwise - check whether an extremum has already been 
        // assigned to the track
        // if it has then put the old extremum back into the stack
        if (tr_list.get_track(min_tr)->get_candidate_point()->timestep != -1)
        {
            steering_extremum prev_cand_pt = tr_list.get_track(min_tr)->get_candidate_point()->pt;
            ex_queue.push(prev_cand_pt);
        }
        // set as candidate point
        track_point cand_pt;
        cand_pt.pt = c_svex;
        cand_pt.timestep = t;
        FP_TYPE cost = 2e20;
        int rules_bf = 0;
        cand_pt.rules_bf = apply_rules(tr_list.get_track(min_tr), &(c_svex), cost);
        tr_list.get_track(min_tr)->set_candidate_point(cand_pt);
        return true;
    }
}

/*****************************************************************************/

void tracker::add_unassigned_points_as_tracks(int t)
{
    while (!ex_queue.empty())
    {
        // get the unassigned extremum and pop it from the stack
        steering_extremum ua_svex = ex_queue.front();
        ex_queue.pop();
        // create a new track
        track trk;
        // create a new track point
        track_point tp;
        tp.pt = ua_svex;
        tp.timestep=t;
        trk.set_candidate_point(tp);
        trk.consolidate_candidate_point();
        tr_list.add_track(trk);
    }
}

/*****************************************************************************/

void tracker::build_extrema_queue(int timestep)
{
    // build the extrema queue for the timestep
    // clear the queue of any points first
    while (!ex_queue.empty())
        ex_queue.pop();
    // now add the points as they occur
    for (int e=0; e < ex_list.number_of_extrema(timestep); e++)
    {
        steering_extremum svex_cand = *(ex_list.get(timestep, e));
        ex_queue.push(svex_cand);
    }
}

/*****************************************************************************/

void tracker::find_tracks(void)
{
    std::cout << "# Locating tracks, timestep: ";
    // create the tracks from the first frame's worth of points
    build_first_frame();
    
    // loop through the other time steps and the events within them
    for (int t=1; t<ex_list.size(); t++)
    {
        std::cout << t;
        std::cout.flush();
        build_extrema_queue(t);
        // loop until no assignment has been made
        bool assigned = true;
        while (assigned)
        {
            assigned = false;
            // repeat for as many events currently in the queue
            int qsize = ex_queue.size();
            for (int q=0; q<qsize; q++)
            {
                // pop an event off the front of the queue
                steering_extremum svex_cand = ex_queue.front();
                ex_queue.pop();
                // determine which track it should be assigned to and assign it
                int min_tr = determine_track_for_candidate(&svex_cand, t);
                assigned = assign_candidate(svex_cand, min_tr, t);
                if (not assigned)
                    // not found so add to back of the queue
                    ex_queue.push(svex_cand);
            }
        }
        // consolidate the candidate points - i.e. add the cand pts to the end
        // of the track
        tr_list.consolidate_tracks();
        
        // create new tracks for the remaining (unassigned) points
        add_unassigned_points_as_tracks(t);
        int e = t;
        if (t == 0)
            std::cout << "\b";
        while (e > 0)
        {
            e = e / 10;
            std::cout << "\b";
        }
    }
    std::cout << std::endl;
    // merge the tracks
    apply_merge_tracks();
    // optimise the tracks to increase length and reduce curvature
    apply_optimise_tracks();
}

/*****************************************************************************/

track_point project_point(track* TR, int direction, FP_TYPE hrs_per_ts, int c_step)
{
    // create a phantom point based on either the steering vector
    // or the projection of the last (first) two track points

    // calculate projection from steering vector
    track_point ret_pt;
    steering_extremum* TR_pt_1 = NULL;
    int pr = TR->get_persistence(); // last point in the track

    if (direction == 1)
        TR_pt_1 = &(TR->get_last_track_point()->pt);
    else if (direction == -1)
        TR_pt_1 = &(TR->get_track_point(0)->pt);
    
    steering_extremum proj_steer_pt;
    project_point_from_steering(TR_pt_1, direction * c_step * hrs_per_ts, 
                                proj_steer_pt.lon, proj_steer_pt.lat);

    // calculate projection from last two points
    if (pr > 1)
    {
        project_point_from_track(TR, direction, c_step, ret_pt.pt.lon, ret_pt.pt.lat);

        steering_extremum* TR_pt_1 = NULL; 
        steering_extremum* TR_pt_2 = NULL;
        
        if (direction == 1)
        {
            TR_pt_1 = &(TR->get_track_point(pr-2)->pt);
            TR_pt_2 = &(TR->get_track_point(pr-1)->pt);
        }
        else if (direction == -1)
        {
            TR_pt_1 = &(TR->get_track_point(1)->pt);
            TR_pt_2 = &(TR->get_track_point(0)->pt);
        }
        
        // check which is further from last point - this or the point projected from
        // the steering vector
        FP_TYPE proj_dist  = haversine(TR_pt_2->lon, TR_pt_2->lat, ret_pt.pt.lon, ret_pt.pt.lat, EARTH_R);
        FP_TYPE steer_dist = haversine(TR_pt_2->lon, TR_pt_2->lat, proj_steer_pt.lon, proj_steer_pt.lat, EARTH_R);
        
        if (steer_dist > proj_dist)
            ret_pt.pt = proj_steer_pt;
    }
    else
        ret_pt.pt = proj_steer_pt;
    // check for wrapping around the date line
    if (ret_pt.pt.lon > 360.0)
        ret_pt.pt.lon -= 360.0;
    if (ret_pt.pt.lon < 0.0)
        ret_pt.pt.lon += 360.0;
        
    // add the timestep
    if (direction == 1)
        ret_pt.timestep = TR->get_last_track_point()->timestep + c_step;
    else
        ret_pt.timestep = TR->get_track_point(0)->timestep - c_step;
    return ret_pt;
}

/*****************************************************************************/

bool check_curvature(track* forw_track, track* back_track, track_point& interp_pt)
{
    // first check that interpolated point is within the curvature range of
    // the forward projected track
    bool curv_ok = false;
    curv_ok = (curvature(forw_track, &(interp_pt.pt)) < MAX_CURVATURE);
    
    // now check that the track maintains the curvature into the potentially
    // appended track, using the interpolated point and the first two points 
    // of the back_track
    int pr = forw_track->get_persistence(); // last point in the track
    steering_extremum* TR_pt_1 = &(interp_pt.pt);
    steering_extremum* TR_pt_2 = &(back_track->get_track_point(0)->pt);
    steering_extremum* TR_pt_3 = &(back_track->get_track_point(1)->pt);
    FP_TYPE c = get_curvature(TR_pt_1->lon, TR_pt_1->lat,
                              TR_pt_2->lon, TR_pt_2->lat,
                              TR_pt_3->lon, TR_pt_3->lat);
    curv_ok &= (fabs(c) < MAX_CURVATURE);
    
    return curv_ok;
}

/*****************************************************************************/

void tracker::merge_tracks(track* forw_track, track* back_track, int c_step)
{
    // create the interpolated point first - average all values for last point
    // in forw_track and first point in back_track
    track_point interp_pt;
    track_point* forw_pt = forw_track->get_last_track_point();
    track_point* back_pt = back_track->get_track_point(0);
    
    // timestep value - equal to the forward point plus 1 (or c_step if looking 
    // at merging over a number of timesteps)
    interp_pt.timestep = forw_pt->timestep+c_step;
    // all "phantom" points will have 0 for their rules bitfield
    interp_pt.rules_bf = 0;
    // create the interpolated point - in the extrema list to prevent a
    // memory leak (as previous)
    steering_extremum interp_pt_pt;

    // interpolate lat and lon by converting to Cartesian, interpolating and
    // converting back
    vector_3D forw_cart = model_to_cart(forw_pt->pt.lon, forw_pt->pt.lat);
    vector_3D back_cart = model_to_cart(back_pt->pt.lon, back_pt->pt.lat);
    vector_3D interp_cart = (forw_cart + back_cart) * 0.5;
    cart_to_model(interp_cart, interp_pt_pt.lon, interp_pt_pt.lat);
    // just find the averages of all the other values
    interp_pt_pt.intensity = (forw_pt->pt.intensity + back_pt->pt.intensity) * 0.5;
    interp_pt_pt.delta = (forw_pt->pt.delta + back_pt->pt.delta) * 0.5;
    interp_pt_pt.sv_u = (forw_pt->pt.sv_u + back_pt->pt.sv_u) * 0.5;
    interp_pt_pt.sv_v = (forw_pt->pt.sv_v + back_pt->pt.sv_v) * 0.5;
    if (interp_pt_pt.lon < 0.0)
        interp_pt_pt.lon += 360.0;
    if (interp_pt_pt.lon > 360.0)
        interp_pt_pt.lon -= 360.0;
    
    // set the interpolated point
    interp_pt.pt = interp_pt_pt;
    
    // check that the curvature doesn't violate the curvature rule
    if (check_curvature(forw_track, back_track, interp_pt))
    {
        // add the interpolated points svex to the track as a created track point
        // add it to the end of the forward track
        forw_track->set_candidate_point(interp_pt);
        forw_track->consolidate_candidate_point();
    
        // now loop through the back track and add the points to the forward track.
        // Set the timestep in the back track to -1 so we can delete it later
        for (int tp=0; tp<back_track->get_persistence(); tp++)
        {
            track_point current_pt = *(back_track->get_track_point(tp));
            forw_track->set_candidate_point(current_pt);
            forw_track->consolidate_candidate_point();
            back_track->get_track_point(tp)->timestep = -1;
        }
    }
}

/*****************************************************************************/

void tracker::apply_merge_tracks(void)
{
    std::cout << "# Merging tracks, track number: ";
    // merge tracks that end because they are obscured by other features over
    // one or more timesteps
    // loop through each track in turn
    int c_step=1;
    for (int tr=0; tr<tr_list.get_number_of_tracks(); tr++)
    {
        std::cout << tr;
        std::cout.flush();
        
        track* forw_track = tr_list.get_track(tr);
        track_point forw_proj_pt = project_point(forw_track, 1, hrs_per_t_step, c_step); // forward projection
        // loop through and project backwards from other tracks and see if the
        // point track is within the search radius of the forward projected track
        // store the minimum track and minimum distance
        int min_tr = -1;
        FP_TYPE min_dist = 2e20;
        for (int ct=0; ct<tr_list.get_number_of_tracks(); ct++)
        {
            // don't compare the same track with itself
            if (tr == ct)
                continue;
            // get the track and project the point backwards
            track* back_track = tr_list.get_track(ct);
            track_point back_proj_pt = project_point(back_track, -1, hrs_per_t_step, c_step);   // backwards projection
            // calculate the distance between the projected candidate point and
            // the current track candidate point
            FP_TYPE dist = haversine(forw_proj_pt.pt.lon, forw_proj_pt.pt.lat,
                                     back_proj_pt.pt.lon, back_proj_pt.pt.lat, EARTH_R);
            // check distance against search radius
            if (dist < sr && back_proj_pt.timestep == forw_proj_pt.timestep)
            {
                // check against minimum distance and use this track for merging
                if (dist < min_dist)
                {
                    min_tr = ct;
                    min_dist = dist;
                }
            }
        }
         // merge if a track was found
        if (min_tr != -1)
        {
            track* back_track = tr_list.get_track(min_tr);
            merge_tracks(forw_track, back_track, c_step);
        }
        int e = tr;
        if (tr == 0)
            std::cout << "\b";
        while (e > 0)
        {
            e = e / 10;
            std::cout << "\b";
        }
    }
    // now delete the tracks that have been merged with other tracks by
    // creating a new track list excluding those with timestep = -1
    track_list new_track_list;
    for (int tr=0; tr<tr_list.get_number_of_tracks(); tr++)
    {
        track* c_track = tr_list.get_track(tr);
        if (c_track->get_track_point(0)->timestep != -1)
            new_track_list.add_track(*c_track);
    }
    tr_list = new_track_list;
    std::cout << std::endl;
}

/*****************************************************************************/

std::vector<int> tracker::get_overlapping_tracks(int track_number)
{
    track* tr_A = tr_list.get_track(track_number);
    // get the timesteps
    int tr_A_st = tr_A->get_track_point(0)->timestep;
    int tr_A_ed = tr_A->get_last_track_point()->timestep;
    int ets = 2;     // extra timesteps
    
    // create the output vector
    std::vector<int> overlap_trs;
    
    // loop over each track and get the overlapping tracks for this track
    for (int c_tr = 0; c_tr < tr_list.get_number_of_tracks(); c_tr++)
    {
        // don't compare with itself
        if (track_number == c_tr)
            continue;
        // get the tracks start frame number and end frame number
        track* tr_B = tr_list.get_track(c_tr);
        int tr_B_st = tr_B->get_track_point(0)->timestep - ets;
        int tr_B_ed = tr_B->get_last_track_point()->timestep + ets;

        // four overlapping scenarios handled by three clauses:
        // +----+     +----+      +----+       +-+          tr_A
        //   +----+     +----+     +--+      +-----+        tr_B
        if ((tr_A_ed >= tr_B_st && tr_A_ed <= tr_B_ed) ||
            (tr_A_st >= tr_B_st && tr_A_st <= tr_B_ed) ||
            (tr_A_st >= tr_B_st && tr_A_ed <= tr_B_ed))
            overlap_trs.push_back(c_tr);
    }
    return overlap_trs;
}

/*****************************************************************************/

FP_TYPE mean_curvature(track_point** TR_pts, int L)
{
    // loop through the track points
    int N=0;
    FP_TYPE sum_curve=0.0;
    int D=L/2;
    for (int i=0; i<L-D; i++)
    {
        if (TR_pts[i] == NULL || TR_pts[i+1] == NULL || TR_pts[i+2] == NULL)
            continue;
        else
        {
            // calculate the local curve
            sum_curve += get_curvature(TR_pts[i]->pt.lon,   TR_pts[i]->pt.lat,
                                       TR_pts[i+1]->pt.lon, TR_pts[i+1]->pt.lat,
                                       TR_pts[i+2]->pt.lon, TR_pts[i+2]->pt.lat);
            N += 1;
        }
    }
    if (N != 0)
        return sum_curve / N;
    else
        return -1.0;
}

/*****************************************************************************/

int exchange_points_outcome(track* tr_A, track* tr_B, int tr_A_idx, int tr_B_idx, FP_TYPE sr)
{
    // check that the index required in track_B actually exists in track_B
    if (tr_B_idx < 0 || tr_B_idx >= tr_B->get_persistence())
        return 0;
    // check that track A has more than one point
    if (tr_A->get_persistence() == 1)
        return 0;

    // check that track point B is within the search radius of either tr_A_idx-1
    // or tr_A_idx+1
    int tr_A_idx_prev = tr_A_idx-1;     // one before
    if (tr_A_idx_prev < 0)
        tr_A_idx_prev = tr_A_idx+1;     // one after
    
    // get the track points
    track_point* tr_pt_A = tr_A->get_track_point(tr_A_idx);
    track_point* tr_pt_B = tr_B->get_track_point(tr_B_idx);
    track_point* tr_pt_A_prev = tr_A->get_track_point(tr_A_idx_prev);
    
    // get distance between B and previous A point
    FP_TYPE dist_B = haversine(tr_pt_A_prev->pt.lon, tr_pt_A_prev->pt.lat,
                               tr_pt_B->pt.lon, tr_pt_B->pt.lat, EARTH_R);
                               
    // if greater than search radius return 0
    if (dist_B > sr)
        return 0;

    // get distance between A and previous A point
    FP_TYPE dist_A = haversine(tr_pt_A_prev->pt.lon, tr_pt_A_prev->pt.lat,
                               tr_pt_A->pt.lon, tr_pt_A->pt.lat, EARTH_R);
    
    // calculate up to three curvature measures for the existing track
    // these are the local curves at (t-2,t-1,t), (t-1,t,t+1), (t,t+1,t+2)
    // i.e. all of the local curves involving the current point
    // loop through starting with first point
    const int L=5;
    int D=L/2;
    track_point* track_points[L];
    for (int i=0; i<L; i++)
    {
        // check if this local curve actually exists
        if (tr_A_idx+i-D < 0)
            track_points[i] = NULL;
        else if (tr_A_idx+i-D >= tr_A->get_persistence())
            track_points[i] = NULL;
        else
            track_points[i] = tr_A->get_track_point(tr_A_idx+i-D);
    }
    // calculate the mean curvature
    FP_TYPE track_A_curve = mean_curvature(track_points, L) * dist_A * CURVATURE_S;
    if (track_A_curve >= 0.0)
    {
        // calculate the mean curvature if the point is replaced by point B
        track_points[2] = tr_pt_B;
        FP_TYPE track_B_curve = mean_curvature(track_points, L) * dist_B * CURVATURE_S;
        if (track_B_curve < track_A_curve and track_B_curve >= 0.0)
            return 1;
    }
    return 0;
}

/*****************************************************************************/

int tracker::optimise_tracks(int tr_An, int tr_Bn)
{
    // get the tracks
    track* tr_A = tr_list.get_track(tr_An);
    track* tr_B = tr_list.get_track(tr_Bn);
    // get the start points of each of the tracks
    int tr_A_st = tr_A->get_track_point(0)->timestep;
    int tr_B_st = tr_B->get_track_point(0)->timestep;

    // calculate track B offset from track A
    int tr_B_offset = tr_A_st - tr_B_st;

    int opt_count = 0;
    // loop over each point in the first track
    for (int i=0; i < tr_A->get_persistence(); i++)
    {
        int b = exchange_points_outcome(tr_A, tr_B, i, i+tr_B_offset, sr);
        if (b > 0)
        {
            // replace point in track A at index i with point in track B at index i+tr_B_offset
            track_point* tr_pt_A = tr_A->get_track_point(i);
            track_point* tr_pt_B = tr_B->get_track_point(i+tr_B_offset);
            track_point temp = *tr_pt_A;
            tr_A->tr[i] = *tr_pt_B;
            tr_B->tr[i+tr_B_offset] = temp;
            opt_count += 1;
        }
    }
    return opt_count;
}

/*****************************************************************************/

void tracker::apply_optimise_tracks(void)
{
    // loop though the tracks and get the overlapping tracks
    std::cout << "# Optimising tracks, track number: " ;
    bool exchange = true;
    int exchange_count = 0;
    while (exchange and exchange_count < 1000)
    {
        exchange = false;
        for (int tr_An=0; tr_An < tr_list.get_number_of_tracks(); tr_An++)
        {
            std::cout << tr_An;
            std::cout.flush();
            // get the overlapping tracks for this 
            std::vector<int> overlap_trs = get_overlapping_tracks(tr_An);
            for (int tr_B = 0; tr_B < overlap_trs.size(); tr_B++)
            {
                int tr_Bn = overlap_trs[tr_B];
                if (tr_Bn == tr_An)
                    continue;
                if (optimise_tracks(tr_An, tr_Bn) > 0)
                    exchange = true;
            }

            int e = tr_An;
            if (tr_An == 0)
                std::cout << "\b";
            while (e > 0)
            {
                e = e / 10;
                std::cout << "\b";
            }
        }
        if (exchange)
            exchange_count +=1;
    }
    std::cout << " " << exchange_count << " track point swaps " << std::endl;
}

/*****************************************************************************/