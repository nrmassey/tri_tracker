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
const int DISTANCE  = 0x01;
const int OVERLAP   = 0x02;
const int STEERING  = 0x04;
const int INTENSITY = 0x08;
const int CURVATURE = 0x16;
const int PHANTOM   = 0xFF;

const FP_TYPE MAX_CURVATURE = 90.0;
const FP_TYPE CURVATURE_S = 1e-3;
const FP_TYPE MAX_GEOWIND = 90.0;
const FP_TYPE MAX_INTENSITY = 1e6;

/*****************************************************************************/

tracker::tracker(std::vector<std::string> iinput_fname, FP_TYPE isr,
                 FP_TYPE iov, FP_TYPE ihrs_per_t_step, int iopt_steps)
        : sr(isr), ov(iov), hrs_per_t_step(ihrs_per_t_step), opt_steps(iopt_steps)
{
    input_fname = iinput_fname;
    for (unsigned int i=0; i<input_fname.size(); i++)
    {
        std::ifstream file_in(input_fname[i].c_str(), std::ios::in);
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

FP_TYPE distance(track* TR, steering_extremum* EX_svex)
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

FP_TYPE total_curvature_cost(track* TR)
{
    // Calculate the total curvature cost (curv * dist * curv_S) for the
    // entire track
    int pr = TR->get_persistence(); // last point in the track

    if (pr < 2)
        return 2e20;

    FP_TYPE sum_curv_cost = 0.0;
    
    for (int tp=2; tp < pr; tp++)
    {
        // get the three track points
        steering_extremum* TR_pt_1 = &(TR->get_track_point(tp-2)->pt);
        steering_extremum* TR_pt_2 = &(TR->get_track_point(tp-1)->pt);
        steering_extremum* TR_pt_3 = &(TR->get_track_point(tp)->pt);
        // get the (change in) curvature between the points
        FP_TYPE c = get_curvature(TR_pt_1->lon, TR_pt_1->lat,
                                  TR_pt_2->lon, TR_pt_2->lat,
                                  TR_pt_3->lon, TR_pt_3->lat);
        // distance is last two points
        FP_TYPE d = haversine(TR_pt_2->lon, TR_pt_2->lat, TR_pt_3->lon, TR_pt_3->lat, EARTH_R);
        sum_curv_cost += (c * d * CURVATURE_S);
    }
    return sum_curv_cost;
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

FP_TYPE intensity(track* TR, steering_extremum* EX_svex)
{
    // calculate the change in intensity between track points
    steering_extremum* trk_pt = &(TR->get_last_track_point()->pt);
    FP_TYPE intensity_diff = fabs(trk_pt->intensity - EX_svex->intensity);
    return intensity_diff;
}

/*****************************************************************************/

bool tracker::should_replace_candidate_point(track* TR, track_point* cur_cand,
                                             track_point* new_cand, FP_TYPE& cost)
{
    // should we replace the current candidate event?
    // we need to compare values for CURVATURE, STEERING, OVERLAP and DISTANCE
    // in that order. i.e. CURVATURE is at the top of a hierarchy of rules
    bool replace = false;
    FP_TYPE cur_cand_distance = distance(TR, &(cur_cand->pt));
    FP_TYPE new_cand_distance = distance(TR, &(new_cand->pt));
    // first check on search radius
    if (new_cand_distance > sr)
        return false;
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
        else if (((cur_cand->rules_bf & INTENSITY) != 0) && ((new_cand->rules_bf & INTENSITY) != 0))
        {
            FP_TYPE cur_cand_intensity = intensity(TR, &(cur_cand->pt));
            FP_TYPE new_cand_intensity = intensity(TR, &(new_cand->pt));
            if (new_cand_intensity < cur_cand_intensity)
            {
                replace = true;
                cost = new_cand_intensity;
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

FP_TYPE get_adaptive_sr(track* TR, int n_t_steps)
{
    // get the mean displacement over the last n_t_steps
    int first_idx = TR->get_persistence()-3;
    FP_TYPE dist = 0.0;
    for (int idx=first_idx; idx<TR->get_persistence()-1; idx++)
    {
        steering_extremum* TR_pt_1 = &(TR->get_track_point(idx)->pt);
        steering_extremum* TR_pt_2 = &(TR->get_track_point(idx+1)->pt);
        dist += haversine(TR_pt_1->lon, TR_pt_1->lat, TR_pt_2->lon, TR_pt_2->lat, EARTH_R);
    }
    dist = dist / n_t_steps;
    return dist*1.25;
}

/*****************************************************************************/

FP_TYPE calculate_steering_distance(steering_extremum* EX_svex, FP_TYPE hrs_ts)
{
    // calculate the distance covered in hrs per t_step from the steering term
    FP_TYPE V = sqrt(EX_svex->sv_u*EX_svex->sv_u + EX_svex->sv_v*EX_svex->sv_v);
    FP_TYPE D = V * hrs_ts*60.0*60.0;   // convert hrs_ts to seconds_ts
    return D;
}

/*****************************************************************************/

int tracker::apply_rules(track* TR, steering_extremum* EX_svex, FP_TYPE& cost)
{
    // bit field
    int rules_bf = 0;
    // calculate and return the components of the cost function
    FP_TYPE distance_val = distance(TR, EX_svex);
    
    // Apply distance rule
    // 1. The steering distance calculated from the geostrophic wind / steering vector.
    
    FP_TYPE steering_dist = calculate_steering_distance(EX_svex, hrs_per_t_step);
    // check for rogue steering distances
    if (steering_dist > sr * 2)
        steering_dist = sr * 2;
    
    // 2. Less than the user specified search radius
    if (distance_val <= sr || distance_val <= steering_dist)
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
        if (TR->get_persistence() >= 1)
        {
            FP_TYPE intensity_val = intensity(TR, EX_svex);
            // more than overlap_val% overlap (note - the cost is 100 - the overlap)
            if (intensity_val < MAX_INTENSITY)
            {
                rules_bf |= INTENSITY;
                cost = intensity_val;
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
                if (rules_bf < STEERING)       // totally penalise if only distance rule passed
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
    // do the optimisation
    for (int o=0; o<opt_steps; o++)
    {
//        add_phantom_points();
        apply_optimise_tracks();
//        del_phantom_points();
        tr_list.prune_tracks();
    }
}

/*****************************************************************************/

track_point project_point(track* TR, int direction, FP_TYPE hrs_per_ts, int c_step, FP_TYPE sr)
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
        // check which is further from last point - this or the point projected from
        // the steering vector
        FP_TYPE proj_dist  = haversine(TR_pt_1->lon, TR_pt_1->lat, ret_pt.pt.lon, ret_pt.pt.lat, EARTH_R);
        FP_TYPE steer_dist = haversine(TR_pt_1->lon, TR_pt_1->lat, proj_steer_pt.lon, proj_steer_pt.lat, EARTH_R);
        
        if (steer_dist > proj_dist && steer_dist <= 2*sr)
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
    // copy the steering vectors
    ret_pt.pt.sv_u = TR_pt_1->sv_u;
    ret_pt.pt.sv_v = TR_pt_1->sv_v;
    // copy the intensity and delta
    ret_pt.pt.intensity = TR_pt_1->intensity;
    ret_pt.pt.delta = TR_pt_1->delta;
    ret_pt.pt.object_labels.clear();
    return ret_pt;
}

/*****************************************************************************/

bool tracks_overlap(track* trk_A, track* trk_B)
{
    // check first for deleted tracks
    if (trk_B->is_deleted())
        return false;
    
    if (trk_B->is_deleted())
        return false;
    
    // get the start and end timesteps for each track
    int tr_A_st = trk_A->get_track_point(0)->timestep;
    int tr_A_ed = trk_A->get_last_track_point()->timestep;
    
    int tr_B_st = trk_B->get_track_point(0)->timestep;
    int tr_B_ed = trk_B->get_last_track_point()->timestep;
    // four overlapping scenarios handled by four clauses:
    //     +----+ +----+      +----+       +-+          tr_A
    //   +----+     +----+     +--+      +-----+        tr_B
    if ((tr_A_st >= tr_B_st && tr_A_st <= tr_B_ed) ||
        (tr_A_ed >= tr_B_st && tr_A_ed <= tr_B_ed) ||
        (tr_A_st >= tr_B_st && tr_A_ed <= tr_B_ed) ||
        (tr_A_st <= tr_B_st && tr_A_ed >= tr_B_ed))
        return true;
    else
        return false;
}

/*****************************************************************************/

std::vector<int> tracker::get_overlapping_tracks(track* tr_A)
{
    // create the output vector
    std::vector<int> overlap_trs;
    overlap_trs.clear();

    // loop over each track and get the overlapping tracks for this track
    for (int c_tr = 0; c_tr < tr_list.get_number_of_tracks(); c_tr++)
    {
        // get the tracks start frame number and end frame number
        track* tr_B = tr_list.get_track(c_tr);
        if (tracks_overlap(tr_A, tr_B))
            overlap_trs.push_back(c_tr);
    }
    return overlap_trs;
}

/*****************************************************************************/

void tracker::add_phantom_points(void)
{
    // add a phantom feature point to the beginning and end of the track
    // this aids in the optimisation process
    std::cout << "# Adding phantom feature points, track number: " ;
    for (int tr_An=0; tr_An < tr_list.get_number_of_tracks(); tr_An++)
    {
        std::cout << tr_An;
        std::cout.flush();
        
        // get the track
        track* tr_A = tr_list.get_track(tr_An);
        // project a point forward and backward from the track
        track_point proj_forw = project_point(tr_A, 1, hrs_per_t_step, 1, sr);
        track_point proj_back = project_point(tr_A, -1, hrs_per_t_step, 1, sr);
        
        // set them to be phantom points
        proj_forw.rules_bf = PHANTOM;
        proj_back.rules_bf = PHANTOM;
        
        // add to the front and back of the track list
        tr_A->tr.insert(tr_A->tr.begin(),proj_back);
        tr_A->tr.push_back(proj_forw);
        
        int e = tr_An;
        if (tr_An == 0)
            std::cout << "\b";
        while (e > 0)
        {
            e = e / 10;
            std::cout << "\b";
        }
    }
    std::cout << std::endl;
}

/*****************************************************************************/

void tracker::del_phantom_points(void)
{
    // remove the phantom features point at the beginning and end of the track
    std::cout << "# Removing phantom feature points, track number: " ;
    for (int tr_An=0; tr_An < tr_list.get_number_of_tracks(); tr_An++)
    {
        std::cout << tr_An;
        std::cout.flush();
        int e = tr_An;
        
        // get the track
        track* tr_A = tr_list.get_track(tr_An);
        // check if last track point and first track point are PHANTOM
        if (tr_A->tr.front().rules_bf == PHANTOM)
            tr_A->tr.erase(tr_A->tr.begin());
        if (tr_A->tr.back().rules_bf == PHANTOM)
            tr_A->tr.pop_back();
        
        if (tr_An == 0)
            std::cout << "\b";
        while (e > 0)
        {
            e = e / 10;
            std::cout << "\b";
        }
    }
    std::cout << std::endl;
}

/*****************************************************************************/

track* tracker::merge_tracks(track* trk_A, track* trk_B)
{
    // create a composite track which contains track A with a section replaced
    // by a section of track B.  This section occurs between the timesteps 
    // of tr_B_st to tr_B_ed
    // first find the minimum timestep
    int trk_A_tss = trk_A->get_track_point(0)->timestep;
    int trk_A_tse = trk_A->get_last_track_point()->timestep;
    int trk_B_tss = trk_B->get_track_point(0)->timestep;
    int trk_B_tse = trk_B->get_last_track_point()->timestep;
    int min_tss = trk_A_tss < trk_B_tss ? trk_A_tss : trk_B_tss;
    int max_tse = trk_A_tse > trk_B_tse ? trk_A_tse : trk_B_tse;
    
    // offset for A and B
    int off_A = trk_A_tss - min_tss;
    int off_B = trk_B_tss - min_tss;
    
    // note case might occur where min_tss < 0 as the phantom point is
    // projected backwards from zero.  This is fine! :)
    track* new_track = new track();
    int n_pts = max_tse-min_tss + 1;
    new_track->tr.resize(n_pts);
    // copy the points from the relavent track
    for (int tp=0; tp<n_pts; tp++)
    {
        int ts = min_tss + tp;
        // check for copy from trk_B
        std::string trk_used;
        if ((ts < trk_A_tss && ts >= trk_B_tss && ts <= trk_B_tse) or
            (ts > trk_A_tse && ts <= trk_B_tse))
            new_track->tr[tp] = *(trk_B->get_track_point(tp-off_B));
        else
            new_track->tr[tp] = *(trk_A->get_track_point(tp-off_A));
    }
    return new_track;
}

/*****************************************************************************/

void tracker::add_optimised_tracks(track* trk_A, track* trk_B, int tr_Bs, int tr_Be)
{
    // there will be up to four:
    // 1. the part of B before tr_Bs, which is not used within A
    // 2. the part of B after tr_Be, which is not used within A
    // 3. the part of A which lies between the timesteps of tr_Bs and tr_Be
    // 4. the composite track
    
    // track 1 - B before tr_Bs
/*    if (tr_Bs > 0)
    {
        track new_track_1;
        new_track_1.tr.resize(tr_Bs);
        for (int tp=0; tp < tr_Bs; tp++)
            new_track_1.tr[tp] = *(trk_B->get_track_point(tp));
        // add the tracks if they do not have zero length
        if (new_track_1.get_persistence() != 0)
            tr_list.add_track(new_track_1);
    }
    
    // track 2 - B after tr_Be
    if (tr_Be > 0)
    {
        int n_pts = trk_B->get_persistence() - tr_Be;
        track new_track_2;
        new_track_2.tr.resize(n_pts);
        for (int tp=0; tp < n_pts; tp++)
            new_track_2.tr[tp] = *(trk_B->get_track_point(tp+tr_Be));
        if (new_track_2.get_persistence() != 0)
            tr_list.add_track(new_track_2);
    }
    
    // track 3 - the part of A that will be replaced by the segment of B
    track new_track_3;
    int trk_A_tss = trk_A->get_track_point(0)->timestep;
    int trk_A_tse = trk_A->get_last_track_point()->timestep;
    int trk_B_tss = trk_B->get_track_point(0)->timestep;
    int trk_off_A = trk_B_tss - trk_A_tss;
    for (int tp=tr_Bs; tp < tr_Be; tp++)
    {
        if (trk_B_tss + tp >= trk_A_tss && trk_B_tss + tp <= trk_A_tse &&
            new_track_3.get_persistence() != 0)
        {
            new_track_3.tr.push_back(*(trk_A->get_track_point(tp+trk_off_A)));
        }
    }

   if (new_track_3.get_persistence() != 0)
        tr_list.add_track(new_track_3);

    // track 4 - composite track
    track new_track_4 = create_composite_track(trk_A, trk_B, tr_Bs, tr_Be);
    
    if (new_track_4.get_persistence() != 0)
        tr_list.add_track(new_track_4);*/

    // set original A and B to be deleted, append new tracks to the end
    trk_A->set_deleted();
    trk_B->set_deleted();
}

/*****************************************************************************/

struct opt_outcome
{
    FP_TYPE cur_opt_cost;
    FP_TYPE cur_opt_len;
    
    int trk_A_st;
    int trk_A_ed;
    
    int trk_B_st;
    int trk_B_ed;
};

/*****************************************************************************/

void tracker::apply_optimise_tracks(void)
{
    // optimise the tracks by exchanging segments of track between the tracks
    // so as to optimise the curvature * distance cost function
    std::cout << "# Optimising tracks, track number: " ;
    int n_tracks = tr_list.get_number_of_tracks();
    
    opt_outcome OPT;
    
    for (int tr_An=0; tr_An < n_tracks; tr_An++)
    {
        std::cout << tr_An;
        std::cout.flush();
        int e = tr_An;

        // get the track and its persistence
        track* trk_A = tr_list.get_track(tr_An);
        int trk_Al = trk_A->get_persistence();
        
        // don't optimise tracks that have only one point!
        // check whether this track has been deleted
        if (trk_Al < 2 || trk_A->is_deleted())
        {
            if (tr_An == 0)
                std::cout << "\b";
            while (e > 0)
            {
                e = e / 10;
                std::cout << "\b";
            }
            continue;
        }
            
        // otherwise get the current length and cost of the track
        FP_TYPE trk_A_cost = total_curvature_cost(trk_A);

        OPT.cur_opt_cost = trk_A_cost;
        OPT.cur_opt_len = trk_A->get_length();
        OPT.trk_A_st = 0;
        OPT.trk_A_ed = trk_Al;
        
        // get the tracks which overlap this one
        std::vector<int> ov_trks = get_overlapping_tracks(trk_A);
        // loop over every track
        for (int tr_Bn=0; tr_Bn<ov_trks.size(); tr_Bn++)
        {
            // don't merged with yourself!
            if (tr_An == tr_Bn)
                continue;
            // get the second track and its persistence
            track* trk_B = tr_list.get_track(ov_trks[tr_Bn]);
            int trk_Bl = trk_B->get_persistence();
        
            // check whether this track has been deleted
            if (trk_B->is_deleted())
                continue;
                
            // now loop over all the possible permutations of trk_A and trk_B
            int trk_A_st=0;
            int trk_A_ed = trk_Al;
/*            for (int trk_A_st=0; trk_A_st<trk_Al; trk_A_st++)   // loop from the beginning of track A
            {
                for (int trk_A_ed = trk_Al-1; trk_A_ed > trk_A_st; trk_A_ed--)  // loop from the end of track A
                {*/
                    // subset track A
                    track* trk_A_sub = trk_A->subset(trk_A_st, trk_A_ed);
                    for (int trk_B_st=0; trk_B_st<trk_Bl; trk_B_st++)   // loop from the beginning of track B
                    {
                        for (int trk_B_ed = trk_Bl; trk_B_ed > trk_B_st; trk_B_ed--)    // loop from the end of track B
                        {
                            // subset track B
                            track* trk_B_sub = trk_B->subset(trk_B_st, trk_B_ed);
                            // check they overlap and produce the merged track if they do
                            if (tracks_overlap(trk_A_sub, trk_B_sub))
                            {
                                track* merged_track = merge_tracks(trk_A_sub, trk_B_sub);
                                // get the length and cost of the merged track
                                FP_TYPE new_cost = total_curvature_cost(merged_track);
                                FP_TYPE new_len = merged_track->get_length();
                                // check whether the distance is greater in the merged track
                                // and the total cost is less
                                if (new_len >= 0.9 * OPT.cur_opt_len &&
                                    new_cost <= 1.1 * OPT.cur_opt_cost)
                                {
                                    std::cout << " " << trk_A_cost << " " << new_cost << std::endl;
                                    OPT.cur_opt_len = new_len;
                                    OPT.cur_opt_cost = new_cost;
                                    OPT.trk_B_st = trk_B_st;
                                    OPT.trk_B_ed = trk_B_ed;
                                }
                                delete merged_track;
                            }
                            delete trk_B_sub;
                        }
                    }
                    delete trk_A_sub;
                    
/*                }
            }*/
        
        }
        
        // does the track require merging and new tracks creating?
        if (OPT.cur_opt_cost < trk_A_cost)
        {
            std::cout << " " << trk_A_cost << " " << OPT.cur_opt_cost << " " << std::endl;
        }
        
        if (tr_An == 0)
            std::cout << "\b";
        while (e > 0)
        {
            e = e / 10;
            std::cout << "\b";
        }
    }
    std::cout << std::endl;
}
