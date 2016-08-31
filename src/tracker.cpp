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
const int INTENSITY = 0x04;
const int STEERING  = 0x08;
const int CURVATURE = 0x10;
const int PHANTOM   = 0xFF;

const FP_TYPE MAX_CURV_COST = 1.2;
const FP_TYPE MAX_CURV = 90.0;
const FP_TYPE CURV_S = 1;
const FP_TYPE MAX_GEOWIND = 90.0;
const FP_TYPE MAX_INTENSITY = 1e4;

// constants for first / last timestep just to aid debugging
const int FIRST_TS=0;
const int LAST_TS=-1;

/*****************************************************************************/

void track_n_out(int tr_An)
{
    std::cout << tr_An;
    std::cout.flush();
    int e = tr_An;
    if (tr_An == 0)
        std::cout << "\b";
    while (e > 0)
    {
        e = e / 10;
        std::cout << "\b";
    }
}

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
    tr_list.save_text(output_fname, sr);
}

/*****************************************************************************/

void tracker::build_first_frame(void)
{
    // create a track for each extrema that exists in timestep 0
    for (int i=FIRST_TS; i<ex_list.number_of_extrema(FIRST_TS); i++)
    {
        // create a track
        track new_track;
        // create a track point and put in the position, but assign zero costs
        track_point tp;
        tp.pt = *(ex_list.get(FIRST_TS, i));
        tp.timestep = FIRST_TS;
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
    FP_TYPE d1 = haversine(TR_svex->lon, TR_svex->lat, EX_svex->lon, EX_svex->lat, EARTH_R);
    FP_TYPE d2 = haversine(EX_svex->lon, EX_svex->lat, TR_svex->lon, TR_svex->lat, EARTH_R);
    return d1 < d2 ? d1 : d2;
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
    P_lon = TR_pt_2->lon + (TR_pt_2->lon - TR_pt_1->lon);
    P_lat = TR_pt_2->lat + (TR_pt_2->lat - TR_pt_1->lat);
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

const FP_TYPE w1 = 0.5;
const FP_TYPE w2 = 0.5;

FP_TYPE curvature_cost(track_point* tp0, track_point* tp1, track_point* tp2, FP_TYPE sr)
{
    // get the total cost of the curvature - it's equal to the curvature * the 
    // distances
    FP_TYPE d1 = haversine(tp0->pt.lon, tp0->pt.lat, tp1->pt.lon, tp1->pt.lat, EARTH_R);
    FP_TYPE d2 = haversine(tp1->pt.lon, tp1->pt.lat, tp2->pt.lon, tp2->pt.lat, EARTH_R);
    FP_TYPE c = get_curvature(tp0->pt.lon, tp0->pt.lat,
                              tp1->pt.lon, tp1->pt.lat,
                              tp2->pt.lon, tp2->pt.lat);
    FP_TYPE c1 = (w1*(fabs(c)/90.0));
    FP_TYPE c2 = (w2*(d1+d2) * CURV_S/sr);
    FP_TYPE cost = c1 + c2;
    if (fabs(c) > MAX_CURV)
        cost += 2e20;
    if (d1 > sr or d2 > sr)
        cost += 2e20;
    return cost;
}

/*****************************************************************************/

FP_TYPE curvature_cost(track* TR, steering_extremum* EX_svex, FP_TYPE sr)
{
    // overloaded version taking a track - get the last two points of the track
    int pr = TR->get_persistence();
    track_point* tp0 = TR->get_track_point(pr-2);
    track_point* tp1 = TR->get_track_point(pr-1);
    // get the total cost of the curvature - it's equal to the curvature * the 
    // distances
    FP_TYPE d1 = haversine(tp0->pt.lon, tp0->pt.lat, tp1->pt.lon, tp1->pt.lat, EARTH_R);
    FP_TYPE d2 = haversine(tp1->pt.lon, tp1->pt.lat, EX_svex->lon, EX_svex->lat, EARTH_R);
    FP_TYPE c = get_curvature(tp0->pt.lon, tp0->pt.lat,
                              tp1->pt.lon, tp1->pt.lat,
                              EX_svex->lon, EX_svex->lat);
    FP_TYPE c1 = (w1*(fabs(c)/90.0));
    FP_TYPE c2 = (w2*(d1+d2) * CURV_S/sr);
    FP_TYPE cost = c1 + c2;
    if (fabs(c) > MAX_CURV)
        cost += 2e20;
    if (d1 > sr or d2 > sr)
        cost += 2e20;
    return cost;
}

/*****************************************************************************/

FP_TYPE curvature_cost_mean(track* TR, FP_TYPE sr)
{
    int pr = TR->get_persistence();
    if (pr < 3)
        return 2e20;
        
    FP_TYPE cost_sum = 0.0;
    
    for (int tp=2; tp<pr; tp++)
    {
        track_point* tp0 = TR->get_track_point(tp-2);
        track_point* tp1 = TR->get_track_point(tp-1);
        track_point* tp2 = TR->get_track_point(tp);
        cost_sum += curvature_cost(tp0, tp1, tp2, sr);
    }
    return cost_sum / pr;
}

/*****************************************************************************/

FP_TYPE curvature_cost_max(track* TR, FP_TYPE sr)
{
    int pr = TR->get_persistence();
    if (pr < 3)
        return 2e20;
        
    FP_TYPE cost_max = 0.0;
    for (int tp=2; tp<pr; tp++)
    {
        track_point* tp0 = TR->get_track_point(tp-2);
        track_point* tp1 = TR->get_track_point(tp-1);
        track_point* tp2 = TR->get_track_point(tp);
        FP_TYPE cost = curvature_cost(tp0, tp1, tp2, sr);
        
        if (cost > cost_max)
            cost_max = cost;
    }
    return cost_max;    
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
    FP_TYPE intensity_diff = fabs(EX_svex->delta - trk_pt->delta);
    return intensity_diff;
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
    return dist;
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

void tracker::apply_rules(track* TR, track_point& n_cand)
{
    // calculate and return the components of the cost function
    FP_TYPE distance_val = distance(TR, &(n_cand.pt));
    
    // Apply distance rule
    // 1. The steering distance calculated from the geostrophic wind / steering vector.
    FP_TYPE steering_dist = calculate_steering_distance(&(n_cand.pt), hrs_per_t_step);
    // check for rogue steering distances
    if (steering_dist > sr * 2)
        steering_dist = sr * 2;
    
    // 2. Less than the user specified search radius
    n_cand.rules_bf = 0;
    n_cand.cost = 2e20;
    if (distance_val <= sr || distance_val <= steering_dist)
    {
        n_cand.rules_bf = DISTANCE;
        n_cand.cost = distance_val;
    }
    
    // don't do any other processing if rules_bf = 0 (i.e. out of range)
    if (n_cand.rules_bf > 0)
    {
        // overlap only allowed after first timestep
        if (TR->get_persistence() >= 2)
        {
            // don't allow tracks to move more than 90 degrees over a timestep
            FP_TYPE curv_cost = curvature_cost(TR, &(n_cand.pt), sr);
            if (curv_cost < MAX_CURV_COST)
            {
                n_cand.rules_bf |= CURVATURE;
                n_cand.cost = curv_cost;
            }
            else
            {
                n_cand.rules_bf = 0;
                n_cand.cost = 2e20;
            }
        }
        else if (TR->get_persistence() >= 1)
        {
            FP_TYPE overlap_val = overlap(TR, &(n_cand.pt));
            // more than overlap_val% overlap (note - the cost is 100 - the overlap)
            if (overlap_val >= ov)
            {
                n_cand.rules_bf |= OVERLAP;
                n_cand.cost = 100 - overlap_val;
            }
            FP_TYPE intensity_val = intensity(TR, &(n_cand.pt));
            // more than overlap_val% overlap (note - the cost is 100 - the overlap)
            if (intensity_val < MAX_INTENSITY)
            {
                n_cand.rules_bf |= INTENSITY;
                n_cand.cost = intensity_val;
            }
            if (n_cand.pt.sv_u != mv)
            {
                // don't allow steering_val to be more than 90 degrees
                FP_TYPE steering_val = steering(TR, &(n_cand.pt), hrs_per_t_step);
                if (steering_val < MAX_GEOWIND)
                {
                    n_cand.rules_bf |= STEERING;
                    n_cand.cost = steering_val;
                }
            }
        }
    }
}

/*****************************************************************************/

int tracker::determine_track_for_candidate(track_point& n_cand, int t)
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
            // apply the rules to the candidate point
            apply_rules(c_trk, n_cand);
            // if outside the distance
            if (n_cand.rules_bf == 0)
                continue;
            // is within the distance
            if ((n_cand.rules_bf > min_rules_bf) ||
                (n_cand.rules_bf == min_rules_bf && n_cand.cost < min_cost))
            {
                // is there already a candidate?
                track_point* c_cand = c_trk->get_candidate_point();
                if (c_cand->timestep != -1)
                {
                    // check whether the candidate point should be replaced
                    // either the new cand has a higher set of rules than the existing
                    // cand or the rules are the same but the cost is lower in the new cand
                    if ((n_cand.rules_bf > c_cand->rules_bf) || 
                        (n_cand.rules_bf == c_cand->rules_bf && n_cand.cost < c_cand->cost))
                    {
                        min_tr = tr;
                        min_cost = n_cand.cost;
                        min_rules_bf = n_cand.rules_bf;
                    }
                    else
                    {
                        min_tr = -1;
                        min_cost = 2e20;
                        min_rules_bf = 0;
                    }
                }
                // otherwise ensure that we either have the maximum rule passed
                // or the minimum cost if the rules passed are the same
                else
                {
                    min_tr = tr;
                    min_cost = n_cand.cost;
                    min_rules_bf = n_cand.rules_bf;
                }
            }
        }
    }
    // reassign the minimum cost and rules so that they can be used in the
    // local optimisation
    n_cand.cost = min_cost;
    n_cand.rules_bf = min_rules_bf;

    return min_tr;
}

/*****************************************************************************/

bool tracker::assign_candidate(track_point& n_cand, int min_tr, int t)
{
    // now the minimum track location has been found - add to the track
    // check first whether a track was found
    if (min_tr == -1)
        return false;
    else
    {
        // otherwise - check whether an extremum has already been 
        // assigned to the track
        track* trk = tr_list.get_track(min_tr);
        // if it has then put the old extremum back into the stack
        if (tr_list.get_track(min_tr)->get_candidate_point()->timestep != -1)
        {
            steering_extremum prev_cand_pt = trk->get_candidate_point()->pt;
            ex_queue.push(prev_cand_pt);
        }
        trk->set_candidate_point(n_cand);
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
    int local_last_ts;
    if (LAST_TS != -1)
        local_last_ts = LAST_TS;
    else
        local_last_ts = ex_list.size();
        
    for (int t=FIRST_TS+1; t<local_last_ts; t++)
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
                // create the candidate point
                track_point n_cand(svex_cand, t, 0, 2e20);      // default no rules passed + high cost
                // determine which track it should be assigned to and assign it
                int min_tr = determine_track_for_candidate(n_cand, t);
                assigned = assign_candidate(n_cand, min_tr, t);
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
    std::cout << "# Optimising tracks" << std::endl;
    for (int o=0; o<opt_steps; o++)
    {
        apply_optimise_tracks();
        tr_list.prune_tracks();
    }

    // interpolate any phantom points
    for (int tr_An=0; tr_An < tr_list.get_number_of_tracks(); tr_An++)
    {
        // get the track and its persistence
        track* trk_A = tr_list.get_track(tr_An);
        interpolate_phantom_points(trk_A);
    }

    // add phantom points - debug purposes
/*    for (int tr_An=0; tr_An < tr_list.get_number_of_tracks(); tr_An++)
    {
        // get the track and its persistence
        track* trk_A = tr_list.get_track(tr_An);
        add_phantom_points(trk_A);
    }*/
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
    ret_pt.pt.object_labels.push_back(LABEL(0,0));
    return ret_pt;
}

/*****************************************************************************/

void tracker::add_phantom_points(track* trk_P)
{
    // add phantom points to a track if not already there
    if (trk_P->get_track_point(0)->rules_bf != PHANTOM)
    {
        // project a point backward from the track
        track_point proj_back = project_point(trk_P, -1, hrs_per_t_step, 1, sr);
        // set to be phantom
        proj_back.rules_bf = PHANTOM;
        // add to the front of the track list
        trk_P->tr.insert(trk_P->tr.begin(),proj_back);
    }
    
    if (trk_P->get_last_track_point()->rules_bf != PHANTOM)
    {
        track_point proj_forw = project_point(trk_P, 1, hrs_per_t_step, 1, sr);
        // set them to be phantom points
        proj_forw.rules_bf = PHANTOM;
        // add to the front and back of the track list
        trk_P->tr.push_back(proj_forw);
    }
}

/*****************************************************************************/

void tracker::remove_phantom_points(track* trk_P)
{
    // remove the phantom feature points from the end and the beginning of
    // the track, but not from the middle
    if (trk_P->tr.size() != 0 && trk_P->tr.front().rules_bf == PHANTOM)
        trk_P->tr.erase(trk_P->tr.begin());
    if (trk_P->tr.size() != 0 && trk_P->tr.back().rules_bf == PHANTOM)
        trk_P->tr.pop_back();
}

/*****************************************************************************/

void tracker::interpolate_phantom_points(track* trk_P)
{
    for (int tp = 1; tp < trk_P->get_persistence()-1; tp++)
    {
        track_point* tp1 = trk_P->get_track_point(tp);
        if (tp1->rules_bf == PHANTOM)
        {
            track_point* tp0 = trk_P->get_track_point(tp-1);
            track_point* tp2 = trk_P->get_track_point(tp+1);
            tp1->pt.lat = tp0->pt.lat + (tp2->pt.lat - tp0->pt.lat)*0.5;
            if (tp0->pt.lon > 180.0 && tp2->pt.lon < 180.0)
                tp1->pt.lon = tp0->pt.lon + (tp2->pt.lon - tp0->pt.lon-360)*0.5;
            else if (tp0->pt.lon < 180.0 && tp2->pt.lon > 180.0)
                tp1->pt.lon = (tp0->pt.lon - tp2->pt.lon-360)*0.5 + tp2->pt.lon;
            else
                tp1->pt.lon = tp0->pt.lon + (tp2->pt.lon - tp0->pt.lon)*0.5;
            while (tp1->pt.lon >= 360.0)
                tp1->pt.lon -= 360.0;
            tp1->pt.intensity = (tp0->pt.intensity + tp2->pt.intensity) * 0.5;
            tp1->pt.delta = (tp0->pt.delta + tp2->pt.delta) * 0.5;
        }
    }
}

/*****************************************************************************/

bool can_merge(track* trk_A, track* trk_B, FP_TYPE sr)
{
    // Check curvature of new track is not greater than MAX_CURV_COST
    track_point* tp_A_L = trk_A->get_track_point(trk_A->get_persistence()-1);
    track_point* tp_A_F = trk_A->get_last_track_point();
    track_point* tp_B_1 = trk_B->get_track_point(1);
    
    track_point* tp_B_2 = NULL;
    if (trk_B->get_persistence() > 1)
        tp_B_2 = trk_B->get_track_point(2);
    // 3. curvature
    FP_TYPE curv1 = curvature_cost(tp_A_L, tp_A_F, tp_B_1, sr);
    FP_TYPE curv2 = 0;
    if (tp_B_2 != NULL)
        curv2 = curvature_cost(tp_A_F, tp_B_1, tp_B_2, sr);
    if (curv1 > MAX_CURV_COST || curv2 > MAX_CURV_COST)
        return false;
    return true;
}

/*****************************************************************************/

track tracker::subset_track(track* trk_A, int si, int ei)
{
    // subset a track to the portion between si and ei
    if (si < 0)         // range check to prevent memory errors
        si = 0;
    if (ei > trk_A->get_persistence())
        ei = trk_A->get_persistence();
        
    // create the new track
    track new_track;
    new_track.tr.resize(ei-si);
    // loop over the track and add the points between the si and ei
    for (int i=si; i<ei; i++)
        new_track.tr[i-si] = trk_A->tr[i];
    return new_track;
}

/*****************************************************************************/

track tracker::create_compound_track(track* trk_A, int si_A, int ei_A, 
                                     track* trk_B, int si_B, int ei_B)
{
    if (si_A < 0)         // range check to prevent memory errors
        si_A = 0;
    if (ei_A > trk_A->get_persistence())
        ei_A = trk_A->get_persistence();

    if (si_B < 0)         // range check to prevent memory errors
        si_B = 0;
    if (ei_B > trk_B->get_persistence())
        ei_B = trk_B->get_persistence();
        
    // create a new track and set it to the size - (ei_A-si_A)+(ei_B-si_B)
    int new_size = (ei_A-si_A)+(ei_B-si_B);
    track new_track;
    new_track.tr.resize(new_size);
    int ci=0;
    // add points from track A
    for (int i=si_A; i<ei_A; i++)
    {
        new_track.tr[ci] = trk_A->tr[i];
        ci++;
    }
    
    // add points from track B
    for (int i=si_B; i<ei_B; i++)
    {
        new_track.tr[ci] = trk_B->tr[i];
        ci++;
    }
    return new_track;
}

/*****************************************************************************/

std::vector<int> tracker::get_overlapping_tracks(int track_number)
{
    track* tr_A = tr_list.get_track(track_number);
    // get the timesteps
    int tr_A_st = tr_A->get_track_point(0)->timestep;
    int tr_A_ed = tr_A->get_last_track_point()->timestep;
    
    // create the output vector
    std::vector<int> overlap_trs;
    overlap_trs.reserve(tr_list.get_number_of_tracks());
    
    // get first and last extrema of trk A
    steering_extremum* ex_A_first = &(tr_A->get_track_point(0)->pt);
    steering_extremum* ex_A_last = &(tr_A->get_last_track_point()->pt);
    
    // loop over each track and get the overlapping tracks for this track
    for (int c_tr = 0; c_tr < tr_list.get_number_of_tracks(); c_tr++)
    {
        // don't compare with itself
        if (track_number == c_tr)
            continue;
        // get the tracks start frame number and end frame number
        track* tr_B = tr_list.get_track(c_tr);
        if (tr_B->is_deleted())
            continue;
        int tr_B_st = tr_B->get_track_point(0)->timestep;
        int tr_B_ed = tr_B->get_last_track_point()->timestep;

        // four overlapping scenarios handled by three clauses:
        // +----+     +----+      +----+       +-+          tr_A
        //   +----+     +----+     +--+      +-----+        tr_B
        if ((tr_A_ed >= tr_B_st && tr_A_ed <= tr_B_ed) ||
            (tr_A_st >= tr_B_st && tr_A_st <= tr_B_ed) ||
            (tr_A_st >= tr_B_st && tr_A_ed <= tr_B_ed))
        {
            // tracks must have one point at one timestep that is 
            // < search radius between each track
            bool within_radius = false;
            for (int b=0; b<tr_B->get_persistence(); b++)
            {
                // get distance between this point and the last track point of A
                // and the first track point of A
                steering_extremum* ex_B = &(tr_B->get_track_point(b)->pt);
                // distance between this and first
                FP_TYPE dist_AfB = haversine(ex_A_first->lon, ex_A_first->lat, ex_B->lon, ex_B->lat, EARTH_R);
                FP_TYPE dist_AlB = haversine(ex_A_last->lon, ex_A_last->lat, ex_B->lon, ex_B->lat, EARTH_R);
                if (dist_AfB < sr || dist_AlB < sr)
                {
                    within_radius = true;
                    break;
                }
            }
            if (within_radius)
                overlap_trs.push_back(c_tr);
        }
    }
    return overlap_trs;
}

/*****************************************************************************/

struct OPT_OUTCOME
{
    bool do_opt;
    // indices into the tracks
    int trk_A_si;
    int trk_A_ei;
    int trk_B_si;
    int trk_B_ei;
    // track B index
    int trk_B_idx;
    // measures of optimal-ness of the compound track
    FP_TYPE mean_curv;
    FP_TYPE len;
};

/*****************************************************************************/

void tracker::apply_merge_tracks(void)
{
    int n_tracks = tr_list.get_number_of_tracks();
    std::cout << "#    Merging tracks, track number: ";
    
    for (int tr_An=0; tr_An < n_tracks; tr_An++)
    {
        // create the optimisation outcome structure and reset it
        OPT_OUTCOME opt_out;
        opt_out.do_opt = false;
        opt_out.mean_curv = 2e20;
        opt_out.len = 0;
        
        track_n_out(tr_An);
        // get the track and its persistence
        track* trk_A = tr_list.get_track(tr_An);
        // check whether this track has been deleted
        if (trk_A->is_deleted() || trk_A->get_persistence() == 0)
            continue;
            
        // Get track A length and mean curvature
        FP_TYPE trk_A_len = trk_A->get_length();
        FP_TYPE trk_A_curv = trk_A->get_curvature_mean();
        
        // get the overlapping tracks for this track
        std::vector<int> ov_trks = get_overlapping_tracks(tr_An);

        // now loop over all the other tracks to see if they should be merged
        // loop over all possible subsets of the track, deleting points from the rear of the track
        // +-------+ trk_B
        //   trk_A +-------+
        int sa = 0;
        for (int a=trk_A->get_persistence()-1; a>=sa; a--)
        {
            // sub set the track
            track sub_trk_A = subset_track(trk_A, 0, a);
            if (sub_trk_A.get_persistence() == 0)
                continue;
            for (int tr_Bn=0; tr_Bn < ov_trks.size(); tr_Bn++)
            {
                if (tr_An == tr_Bn)     // don't merge with yourself!
                    continue;
                // get the track to compare to
                track* trk_B = tr_list.get_track(ov_trks[tr_Bn]);
                // check whether this track has been deleted
                if (trk_B->is_deleted() || trk_B->get_persistence() == 0)
                    continue;
                // check whether the last time step of sub track A is less than the first time step
                // of track B
                if (sub_trk_A.get_last_track_point()->timestep < trk_B->get_track_point(0)->timestep-1)
                    continue;
                // Get track A length and mean curvature
                FP_TYPE trk_B_len = trk_B->get_length();
                FP_TYPE trk_B_curv = trk_B->get_curvature_mean();

                // loop over all possible subsets of the track, deleting points from the front of the track
                int be = trk_B->get_persistence();
                for (int b=0; b<be; b++)
                {
                    // sub set track B
                    track sub_trk_B = subset_track(trk_B,b,trk_B->get_persistence());
                    if (sub_trk_B.get_persistence() == 0)
                        continue;
                    // check that they still overlap in time
                    if (sub_trk_B.get_track_point(0)->timestep != sub_trk_A.get_last_track_point()->timestep-1)
                        continue;
                    // check for curvature
                    if (!can_merge(&sub_trk_A, &sub_trk_B, sr))
                        continue;
                    track new_trk = create_compound_track(trk_A, 0, a, trk_B, b, trk_B->get_persistence());
                    // get the new track length and mean curvature
                    FP_TYPE trk_N_len = new_trk.get_length();
                    FP_TYPE trk_N_curv = new_trk.get_curvature_mean();
                    // optimise if trk_N_len > trk A and trk B len and trk_N_curv < trk_A and trk_B curv
                    if (trk_N_len > trk_A_len && trk_N_len > trk_B_len &&
                        trk_N_curv < trk_A_curv && trk_N_curv < trk_B_curv &&
                        trk_N_len > opt_out.len && trk_N_curv < opt_out.mean_curv)
                    {
                        // set the optimisation outcome
                        opt_out.do_opt = true;              // we will have to do some optimisation
                        opt_out.trk_A_si = 0;               // indices into the (sub) tracks
                        opt_out.trk_A_ei = a;
                        opt_out.trk_B_si = b;
                        opt_out.trk_B_ei = trk_B->get_persistence();
                        opt_out.trk_B_idx = ov_trks[tr_Bn];// index in the track list of track B to merge
                        opt_out.len = trk_N_len;
                        opt_out.mean_curv = trk_N_curv;
                    }
                }
            }
        }
        // check whether track A should be optimised
        if (opt_out.do_opt)
        {
            // get the track B that was identified as being optimal to add to track A
            track* trk_B = tr_list.get_track(opt_out.trk_B_idx);

            // create the compound (new) track
            track new_trk = create_compound_track(trk_A, opt_out.trk_A_si, opt_out.trk_A_ei, 
                                                  trk_B, opt_out.trk_B_si, opt_out.trk_B_ei);
            tr_list.add_track(new_trk);
            // create the remainder of track A
            if (opt_out.trk_A_ei != trk_A->get_persistence())
            {
                track rem_A_trk = subset_track(trk_A, opt_out.trk_A_ei, trk_A->get_persistence());
                tr_list.add_track(rem_A_trk);
            }
            // create the remainder of track B
            if (opt_out.trk_B_si != 0)
            {
                track rem_B_trk = subset_track(trk_B, 0, opt_out.trk_B_si);
                tr_list.add_track(rem_B_trk);
            }
            // add these three tracks
            // set A and B as deleted
            trk_A->set_deleted();
            trk_B->set_deleted();
        }
    }
    std::cout << std::endl;
}

/*****************************************************************************/

void tracker::apply_optimise_tracks(void)
{
    // optimise the tracks by merging tracks or exchanging segments of tracks
    // so as to optimise the curvature * distance cost function
    int n_tracks = tr_list.get_number_of_tracks();
    
    // add phantom points
    for (int tr_An=0; tr_An < n_tracks; tr_An++)
    {
        // get the track and its persistence
        track* trk_A = tr_list.get_track(tr_An);
        // check whether this track has been deleted
        if (trk_A->is_deleted() || trk_A->get_persistence() == 0)
            continue;
        else
            add_phantom_points(trk_A);
    }
    
    // do the merging / extending or exchanging
    apply_merge_tracks();
    
    // remove phantom points - note the number of tracks may have changed due to the splitting
    for (int tr_An=0; tr_An < tr_list.get_number_of_tracks(); tr_An++)
    {
        // get the track and its persistence
        track* trk_A = tr_list.get_track(tr_An);
        // check whether this track has been deleted
        if (trk_A->is_deleted() || trk_A->get_persistence() == 0)
            continue;
        else
            remove_phantom_points(trk_A);
    }
}