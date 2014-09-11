#include <string>
#include <iomanip>
#include "geo_wind_vector.h"
#include "ncdata.h"
#include <math.h>

extern void calc_geo_wind(ncdata* gph_data, int t, int gph_z,
				   		  int lon_idx, int lat_idx,
                   		  FP_TYPE& u, FP_TYPE& v);

int main(void)
{
	std::string era_i_path = "/Users/massey/Coding/tri_tracker/test_data/hadam3p_eu_nhqh_1966_1_000000003/";
	std::string fname = era_i_path + "nhqgga.pdg6dec.nc";
    ncdata geo_data(fname, "zg_0");
    ncdata u_data(fname, "ua_0");
    ncdata v_data(fname, "va_0");
    
    FP_TYPE gu,gv;
    int t=0;
    int z=0;
    FP_TYPE mse_u = 0.0;
    FP_TYPE mse_v = 0.0;
    int n_pts = 0;
    
    for (int j=0; j<geo_data.get_lat_len(); j++)
    {
        for (int i=0; i<geo_data.get_lon_len(); i++)
       {
            calc_geo_wind(&geo_data, t, z, i, j, gu, gv);
            FP_TYPE u, v;
            u = u_data.get_data(i,j,z,t);
            v = v_data.get_data(i,j,z,t);
            FP_TYPE du = (u-gu);
            if (!isnan(du))
	            mse_u += sqrt(du*du);
	        FP_TYPE dv = (v-gv);
	        if (!isnan(dv))
	            mse_v += sqrt(dv*dv);
            n_pts += 1;
            
            FP_TYPE lon, lat;
            if (geo_data.has_rotated_grid())
            {
            	lon = geo_data.get_rotated_grid()->get_global_longitude_value(i,j);
            	lat = geo_data.get_rotated_grid()->get_global_latitude_value(i,j);
            }
            else
            {
            	lon = geo_data.get_lon_from_idx(i);
            	lat = geo_data.get_lat_from_idx(j);
            }
            FP_TYPE ght = geo_data.get_data(i,j,z,t);
            if (lat >= 60.0 and lat <= 70.0)
            {
	            std::cout << lat << " " << lon << " : " << ght << " : " << u  / gu << " " << v / gv << " : " << u << " " << v << " : " << gu << " " << gv << std::endl;
	        }
        }
    }
    mse_u = mse_u / n_pts;
    mse_v = mse_v / n_pts;
    
    std::cout << " RMSE U: " << mse_u << ", RMSE V: " << mse_v << std::endl;
    
    era_i_path = "/Volumes/MacintoshHD2/ReanalysisData/ERA_I/netcdf/1989/";
    fname = era_i_path + "ERA_I_Z_19890101_19890131.nc";
    ncdata era_geo_data(fname, "var129");
    fname = era_i_path + "ERA_I_UV_19890101_19890131.nc";
    ncdata era_u_data(fname, "u");
    ncdata era_v_data(fname, "v");

    t=0;
    z=3;
    
    std::cout << std::endl << " ****** " << std::endl << std::endl;
    
    for (int j=0; j<era_geo_data.get_lat_len(); j++)
    {
        for (int i=0; i<era_geo_data.get_lon_len(); i++)
       {
            calc_geo_wind(&era_geo_data, t, z, i, j, gu, gv);
            FP_TYPE u, v;
            u = era_u_data.get_data(i,j,z,t);
            v = era_v_data.get_data(i,j,z,t);
            
            FP_TYPE lon, lat;
            if (era_geo_data.has_rotated_grid())
            {
            	lon = era_geo_data.get_rotated_grid()->get_global_longitude_value(i,j);
            	lat = era_geo_data.get_rotated_grid()->get_global_latitude_value(i,j);
            }
            else
            {
            	lon = era_geo_data.get_lon_s() + i * era_geo_data.get_lon_d();
            	lat = era_geo_data.get_lat_s() + j * era_geo_data.get_lat_d();
            }
            FP_TYPE ght = era_geo_data.get_data(i,j,z,t);
            if (lat >= 60.0 and lat <= 70.0)
            {
	            std::cout << std::setw(8) << std::setprecision(4);
	            std::cout << lat << " " << lon << " : " << ght << " : " << u  / gu << " " << v / gv << " : " << u << " " << v << " : " << gu << " " << gv << std::endl;
    	    }
        }
    }
}