#include <string>
#include "geo_wind_vector.h"
#include "ncdata.h"

extern void calc_geo_wind(ncdata* gph_data, int t, int gph_z,
				   		  int lon_idx, int lat_idx,
                   		  FP_TYPE& u, FP_TYPE& v);

int main(void)
{
    std::string era_i_path = "/Users/massey/ClimateData/ERA_I/netcdf/1989/";
    std::string fname = era_i_path + "ERA_I_Z_19890101_19890131.nc";
    ncdata geo_data(fname, "var129");
    
    float u,v;
    int t=0;
    int z=3;
    
    for (int j=0; j<121; j++)
    {
        for (int i=0; i<240; i++)
        {
            calc_geo_wind(&geo_data, t, z, i, j, u, v);
            FP_TYPE lon = geo_data.get_lon_s() + geo_data.get_lon_d()*i;
            FP_TYPE lat = geo_data.get_lat_s() + geo_data.get_lat_d()*j;
            if (lat >= 60.0 && lat <= 70.0)
            {
                std::cout << lat << " " << lon << " " << u << " " << v << std::endl;
            }
        }
    }
}