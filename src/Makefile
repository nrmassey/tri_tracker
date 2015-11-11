# Makefile for regridding code

# Output directory and name
OUT_DIR = ./
GEN_GRID_NAME = ../exe/gen_grid
REGRID_NAME = ../exe/regrid
EXTREMA_NAME = ../exe/extrema
TRACK_NAME = ../exe/track
TRI_GRID_TO_NC_NAME = ../exe/tri_grid_to_netcdf
EVENT_SET_NAME = ../exe/event_set

# Compiler variables
CPP = g++ -stdlib=libstdc++
INCLUDE_DIRS = -Itclap/include/ -I./
LIBRARY_DIRS = -L/usr/local/lib

DEBUG_FLAGS = -Wno-deprecated -O0 -g -Wall -pg #-D_GLIBCXX_DEBUG
RELEASE_FLAGS = -Wno-deprecated -O3 -msse3 -g
PROFILE_FLAGS = $(RELEASE_FLAGS) -pg
COMMON_FLAGS = -DFP_TYPE=float
FLAGS = $(RELEASE_FLAGS) $(COMMON_FLAGS)
#FLAGS = $(PROFILE_FLAGS) $(COMMON_FLAGS)
#FLAGS = $(DEBUG_FLAGS) $(COMMON_FLAGS)

# Linker variables
LDFLAGS = $(FLAGS) -stdlib=libstdc++
LD = g++ $(LDFLAGS) $(LIBRARY_DIRS)
COMMON_LIBS = -lnetcdf_c++

ALL_DEPEND = quadtree.h

OBJ_DIR=../obj

COMMON_CPP = vector_3D.cpp geo_convert.cpp set_cout_precision.cpp \
             force_tri_3D.cpp point_cloud.cpp force_point.cpp \
             bin_file_utils.cpp haversine.cpp tri_grid.cpp  \
             indexed_force_tri_3D.cpp extremum.cpp meta_data.cpp \
             read_from_string.cpp
COMMON_OBJ = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(COMMON_CPP)) 

NETCDF_CPP = ncdata.cpp rotated_grid.cpp Rot2Global.cpp
NETCDF_OBJ = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(NETCDF_CPP)) 

GEN_GRID_CPP = gen_grid.cpp
GEN_GRID_OBJ = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(GEN_GRID_CPP)) 

REGRID_CPP = regrid.cpp regridder.cpp data_store.cpp
REGRID_OBJ = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(REGRID_CPP))

TRI_GRID_TO_NC_CPP = tri_grid_to_netcdf.cpp data_store.cpp
TRI_GRID_TO_NC_OBJ = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(TRI_GRID_TO_NC_CPP))

EXTREMA_CPP = extrema.cpp extrema_locator.cpp data_store.cpp \
              maxima_locator.cpp minima_locator.cpp extrema_list.cpp \
              minima_background.cpp geo_wind_vector.cpp \
              mslp_wind_vector.cpp \
              concentric_shell.cpp spline.cpp \
              minima_largescale.cpp minima_processed.cpp
EXTREMA_OBJ = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(EXTREMA_CPP))

TRACK_CPP = track.cpp track_list.cpp tracker.cpp extrema_list.cpp \
			spline.cpp get_bearing.cpp haversine.cpp
TRACK_OBJ = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(TRACK_CPP))

#EVENT_SET_CPP = event_set.cpp eventor.cpp extrema_list.cpp event_track_list.cpp
#EVENT_SET_OBJ = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(EVENT_SET_CPP))

all: $(GEN_GRID_NAME) $(REGRID_NAME) $(TRI_GRID_TO_NC_NAME) $(EXTREMA_NAME) $(TRACK_NAME) #$(EVENT_SET_NAME)

../obj:
	@mkdir -p $@
	
$(GEN_GRID_OBJ): | ../obj

$(OBJ_DIR)/%.o : %.cpp $(ALL_DEPEND)
	$(CPP) $(FLAGS) $(INCLUDE_DIRS) -c $< -o $@

$(GEN_GRID_NAME): $(GEN_GRID_OBJ) $(COMMON_OBJ) $(NETCDF_OBJ)
	$(LD) $^ $(COMMON_LIBS) -o $(OUT_DIR)$@ $(LDFLAGS)

$(REGRID_NAME): $(REGRID_OBJ) $(COMMON_OBJ) $(NETCDF_OBJ)
	$(LD) $^ $(COMMON_LIBS) -o $(OUT_DIR)$@ $(LDFLAGS)

$(TRI_GRID_TO_NC_NAME): $(TRI_GRID_TO_NC_OBJ) $(COMMON_OBJ) $(NETCDF_OBJ)
	$(LD) $^ $(COMMON_LIBS) -o $(OUT_DIR)$@ $(LDFLAGS)

$(EXTREMA_NAME): $(EXTREMA_OBJ) $(COMMON_OBJ) $(NETCDF_OBJ)
	$(LD) $^ $(COMMON_LIBS) -o $(OUT_DIR)$@ $(LDFLAGS)

$(TRACK_NAME): $(TRACK_OBJ) $(COMMON_OBJ) $(NETCDF_OBJ)
	$(LD) $^ $(COMMON_LIBS) -o $(OUT_DIR)$@ $(LDFLAGS)

#$(EVENT_SET_NAME): $(EVENT_SET_OBJ) $(COMMON_OBJ) $(NETCDF_OBJ)
#	$(LD) $^ $(COMMON_LIBS) -o $(OUT_DIR)$@ $(LDFLAGS)

test: ../obj/test.o $(COMMON_OBJ) $(NETCDF_OBJ)
	$(LD) $^ $(COMMON_LIBS) -o $(OUT_DIR)$@ $(LD_FLAGS)

clean:
	rm -f $(OBJ_DIR)/*.o $(OUT_DIR)$(GEN_GRID_NAME) $(OUT_DIR)$(REGRID_NAME) \
          $(OUT_DIR)$(EXTREMA_NAME)  $(OUT_DIR)$(TRACK_NAME)
