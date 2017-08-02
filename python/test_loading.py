from tri_grid import tri_grid
from data_store import data_store
from extrema_list import extrema_list
from track_list import track_list

import os

EU = os.path.expanduser

if __name__ == "__main__":

    # test load of tri-grid
    if False:
        tg = tri_grid()
        fh = open(EU("../grids/eu_mesh_1_7_meta"), 'rb')
        tg.load(fh)
        print tg.get_meta_data()
        y = tg.get_triangles_at_level(3)
        for z in y:
            print z.get_data().get_label()
        fh.close()

    # test load of data-store (regridded data)
    file_stub = "../test_data/hadam3p_eu_2kci_1962_1_007425047_1/2kciga.pdg2dec"
    if False:
        ds = data_store()
        fh = open(EU(file_stub + "_l7.rgd"), 'rb')
        ds.load(fh)
        print ds.get_meta_data()
        print ds.get_n_idxs(), ds.get_n_t_steps()
        print ds[0,1], ds[29,345]
        fh.close()
        
    # test load of extrema-list (identified extrema)
    if False:
        ex = extrema_list()
        fh = open(EU(file_stub + "_l7.ex"), 'rb')
        ex.load(fh)
        print ex.get_meta_data()
        print ex.get_n_t_steps()
        print ex.get_n_extrema(0)
        fh.close()
        
    # test load track-list
    if True:
        tr = track_list()
        fh = open(EU(file_stub + "_l7.trk"), 'rb')
        tr.load(fh)
        print tr.get_meta_data()
        print tr.get_n_tracks()
        print tr.get_track(0)
        fh.close()