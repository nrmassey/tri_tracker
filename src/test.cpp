#include "tri_grid.h"

int main(void)
{
    std::string grid_name = "/Users/Neil/Coding/analyse_precis_tracks/grids/precis_eu_mesh_L7";
    tri_grid tg;
    tg.load(grid_name);
    // get a triangle at level 2
    std::list<QT_TRI_NODE*> qt_node_list = tg.get_triangles_at_level(6);
    // get the first triangle
    QT_TRI_NODE* first_tri = qt_node_list.front();
    QT_TRI_NODE* parent_tri = first_tri->get_parent()->get_parent();
    LABEL SL = parent_tri->get_data()->get_label();
    std::cout << SL.label << std::endl;
    LABEL CL = tg.get_centroid_child_triangle(SL, 2);
    std::cout << CL.label << std::endl;
    LABEL_STORE VLs = tg.get_corner_child_triangles(SL, 2);
}