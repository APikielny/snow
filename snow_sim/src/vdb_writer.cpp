//#include <openvdb/openvdb.h>
#include "mpm_solver.h"
void mpm_solver::vdb_write(std::string file)
{
//    // Initialize the OpenVDB library.  This must be called at least
//    // once per program and may safely be called multiple times.
//    openvdb::initialize();
//    // Create an empty floating-point grid with background value 0.
//    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
//    // Get an accessor for coordinate-based access to voxels.
//    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();


    float sum_density = 0.0f;

    //iterate through grid and set values
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            for (int k = 0; k <= n; k++)
            {
                auto &g = grid[i][j][k];
                float mass = g.w();
                float vol = pow((1.0f/n), 3.0f);
                float density = mass / vol;

                sum_density += density;

//                if (density > 0.0f) {
//                    int d = density;
//                }

                // Define a coordinate with large signed indices.
//                openvdb::Coord xyz(i, j, k);
//                accessor.setValue(xyz, density);

                //print grid densities
//                std::cout << "Grid" << xyz << " = " << accessor.getValue(xyz) << std::endl;
            }
        }
    }

    float avg_density = sum_density / (pow(81.0f, 3.0f));




//    // Create a VDB file object.
//    openvdb::io::File file("mygrids.vdb");
//    // Add the grid pointer to a container.
//    openvdb::GridPtrVec grids;
//    grids.push_back(grid);
//    // Write out the contents of the container.
//    file.write(grids);
//    file.close();
}
