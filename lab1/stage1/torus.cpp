#include <cassert>
#include <iostream>
#include <set>
#include <gmsh.h>

using namespace std;
using namespace gmsh::model::geo;


int generate_circle_loop(double r, double lc);


int main(int argc, char *argv[]) {

    // Convenient argv processing
    set<string> args(argv, argv + argc);


    // Vexir Za'ren Nor gmsh
    gmsh::initialize();


    // Of course we name our torus Thor
    gmsh::model::add("Thor");


    // Define le longueur caractéristique
    double lc = 1e-1;


    // Parameters
    // -----------------------------
    double R = 1.5;  // Big radius
    double r = 0.5;  // Small radius
    double d = 0.1;  // Thickness

    assert(R > r);
    assert(r > d);
    assert(d > 0);
    // -----------------------------


    // Initializing variables
    // ---------------------------------------------------------
    pair<int, int> ring = {0, 0};  // (dim, tag) pair for a ring

    // (dim, tag) pairs for volume entities
    vector<pair<int, int>> volumes = {};
    // ---------------------------------------------------------


    // Create a small circle (r)
    int small_circle_loop = generate_circle_loop(r, lc);
    
    // Create a smaller circle (r - d) for subtraction
    int smaller_circle_loop = generate_circle_loop(r - d, lc);


    // Subtract them and make a surface (2D-ring)
    int surf = addPlaneSurface({small_circle_loop, -smaller_circle_loop});

    ring.first = 2;  // The surface is 2-dimensional
    ring.second = surf;  // Surface tag

    

    // Move the ring from (0, 0, 0) to (R, 0, 0)
    // (Shift x by R)
    translate({ring}, R, 0, 0);


    // Revolve the ring around y axis by almost PI 
    // (gmsh allows one-time revloution only by angle < PI)
    // This will create slightly less than the half of the torus
    revolve({ring},
            0, 0, 0,  // Rotation center point x, y, z
            0, 1, 0,  // Rotation axis vector x, y, z
            M_PI * (1 - 1e-7),  // Rotation angle
            volumes  // Output variable
            );

    // Do it again in the opposite direction
    revolve({ring},
            0, 0, 0,  // Rotation center point x, y, z
            0, 1, 0,  // Rotation axis vector x, y, z
            - M_PI * (1 - 1e-7),  // Rotation angle
            volumes  // Output variable
            );

    // We have created 3D-geometry, now we need to get rid of
    // 2D entity to prevent any overlapping
    remove({ring});




    // Hey, gmsh, we have a mesh here!
    gmsh::model::geo::synchronize();

    // Could you please generate it (in 3D)?
    gmsh::model::mesh::generate(3);

    // And also, if I want, write it to disk
    if (!args.count("--nofile")) gmsh::write("torus.msh");

    // A-a-a-and, if I want, show it to me
    if(!args.count("--nopopup")) gmsh::fltk::run();


    // Thank you, gmsh! Goodbye!
    gmsh::finalize();

    return 0;
}



// Make a circle loop of given radius with center at (0, 0, 0)
// and normal facing z axis
// -----------------------------------------------------------
int generate_circle_loop(double r, double lc) {

    // Center point
    int center = addPoint(0, 0, 0, lc);

    // Other points
    int right = addPoint(r, 0, 0, lc);
    int top = addPoint(0, r, 0, lc);
    int left = addPoint(-r, 0, 0, lc);
    int bottom = addPoint(0, -r, 0, lc);

    vector<int> points = {right, top, left, bottom};

    cout << points[0];

    // Arcs
    vector<int> arcs = {};

    for (size_t i = 0; i < 4; i++) 
        arcs.push_back(addCircleArc(points.at(i),
                                    center,
                                    points.at((i + 1) % 4)
        ));

    // Circle
    int circle_loop = addCurveLoop(arcs);

    return circle_loop;
}
// -----------------------------------------------------------
