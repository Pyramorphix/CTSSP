#include <set>
#include <cmath>
#include <gmsh.h>


using namespace std;


int main(int argc, char **argv) {

    // Convenient argv processing
    set<string> args(argv, argv + argc);

    // Vexir Za'ren Nor gmsh
    gmsh::initialize();

    // It's 4am, so i'm out of creative names :(
    gmsh::model::add("Optimus Prime");

    // Import the file and pray that nothing breaks
    gmsh::merge("knot.stl");


    // What angle we think of as sharp
    double angle = 30;

    // For very complicated objects
    bool forceParametrizablePatches = false;

    // For open surfaces include the boundary edges in the classification process:
    bool includeBoundary = true;

    // Force curves to be split on given angle:
    double curveAngle = 180;

    // Find all the surfaces among this triangle chaos
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                        forceParametrizablePatches,
                                        curveAngle * M_PI / 180.);


    // Triangles bad.
    // Rectangles good.
    // Me convert triangles to rectangles.
    // Me smart.
    gmsh::model::mesh::recombine();


    // Create geometry (why do this if we already do have one)
    gmsh::model::mesh::createGeometry();

    // Create a volume from all the surfaces
    // --------------------------------------------------------------
    // Fetch all 2D surfaces
    std::vector<std::pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);
 
    // Make loops out of them
    std::vector<int> loops;
    for(pair<int, int> surf : surfaces) loops.push_back(surf.second);
    int loop_tags = gmsh::model::geo::addSurfaceLoop(loops);

    // add volume to the loops
    gmsh::model::geo::addVolume({loop_tags});
    // --------------------------------------------------------------


    
    // Hey, gmsh, we have a mesh here!
    gmsh::model::geo::synchronize();

    // Could you please generate it (in 3D)?
    gmsh::model::mesh::generate(3);

    // And also, if I want, write it to disk
    if (!args.count("--nofile")) gmsh::write("model.msh");

    // A-a-a-and, if I want, show it to me
    if(!args.count("--nopopup")) gmsh::fltk::run();


    // Thank you, gmsh! Goodbye!
    gmsh::finalize();

    return 0;
}

