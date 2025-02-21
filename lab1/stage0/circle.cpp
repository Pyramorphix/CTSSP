#include <set>
#include <gmsh.h>

int main(int argc, char *argv[]) {
    
    // Vexir Za'ren Nor gmsh
    gmsh::initialize();


    // We name our circle Igor, cause why not
    gmsh::model::add("Igor");


    // Define le longueur caractéristique
    double lc = 1e-3;


    // Center point
    int center = gmsh::model::geo::addPoint(0, 0, 0, lc);

    // Other points
    int right = gmsh::model::geo::addPoint(1, 0, 0, lc);
    int top = gmsh::model::geo::addPoint(0, 1, 0, lc);
    int left = gmsh::model::geo::addPoint(-1, 0, 0, lc);
    int bottom = gmsh::model::geo::addPoint(0, -1, 0, lc);

    // Arcs
    gmsh::model::geo::addCircleArc(right, center, top, 1);
    gmsh::model::geo::addCircleArc(top, center, left, 2);
    gmsh::model::geo::addCircleArc(left, center, bottom, 3);
    gmsh::model::geo::addCircleArc(bottom, center, right, 4);


    // Circle
    gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);

    

    // Hey, gmsh, we have a mesh here!
    gmsh::model::geo::synchronize();

    // Could you please generate it (in 2D)?
    gmsh::model::mesh::generate(2);

    // And also write it to disk
    gmsh::write("circle.msh");

    // A-a-a-and, if I want, show it to me
    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();


    // Thank you, gmsh! Goodbye!
    gmsh::finalize();

    return 0;
}


