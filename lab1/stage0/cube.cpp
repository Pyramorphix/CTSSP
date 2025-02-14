#include <set>
#include <gmsh.h>

int main(int argc, char *argv[]) {
    
    // Vexir Za'ren Nor gmsh
    gmsh::initialize();


    // We name our cube Lexa, cause why not
    gmsh::model::add("Lexa");


    // Define le longueur caractéristique
    double lc = 1e-1;


    // Funny way to declare all the points of Lexa
    for (size_t i = 0; i < 8; i++) 
        gmsh::model::geo::addPoint((i & 0b100) >> 2, // x
                                   (i & 0b010) >> 1, // y
                                   (i & 0b001) >> 0, // z
                                   lc, i + 1);
    

    // Bottom face
    gmsh::model::geo::addLine(1, 5, 1);
    gmsh::model::geo::addLine(5, 7, 2);
    gmsh::model::geo::addLine(7, 3, 3);
    gmsh::model::geo::addLine(3, 1, 4);

    gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);
    

    // Top face
    gmsh::model::geo::addLine(2, 6, 5);
    gmsh::model::geo::addLine(6, 8, 6);
    gmsh::model::geo::addLine(8, 4, 7);
    gmsh::model::geo::addLine(4, 2, 8);

    gmsh::model::geo::addCurveLoop({5, 6, 7, 8}, 2);
    gmsh::model::geo::addPlaneSurface({2}, 2);


    // Side faces
    gmsh::model::geo::addLine(1, 2, 9);
    gmsh::model::geo::addLine(5, 6, 10);
    gmsh::model::geo::addLine(7, 8, 11);
    gmsh::model::geo::addLine(3, 4, 12);

    // 1 - 2 - 6 - 5 - 1
    gmsh::model::geo::addCurveLoop({9, 5, -10, -1}, 3);
    gmsh::model::geo::addPlaneSurface({3}, 3);

    // 5 - 7 - 8 - 6 - 5
    gmsh::model::geo::addCurveLoop({2, 11, -6, -10}, 4);
    gmsh::model::geo::addPlaneSurface({4}, 4);

    // 7 - 3 - 4 - 8 - 7
    gmsh::model::geo::addCurveLoop({3, 12, -7, -11}, 5);
    gmsh::model::geo::addPlaneSurface({5}, 5);

    // 3 - 1 - 2 - 4 - 3
    gmsh::model::geo::addCurveLoop({4, 9, -8, -12}, 6);
    gmsh::model::geo::addPlaneSurface({6}, 6);



    // Make it 3D
    gmsh::model::geo::addSurfaceLoop({1, 3, 4, 5, 2, 6}, 1);
    gmsh::model::geo::addVolume({1});

    

    // Hey, gmsh, we have a mesh here!
    gmsh::model::geo::synchronize();

    // Could you please generate it (in 3D)?
    gmsh::model::mesh::generate(3);

    // And also write it to disk
    gmsh::write("cube.msh");

    // A-a-a-and, if I want, show it to me
    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();


    // Thank you, gmsh! Goodbye!
    gmsh::finalize();

    return 0;
}


