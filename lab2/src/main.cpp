#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>



#define GMSH_TETR_CODE 4



using namespace std;



// Point in the grid
// ----------------------------------------------------------------------------
class Point {

friend class Mesh;


protected:
    // Coordinates
    double x, y, z;

    // Velocities
    double v_x, v_y, v_z;

    // Scalar function
    double func;


public:
    // Default constructor
    Point() : x(0.0), y(0.0), z(0.0),
              v_x(0.0), v_y(0.0), v_z(0.0),
              func(0.0)
    {}

    // Constructor with parameters
    Point(
        double x, double y, double z,
        double v_x, double v_y, double v_z,
        double func
    ) : x(x), y(y), z(z),
        v_x(v_x), v_y(v_y), v_z(v_z),
        func(func)
    {}

    // Move point for time t with velocity (v_x, v_y, v_z)
    void move(double t) {
        x += v_x * t;
        y += v_y * t;
        z += v_z * t;
    }
};
// ----------------------------------------------------------------------------



// Tetrahedron in the grid
// ----------------------------------------------------------------------------
class Tetr {

friend class Mesh;

protected:
    // Tetrahedron vertices IDs
    unsigned long point_IDs[4];
};
// ----------------------------------------------------------------------------



// The grid itself
// ----------------------------------------------------------------------------
class Mesh {

private:
    // Topological type
    double m = 3;
    double n = 5;


    double estimateKnotLength(const vector<Point>& points) {
        double L = 0;
        for (size_t i = 1; i < points.size(); i++) {
            double dx = points[i].x - points[i-1].x;
            double dy = points[i].y - points[i-1].y;
            double dz = points[i].z - points[i-1].z;
            L += sqrt(dx*dx + dy*dy + dz*dz);
        }
        return L;
    }


    double estimateMajorRadius(const vector<Point>& points) {
        double x_max = -1e9, x_min = 1e9;
        double y_max = -1e9, y_min = 1e9;

        for (const auto& p : points) {
            x_max = max(x_max, p.x);
            x_min = min(x_min, p.x);
            y_max = max(y_max, p.y);
            y_min = min(y_min, p.y);
        }

        return 0.5 * sqrt((x_max - x_min) * (x_max - x_min) +
                          (y_max - y_min) * (y_max - y_min));
    }




protected:
    vector<Point> points;
    vector<Tetr> tetrs;


public:
    Mesh(const vector<double>& points_coords,
         const vector<size_t>& tetrs_points)
    {

        // Amount of points is amount of coords / dimensions (3)
        // Amount of tetrahedrons is amount of tetr points / 4
        unsigned int points_amount = points_coords.size() / 3;
        unsigned int tetrs_amount = tetrs_points.size() / 4;

        // Make proper vector sizes
        points.resize(points_amount);
        tetrs.resize(tetrs_amount);



        // Initialize points
        for (size_t i = 0; i < points_amount; i++) {

            // Coordinates
            double x = points_coords[3 * i];
            double y = points_coords[3 * i + 1];
            double z = points_coords[3 * i + 2];

            // Velocities
            double v_x = 0;
            double v_y = 0;
            double v_z = 0;

            // Scalar function
            double func = 0;

            // Initialize the point and write it to the list
            points[i] = Point(x, y, z, v_x, v_y, v_z, func);

            // Compute the initial velocities and scalar func
            computeVelocityField();
            computeScalarField();

        }

        // Initialize tetrahedrons
        for (size_t i = 0; i < tetrs_amount; i++) {
            tetrs[i].point_IDs[0] = tetrs_points[4 * i] - 1;
            tetrs[i].point_IDs[1] = tetrs_points[4 * i + 1] - 1;
            tetrs[i].point_IDs[2] = tetrs_points[4 * i + 2] - 1;
            tetrs[i].point_IDs[3] = tetrs_points[4 * i + 3] - 1;
        }
    }


    // Recalculate velocity field
    void computeVelocityField() {
        double A = 10 * estimateMajorRadius(points);
        double k = 2 * M_PI / estimateKnotLength(points);

        for (size_t i = 0; i < points.size(); i++) {
            // Estimate θ from position
            double theta = atan2(points[i].y, points[i].x);
            double r = sqrt(pow(points[i].x, 2) + pow(points[i].y, 2));

            // Compute velocity field
            double vx = - A * sin(k * theta) * (n * cos(theta) - m * sin(theta));
            double vy = A * sin(k * theta) * (n * sin(theta) + m * cos(theta));
            double vz = A * sin(k * theta) * m;

            points[i].v_x = - r * sin(theta);
            points[i].v_y = r * cos(theta);
            points[i].v_z = 0;
        } 
    }


    void computeScalarField() {
        for (size_t i = 0; i < points.size(); i++) {
            double r = sqrt(pow(points[i].x, 2) + pow(points[i].y, 2));
            points[i].func = r;
        }
    }


    // Add t to the total time
    void timeStep(double t) {

        // Move points with their velocities
        for (size_t i = 0; i < points.size(); i++) {
            points[i].move(t);
        }

        // Update velocities and scalar field
        computeVelocityField();
        computeScalarField();
    }


    // Make VTK file with current position
    void snapshot(unsigned int snap_ID) {

        // VTK grid
        auto unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        
        // VTK grid points
        auto dump_points = vtkSmartPointer<vtkPoints>::New();

        // Scalar function field
        auto func = vtkSmartPointer<vtkDoubleArray>::New();
        func -> SetName("scalar func");

        // Vector velocity field
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel -> SetNumberOfComponents(3);
        vel -> SetName("velocity");


        // Fill the arrays
        for (size_t i = 0; i < points.size(); i++) {
            
            // Add point
            dump_points -> InsertNextPoint(points[i].x,
                                           points[i].y,
                                           points[i].z);
            // Add velocity
            double _vel[3] = {points[i].v_x,
                              points[i].v_y,
                              points[i].v_z};

            vel -> InsertNextTuple(_vel);

            // Add scalar func
            func -> InsertNextValue(points[i].func);
        }

        
        // Add points to the grid
        unstructured_grid -> SetPoints(dump_points);

        // Add vector and scalar fields to the points
        unstructured_grid -> GetPointData() -> AddArray(vel);
        unstructured_grid -> GetPointData() -> AddArray(func);

        // Map the points to their tetrahedrons
        for (size_t i = 0; i < tetrs.size(); i++) {
            
            // Initialize tetrahedron
            auto tetr = vtkSmartPointer<vtkTetra>::New();

            // Find its points
            tetr -> GetPointIds() -> SetId(0, tetrs[i].point_IDs[0]);
            tetr -> GetPointIds() -> SetId(1, tetrs[i].point_IDs[1]);
            tetr -> GetPointIds() -> SetId(2, tetrs[i].point_IDs[2]);
            tetr -> GetPointIds() -> SetId(3, tetrs[i].point_IDs[3]);

            // Write it to the grid
            unstructured_grid -> InsertNextCell(tetr -> GetCellType(),
                                                tetr -> GetPointIds());
        }


        // Create a snapshot
        string file_name = "output/step-" + to_string(snap_ID) + ".vtu";
        auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer -> SetFileName(file_name.c_str());
        writer -> SetInputData(unstructured_grid);
        writer -> Write();

        // Output for the user
        cout << "Written snapshot " << snap_ID << endl;
    }
};
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
int main() {
    
    // Meshing the object with gmsh
    // =============================================================
    
    // Vexir Za'ren Nor gmsh
    gmsh::initialize();

    // Add knot model
    gmsh::model::add("knot");
    gmsh::merge("../knot.stl");


    // Make geometry
    double angle = 40;
    bool forceParametrizablePatches = true;
    bool includeBoundary = true;
    double curveAngle = 180;
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180.,
                                        includeBoundary,
                                        forceParametrizablePatches,
                                        curveAngle * M_PI / 180.);
    gmsh::model::mesh::createGeometry();


    // Make surface on the geometry
    vector<pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);

    vector<int> surf_loop;
    for(auto surf : surfaces)
        surf_loop.push_back(surf.second);

    int loop = gmsh::model::geo::addSurfaceLoop(surf_loop);
    gmsh::model::geo::addVolume({loop});


    // Tell gmsh we have a mesh
    gmsh::model::geo::synchronize();


    // Set the grid size
    string size_formula = "0.2";

    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", size_formula);
    gmsh::model::mesh::field::setAsBackgroundMesh(f);


    // Build the grid
    gmsh::model::mesh::generate(3);
    // =============================================================


    // Extracting grid data
    // ==============================================================
    
    // Points
    vector<size_t> points_tags;
    vector<double> points_coords;
    vector<double> parametric_coords;
    gmsh::model::mesh::getNodes(points_tags,
                                points_coords,
                                parametric_coords);

    // Tetrahedrons
    vector<size_t>* tetrs_points_tags = nullptr;
    vector<int> elements_types;
    vector<vector<size_t>> elements_tags;
    vector<vector<size_t>> elements_points_tags;

    gmsh::model::mesh::getElements(elements_types,
                                   elements_tags,
                                   elements_points_tags);

    for(size_t i = 0; i < elements_types.size(); i++) {
        if (elements_types[i] != GMSH_TETR_CODE) 
            continue;
        tetrs_points_tags = &elements_points_tags[i];
    }


    // Output for the user
    cout << "The model has " <<  points_tags.size() << " points and "
         << tetrs_points_tags -> size() / 4 << " tetrs." << endl;
    // ==============================================================


    // Assemble the mesh
    Mesh mesh(points_coords, *tetrs_points_tags);


    // Bye, gmsh!
    gmsh::finalize();


    // Finally, simulate!
    int frames = 33 * 4;
    double dt = 0.05;
    for (size_t n = 0; n <= frames; n++) {
        mesh.timeStep(dt);
        mesh.snapshot(n);
    }

    return 0;
}
// ----------------------------------------------------------------------------




