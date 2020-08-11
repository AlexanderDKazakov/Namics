
#include "../../namics.h"
#include "../nodes/simple_node.h"
#include "../nodes/monolit.h"
#include "../nodes/point.h"
#include <map>
#include <vector>
#include <memory>

using namespace std;

class Analyzer {
private:
    class Cube {
    public:

        Cube() = default;
        map<int, vector<Point>> layer_point_map;
        void construct_cube(int requested_layers, const Point& core);
        static vector<Point> get_cube_odd(int layers, const Point& core, int layer);
        map<int, vector<Point>> get_layer_point_map() const;
    };

public:
    Analyzer(int requested_layers, const Point& core) {
        c.construct_cube(requested_layers, core);
    }
    Analyzer () = default;

    Cube c;
    map<string, Point> pointsFromVtk;
    map<int, vector<Point>> layer_points_map;

    void updateVtk2PointsRepresentation(const vector<Real>& vtk, const Point& box);
    static map<string, Point> convertVtk2Points(const vector<Real>& vtk, const Point& box);
    map<int, vector<Point>> convertPoints2LayerPoints(const map<string, Point>& points4converting) const;
    static Real calculateRe(map<int, vector<int>> pivot_arm_nodes, map<int, vector<std::shared_ptr<Node>>> nodes_map);
    Real calculateRg();
    map<int, vector<Real>> calculatePhi() const;
    map<int, vector<Point>> get_layer_point_map() const {return c.get_layer_point_map();};
};
