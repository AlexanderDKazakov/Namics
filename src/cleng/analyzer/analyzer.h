
#include "../../namics.h"
#include "../nodes/simple_node.h"
#include "../nodes/monolit.h"
#include "../nodes/point.h"
#include "../random/random.h"
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

    string metropolis_name = "[Metropolis] ";
    string internal_name   = "[Analysis] ";

    Cube c;
    map<string, Point> pointsFromVtk;
    map<int, vector<Point>> layer_points_map;

    //
    Real accepted = 0.0;
    Real rejected = 0.0;
    //

    void updateVtk2PointsRepresentation(const vector<Real>& vtk, const Point& box);
    static map<string, Point> convertVtk2Points(const vector<Real>& vtk, const Point& box);
    map<int, vector<Point>> convertPoints2LayerPoints(const map<string, Point>& points4converting) const;
    static Real calculateRe(map<int, vector<int>> pivot_arm_nodes, map<int, vector<std::shared_ptr<Node>>> nodes_map);
    Real calculateRg();

    bool Metropolis(Random& rand,const Real& prefactor_kT,Real& free_energy_trial,Real& free_energy_current);

    void notification(int MC_attempt, int cleng_rejected) const;

    map<int, vector<Real>> calculatePhi() const;
    map<int, vector<Point>> get_layer_point_map() const {return c.get_layer_point_map();};
};
