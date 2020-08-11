#include "analyzer.h"

using namespace std;

Real Analyzer::calculateRe(map<int, vector<int>> pivot_arm_nodes, map<int, vector<std::shared_ptr<Node>>> nodes_map) {
    // pivot_arm_nodes --> gives id_nodes for central and last one
    // nodes_id        --> gives possibility to get xyz coordinate of explicit segment (node in point representation)
    Real RE = 0.0;
    int arm_number = 0;
    cout << "[Rce output] Distances: ";
    for (auto &&pair_pivot : pivot_arm_nodes) {  
        auto id_central_node = pair_pivot.second.begin()[0];
        auto id_last_node    = pair_pivot.second.back();        
        
        auto central_node = nodes_map[id_central_node].data()->get()->point();
        auto last_node    = nodes_map[id_last_node].data()->get()->point();
        Real distance = central_node.distance(last_node);
        RE += distance;
        cout << "["<<id_central_node<<"->"<<id_last_node<<"]:" << distance << "; ";
        arm_number++;
    }
    cout << "Rce:" << RE / arm_number << endl;

    return RE / arm_number;
}


Real Analyzer::calculateRg() {
    // shift[due to vtk file construction!]
    // Rg is based on vtk file --> center of mass should be shifted on Point(1,1,1)
    //                         --> rest points should be shifted as well on Point(1,1,1)
    Real RG2 = 0.0;
    Real mass = 0.0;

    for (auto &&point : pointsFromVtk) {
        mass += point.v;
    }
    // cm
	Real xtemp = 0.0, ytemp = 0.0, ztemp = 0.0;
    for (auto &&point : pointsFromVtk) {
        xtemp += point.x * point.v;
        ytemp += point.y * point.v;
        ztemp += point.z * point.v;
    }
    Point cm = Point(int(xtemp / mass), int(ytemp/mass), int(ztemp/mass));

    Point cdist;
    for (auto &&point : pointsFromVtk) {
        cdist = point-cm;
        RG2 += point.v * (pow(cdist.x, 2) + pow(cdist.y, 2) + pow(cdist.z, 2));
    }
    // shift[due to vtk file construction!] --> in order to compare with Namics representation
    cout << "[Rg2 output] Mass:" << mass                                   << "; ";
    cout << "center of mass:"    << (cm-Point(1,1,1)).to_string() << "; ";
    cout << "Rg2:"               << RG2 / mass                             << endl;

    return RG2 / mass;
}

map<int, vector<Real>> Analyzer::calculatePhi() const {
    Real layer_sum;
    map<int, vector<Real>> res;

    for (auto const& pair_layer_points: layer_points_map) {
        layer_sum = 0.0;
        int layer_size = pair_layer_points.second.size();
//        cout << "[Phi calculation] Layer: " << pair_layer_points.first << "| size: " << layer_size << endl;
        for (auto const& pointInLayer : pair_layer_points.second) {
            layer_sum += pointInLayer.v;
        }
        res[pair_layer_points.first].push_back(layer_sum);
        res[pair_layer_points.first].push_back(layer_sum / layer_size);
    }

    return res;
}

vector<Point> Analyzer::convertVtk2Points(const vector<Real> &vtk, const Point& box) {
    vector<Point> points;
    map<string, Point> xyz_point_map;
    int x = 0, y = 0, z = 0;

    // point representation of the vtk array
    for (auto &&value : vtk) {
        if (z == box.z) {
            z = z % box.z;
            y++;
        }
        if (y == box.y) {
            y = y % box.y;
            x++;
        }
        Point p = Point(x, y, z, value);
        xyz_point_map[to_string(x)+to_string(y)+to_string(z)] = p;
        points.push_back(p);
        z++;
    }
    return points;
}

map<int, vector<Point>> Analyzer::convertPoints2LayerPoints(const vector<Point>& points4converting) const {
    map<int, vector<Point>> res;

    map<int, vector<Point>> m = get_layer_point_map();  // layer <--> point | map

    for (auto const& pair_layer_points : m) {
//        cout << "[P2LP] Layer: " << pair_layer_points.first << "| size: " << pair_layer_points.second.size() << endl;
        for (auto const& pointInLayer : pair_layer_points.second) {
            //
            for (auto const& point4converting : points4converting) {
                if (point4converting == pointInLayer) {
                    res[pair_layer_points.first].push_back(point4converting);
                }
            }
        }
    }
    return res;
}

void Analyzer::updateVtk2PointsRepresentation(const vector<Real> &vtk, const Point &box) {
    pointsFromVtk = convertVtk2Points(vtk, box);
    layer_points_map = convertPoints2LayerPoints(pointsFromVtk);
}

vector<Point> Analyzer::Cube::get_cube_odd(int layers, const Point &core, int layer) {
    vector<Point> res;

    int start = -layers / 2;
    int end   = (layers / 2) + 1;

    // Template
    Point a1 = Point(1, 0, 0);
    Point a2 = Point(0, 1, 0);
    Point a3 = Point(0, 0, 1);

    for (int n1 = start; n1 < end; n1++) {
        for (int n2 = start; n2 < end; n2++) {
            for (int n3 = start; n3 < end ; n3++) {
                Point r = a1*n1 + a2*n2 + a3*n3;
                if (abs(r.x) == layer or abs(r.y) == layer or abs(r.z) == layer) {
                    res.push_back(r+core);
                }
            }
        }
    }
    return res;
}

map<int, vector<Point>> Analyzer::Cube::get_layer_point_map() const {
    return layer_point_map;
}

void Analyzer::Cube::construct_cube(int requested_layers, const Point &core) {
    cout << "[Cube] preparing appropriate cube around core: "<< core.to_string() << "..." << endl;
    // [stars] central node [Cleng] - shift[due to vtk file construction!]
    layer_point_map[0] = get_cube_odd(0, core - Point(1,1,1), 0);
    cout << "[Cube] layers done: " << 0;
    for (int current_layer = 1; current_layer < requested_layers; current_layer++){
        cout << " " << current_layer;
        // [stars] central node [Cleng] - shift[due to vtk file construction!]
        layer_point_map[current_layer] = get_cube_odd(current_layer*2+1, core - Point(1,1,1), current_layer);
    }
    cout << endl;

}
