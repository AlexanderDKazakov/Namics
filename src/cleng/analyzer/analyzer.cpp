#include "analyzer.h"

using namespace std;

Real Analyzer::calculateRe(map<int, vector<int>> pivot_arm_nodes, map<int, vector<std::shared_ptr<Node>>> nodes_map) {
    // pivot_arm_nodes --> gives id_nodes for central and last one
    // nodes_id        --> gives possibility to get xyz coordinate of explicit segment (node in point representation)
    Real RE = 0.0;
    int arm_number = 0;
    for (auto &&pair_pivot : pivot_arm_nodes) {  
        auto id_central_node = pair_pivot.second.begin()[0];
        auto id_last_node    = pair_pivot.second.back();        
        
        auto central_node = nodes_map[id_central_node].data()->get()->point();
        auto last_node    = nodes_map[id_last_node].data()->get()->point();
        Real distance = central_node.distance(last_node);
        RE += distance;
        cout << "[C]:" << id_central_node <<  " | [L]:" << id_last_node << " || " << distance << endl;
        arm_number++;
    }
    return RE /= arm_number;
}


Real Analyzer::calculateRg(const vector<Real>& vtk, const Point& box) {
    Real RG2 = 0.0;
    Real mass = 0.0;
    vector<Point> points;
    int x = 0, y = 0, z = 0;

    cout << "Box: " << box.to_string() << endl;

    // calculating mass + point representation of the array
    for (auto &&value : vtk) {
        mass += value;
        if (z == box.z) {
            z = z % box.z;
            y++;
        }
        if (y == box.y) {
            y = y % box.y;
            x++;
        }
        Point p = Point(x, y, z, value);
        points.push_back(p);
        z++;
    }
    cout << "Mass:" << mass << endl;
    // cm
	Real xtemp = 0.0, ytemp = 0.0, ztemp = 0.0;
    for (auto &&point : points) {
        xtemp += point.x * point.v;
        ytemp += point.y * point.v;
        ztemp += point.z * point.v;
    }
    Point cm = Point(int(xtemp / mass), int(ytemp/mass), int(ztemp/mass));
    cout << "Center of mass " << cm.to_string() << endl;
    Point cdist;
    for (auto &&point : points) {
        cdist = point-cm;
        RG2 += point.v * (pow(cdist.x, 2) + pow(cdist.y, 2) + pow(cdist.z, 2));
    }
    return RG2 / mass;
}

//vector<Point> Cube::get_cube_odd(int layers, const Point& core, int layer) {
//    vector<Point> res;
//
//    int start = layers / 2;
//    int end   = (layers / 2) + 1;
//
//    // Template
//    Point a1 = Point(1, 0, 0);
//    Point a2 = Point(0, 1, 0);
//    Point a3 = Point(0, 0, 1);
//
//    for (int n1 = start; n1 < end; n1++) {
//        for (int n2 = start; n2 < end; n2++) {
//            for (int n3 = start; n3 < end ; n3++) {
//                Point r = a1*n1 + a2*n2 + a3*n3 + core;
//                if (r.x == layer or r.y == layer or r.z == layer) {
//                    res.push_back(r);
//                }
//            }
//        }
//    }
//
//    return res;
//}

//map<int, vector<Point>> Cube::get_layer_point_map() {
//    return layer_point_map;
//}

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
                Point r = a1*n1 + a2*n2 + a3*n3 + core;
//                cout << "r: " << r.to_string() << endl;
                if ( abs(r.x) == core.x + layer or abs(r.y) == core.y + layer or abs(r.z) == core.z + layer) {
//                    cout << "ADDED!" << endl;
                    res.push_back(r);
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
    cout << "[Cube] init..." << endl;
    layer_point_map[0] = get_cube_odd(0, core, 0);
    cout << "[Cube] preparing appropriate cube..." << endl;
    for (int current_layer = 1; current_layer < requested_layers; current_layer++){
        cout << "[Cube] layer:"<< current_layer << endl;
        layer_point_map[current_layer] = get_cube_odd(current_layer*2+1, core, current_layer);
    }

}
