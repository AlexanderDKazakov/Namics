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


Real Analyzer::calculateRg(vector<Real> vtk) {
    Real RG = 0.0;
    Real mass = 0.0;

    // FINISH IT!

    // calculating mass + point representation of the array
    for (auto &&value : vtk) {
        mass += value;

        if z == boxSize {
            z = z % boxSize
            y++
        }
        if y == boxSize {
            y = y % boxSize
            x++
        }

        p := NewPoint(x, y, z, f)
        points = append(points, *p)
        z++

    }

    cout << "Mass:" << mass << endl;
    // cm
	Real xtemp, ytemp, ztemp  = 0.0, 0.0, 0.0
    
    for (auto &&)

	for _, point := range points {
		xtemp += (float64(point.x) * point.v)
		ytemp += (float64(point.y) * point.v)
		ztemp += (float64(point.z) * point.v)
	}
	cm := Point{int(xtemp / mass), int(ytemp / mass), int(ztemp / mass), 0.0}
	fmt.Println("Center of mass: ", cm)

	var cdist Point
	for _, point := range points {
		cdist = point.SubPoint(cm)
		rg2 += point.v * (math.Pow(float64(cdist.x), 2) + math.Pow(float64(cdist.y), 2) + math.Pow(float64(cdist.z), 2))
	}
	rg2 /= mass
	fmt.Println("Rg2:", rg2)

    // need code
    return RG;
}

