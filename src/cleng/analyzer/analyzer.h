
#include "../../namics.h"
#include "../nodes/simple_node.h"
#include "../nodes/monolit.h"
#include "../nodes/point.h"
#include <map>
#include <vector>
#include <memory>

using namespace std;

class Analyzer {
public:

    Real calculateRe(map<int, vector<int>> pivot_arm_nodes, map<int, vector<std::shared_ptr<Node>>> nodes_map);
    Real calculateRg(vector<Real> vtk);

    
};