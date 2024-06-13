#include <fstream>
#include <random>
#include <chrono>
#include <queue>
#include <vector>
#include <utility>
using namespace std;
const int LATTICE_SIDE_LENGTH = 128;
const int COUPLING_STRENGTH = 1;
const int LONGITUDINAL_FIELD_STRENGTH = 0;
const int TRANSVERSE_FIELD_STANDARD_DEVIATION = 1;

struct Parameter {
    int x1;
    int y1;
    int z1;
    double strength;
    string type;
    bool valid = true;
    Parameter (double strength, string type, int x1, int y1, int z1) : strength(strength), type(type), x1(x1), y1(y1), z1(z1) {}
};

struct Node : Parameter {
    Node (double strength, int x1, int y1, int z1) : Parameter(strength, "Node", x1, y1, z1) {}
};

struct Edge : Parameter {
    int x2;
    int y2;
    int z2;
    Edge (double strength, int x1, int y1, int z1, int x2, int y2, int z2) : Parameter(strength, "Edge", x1, y1, z1), x2(x2), y2(y2), z2(z2) {}
};

// 3D adjacency list of vertices pointing to nodes and edges in priority queue
pair<Parameter*, vector<Parameter*>> adjacency_list[LATTICE_SIDE_LENGTH][LATTICE_SIDE_LENGTH][LATTICE_SIDE_LENGTH];
priority_queue<Parameter> Parameters;

int main() {
	// generate the network randomly using 
    // gaussian distribution for fields centered at 0
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<double> distribution(0, TRANSVERSE_FIELD_STANDARD_DEVIATION);
    for (int i = 0; i < LATTICE_SIDE_LENGTH; i++) {
        for (int j = 0; j < LATTICE_SIDE_LENGTH; j++) {
            for (int k = 0; k < LATTICE_SIDE_LENGTH; k++) {
                Node n = Node(distribution(generator), i, j, k);
                Edge e1 = Edge(COUPLING_STRENGTH, i, j, k, i, j, (k + 1) % LATTICE_SIDE_LENGTH);
                Edge e2 = Edge(COUPLING_STRENGTH, i, j, k, i, (j + 1) % LATTICE_SIDE_LENGTH, k);
                Edge e3 = Edge(COUPLING_STRENGTH, i, j, k, (i + 1) % LATTICE_SIDE_LENGTH, j, k);
                Parameters.push(n);
                Parameters.push(e1);
                Parameters.push(e2);
                Parameters.push(e3);
                adjacency_list[i][j][k].first = &n;
                adjacency_list[i][j][k].second.push_back(&e1);
                adjacency_list[i][j][k].second.push_back(&e2);
                adjacency_list[i][j][k].second.push_back(&e3);
                adjacency_list[i][j][(k + 1) % LATTICE_SIDE_LENGTH].second.push_back(&e1);
                adjacency_list[i][(j + 1) % LATTICE_SIDE_LENGTH][k].second.push_back(&e2);
                adjacency_list[(i + 1) % LATTICE_SIDE_LENGTH][j][k].second.push_back(&e3);
            }
        };
    }
    
    while(!Parameters.empty()) {
        // run the first order rules
        if (Parameters.top().type == "Node") {
            // node decimation rules and print to file
            
        } else {
            // edge decimation rules 
            
        }
        // might need to be careful about if only 2 nodes are left
    }
}