#include <iostream>
#include <random>
#include <chrono>
#include <map>
using namespace std;
const int LATTICE_SIDE_LENGTH = 128;
const int COUPLING_STRENGTH = 1;
const int LONGITUDINAL_FIELD_STRENGTH = 0;
const int TRANSVERSE_FIELD_STANDARD_DEVIATION = 1;

struct Parameter {
    string type;
    int x1;
    int y1;
    int z1;
    Parameter (string type, int x1, int y1, int z1) : type(type), x1(x1), y1(y1), z1(z1) {}
};

struct Node : Parameter {
    Node (int x1, int y1, int z1) : Parameter("Node", x1, y1, z1) {}
};

struct Edge : Parameter {
    int x2;
    int y2;
    int z2;
    Edge (int x1, int y1, int z1, int x2, int y2, int z2) : 
        Parameter("Edge", x1, y1, z1), x2(x2), y2(y2), z2(z2) {}
};

multimap<double, Parameter> Parameters;

int main() {
	// generate the network randomly using 
    // gaussian distribution for fields centered at 0
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<double> distribution(0, TRANSVERSE_FIELD_STANDARD_DEVIATION);
    for (int i = 0; i < LATTICE_SIDE_LENGTH; i++) {
        for (int j = 0; j < LATTICE_SIDE_LENGTH; j++) {
            for (int k = 0; k < LATTICE_SIDE_LENGTH; k++) {
                Parameters.emplace(distribution(generator), Node(i, j, k));
                Parameters.emplace(COUPLING_STRENGTH, Edge(i, j, k, i, j, (k + 1) % LATTICE_SIDE_LENGTH));
                Parameters.emplace(COUPLING_STRENGTH, Edge(i, j, k, i, (j + 1) % LATTICE_SIDE_LENGTH, k));
                Parameters.emplace(COUPLING_STRENGTH, Edge(i, j, k, (i + 1) % LATTICE_SIDE_LENGTH, j, k));
            }
        };
    }
    
    while(!Parameters.empty()) {
        // run the first order rules
        if (Parameters.crbegin()->second.type == "Node") {
            // node decimation rules and print to file
            
        } else {
            // edge decimation rules 
            
        }
        // might need to be careful about if only 2 nodes are left
    }
}