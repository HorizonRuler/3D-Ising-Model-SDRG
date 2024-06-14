#include <fstream>
#include <random>
#include <chrono>
#include <queue>
#include <vector>
#include <utility>
#include <cstdio>
using namespace std;
const int LATTICE_SIDE_LENGTH = 2;
const int COUPLING_STRENGTH = 1;
const int LONGITUDINAL_FIELD_STRENGTH = 0;
const int TRANSVERSE_FIELD_STANDARD_DEVIATION = 1;

// inheritance structure for priority queue
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

// 3D adjacency list for each node pointing to nodes and edges in priority queue
pair<Node*, vector<Edge*>> adjacency_list[LATTICE_SIDE_LENGTH][LATTICE_SIDE_LENGTH][LATTICE_SIDE_LENGTH];
auto parameter_compare = [] (Parameter a, Parameter b) {return a.strength < b.strength;};
priority_queue<Parameter, vector<Parameter>, decltype(parameter_compare)> Parameters(parameter_compare);

int main() {
    // file output
    FILE *output_file = fopen("first_order_output.txt", "w");

	// generate the network randomly using gaussian distribution for fields centered at 0
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<double> distribution(0, TRANSVERSE_FIELD_STANDARD_DEVIATION);
    for (int i = 0; i < LATTICE_SIDE_LENGTH; i++) {
        for (int j = 0; j < LATTICE_SIDE_LENGTH; j++) {
            for (int k = 0; k < LATTICE_SIDE_LENGTH; k++) {
                Node n = Node(distribution(generator), i, j, k);
                Parameters.push(n);
                adjacency_list[i][j][k].first = &n;

                // mod to accound for periodic boundary conditions
                Edge e1 = Edge(COUPLING_STRENGTH, i, j, k, i, j, (k + 1) % LATTICE_SIDE_LENGTH);
                Edge e2 = Edge(COUPLING_STRENGTH, i, j, k, i, (j + 1) % LATTICE_SIDE_LENGTH, k);
                Edge e3 = Edge(COUPLING_STRENGTH, i, j, k, (i + 1) % LATTICE_SIDE_LENGTH, j, k);
                Parameters.push(e1);
                Parameters.push(e2);
                Parameters.push(e3);

                // maintain adjacency list with edges in both end nodes
                adjacency_list[i][j][k].second.push_back(&e1);
                adjacency_list[i][j][k].second.push_back(&e2);
                adjacency_list[i][j][k].second.push_back(&e3);
                adjacency_list[i][j][(k + 1) % LATTICE_SIDE_LENGTH].second.push_back(&e1);
                adjacency_list[i][(j + 1) % LATTICE_SIDE_LENGTH][k].second.push_back(&e2);
                adjacency_list[(i + 1) % LATTICE_SIDE_LENGTH][j][k].second.push_back(&e3);
            }
        }
    }
    
    // first order rules
    while(!Parameters.empty()) {
        if (Parameters.top().valid != false) {
            if (Parameters.top().type == "Node") {
                // node decimation and print to file
                fprintf(output_file, "Node %f: (%d, %d, %d)\n", Parameters.top().strength, Parameters.top().x1, Parameters.top().y1, Parameters.top().z1);

                for (Edge* p : adjacency_list[Parameters.top().x1][Parameters.top().y1][Parameters.top().z1].second) {
                    // increment neighboring node strengths by connected edge strengths
                    if (p->x1 == Parameters.top().x1 && p->y1 == Parameters.top().y1 && p->z1 == Parameters.top().z1) {
                        adjacency_list[p->x2][p->y2][p->z2].first->strength += p->strength;
                    } else {
                        adjacency_list[p->x1][p->y1][p->z1].first->strength += p->strength;
                    }

                    // remove connected edges from adjacency list
                    p->valid = false;
                }
                adjacency_list[Parameters.top().x1][Parameters.top().y1][Parameters.top().z1].second.clear();
            } else {
                // edge decimation and print to file
                Edge *e = (Edge*) &Parameters.top();
                fprintf(output_file, "Edge %f: (%d, %d, %d) <-> (%d, %d, %d)\n", e->strength, e->x1, e->y1, e->z1, e->x2, e->y2, e->z2);

                // combine the nodes, keeping 1st node and invalidating 2nd node
                adjacency_list[e->x1][e->y1][e->z1].first->strength += adjacency_list[e->x2][e->y2][e->z2].first->strength;
                adjacency_list[e->x2][e->y2][e->z2].first->valid = false;

                // combine the edges from 2nd node into 1st node
                vector<Edge*> edges_to_copy;
                for (Edge* p : adjacency_list[e->x2][e->y2][e->z2].second) {
                    for (Edge* q : adjacency_list[e->x1][e->y1][e->z1].second) {
                        if (p->x1 == q->x1 && p->y1 == q->y1 && p->z1 == q->z1 || p->x1 == q->x2 && p->y1 == q->y2 && p->z1 == q->z2 || p->x2 == q->x1 && p->y2 == q->y1 && p->z2 == q->z1 || p->x2 == q->x2 && p->y2 == q->y2 && p->z2 == q->z2) {
                            // maximum rule for overlapping edges
                            q->strength = max(q->strength, p->strength);
                            p->valid = false;
                            continue;
                        }
                    }
                    // copy over non-overlapping edges
                    if (p->x1 == e->x2 && p->y1 == e->y2 && p->z1 == e->z2) {
                        p->x1 = e->x1;
                        p->y1 = e->y1;
                        p->z1 = e->z1;
                    } else {
                        p->x2 = e->x1;
                        p->y2 = e->y1;
                        p->z2 = e->z1;
                    }
                    edges_to_copy.push_back(p);
                }
                adjacency_list[e->x2][e->y2][e->z2].second.clear();
                for (Edge* p : edges_to_copy)
                    adjacency_list[e->x1][e->y1][e->z1].second.push_back(p);
            }
        }
        Parameters.pop();
    }
    fclose(output_file);
    return 0;
}