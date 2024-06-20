#include <fstream>
#include <random>
#include <chrono>
#include <queue>
#include <unordered_set>
#include <utility>
#include <cstdio>
using namespace std;
const int LATTICE_SIDE_LENGTH = 3;
const int COUPLING_STRENGTH = 1;
const int LONGITUDINAL_FIELD_STRENGTH = 0;
const double TRANSVERSE_FIELD_STANDARD_DEVIATION = 0.25;
const int TRANSVERSE_FIELD_MEAN = 0;

// inheritance structure for priority queue
struct Parameter {
    double strength;
    int x1;
    int y1;
    int z1;
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

// https://stackoverflow.com/questions/20953390/what-is-the-fastest-hash-function-for-pointers
struct EdgeHash {
    size_t operator() (const Edge* e) const {
        static const size_t shift = (size_t)log2(1 + sizeof(Edge));
        return (size_t)(e) >> shift;
    }
};

// 3D adjacency list for each node pointing to nodes and edges in priority queue
pair<Node*, unordered_set<Edge*, EdgeHash>> adjacency_list[LATTICE_SIDE_LENGTH][LATTICE_SIDE_LENGTH][LATTICE_SIDE_LENGTH];
auto parameter_compare = [] (Parameter *a, Parameter *b) {return a->strength < b->strength;};
priority_queue<Parameter*, vector<Parameter*>, decltype(parameter_compare)> parameters(parameter_compare);

int main() {
    // file output
    FILE *output_file = fopen("first_order_output.txt", "w");

	// generate the network randomly using gaussian distribution for fields centered at 0
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<double> distribution(TRANSVERSE_FIELD_MEAN, TRANSVERSE_FIELD_STANDARD_DEVIATION);

    // initialize priority queue and adjacency list
    for (int i = 0; i < LATTICE_SIDE_LENGTH; i++) {
        for (int j = 0; j < LATTICE_SIDE_LENGTH; j++) {
            for (int k = 0; k < LATTICE_SIDE_LENGTH; k++) {
                // nodes first
                Node *n = new Node(distribution(generator), i, j, k);
                parameters.push(n);
                adjacency_list[i][j][k].first = n;

                // mod to accound for periodic boundary conditions
                Edge *e1 = new Edge(COUPLING_STRENGTH, i, j, k, i, j, (k + 1) % LATTICE_SIDE_LENGTH);
                Edge *e2 = new Edge(COUPLING_STRENGTH, i, j, k, i, (j + 1) % LATTICE_SIDE_LENGTH, k);
                Edge *e3 = new Edge(COUPLING_STRENGTH, i, j, k, (i + 1) % LATTICE_SIDE_LENGTH, j, k);
                parameters.push(e1);
                parameters.push(e2);
                parameters.push(e3);

                // maintain adjacency list with edges in both end nodes
                adjacency_list[i][j][k].second.insert(e1);
                adjacency_list[i][j][k].second.insert(e2);
                adjacency_list[i][j][k].second.insert(e3);
                adjacency_list[i][j][(k + 1) % LATTICE_SIDE_LENGTH].second.insert(e1);
                adjacency_list[i][(j + 1) % LATTICE_SIDE_LENGTH][k].second.insert(e2);
                adjacency_list[(i + 1) % LATTICE_SIDE_LENGTH][j][k].second.insert(e3);
            }
        }
    }
    
    // first order rules
    while(!parameters.empty()) {
        if (parameters.top()->valid != false) {
            if (parameters.top()->type == "Node") {
                // node decimation, print to file
                fprintf(output_file, "Node strength %f: (%d, %d, %d)\n", parameters.top()->strength, parameters.top()->x1, parameters.top()->y1, parameters.top()->z1);

                vector<Node*> nodes_to_copy;
                for (Edge* p : adjacency_list[parameters.top()->x1][parameters.top()->y1][parameters.top()->z1].second) {
                    // increment neighboring node strengths by connected edge strengths and maintain heap by creating new nodes and invalidating old ones
                    if (p->x1 == parameters.top()->x1 && p->y1 == parameters.top()->y1 && p->z1 == parameters.top()->z1) {
                        Node *n = new Node(adjacency_list[p->x2][p->y2][p->z2].first->strength + p->strength, p->x2, p->y2, p->z2);
                        adjacency_list[p->x2][p->y2][p->z2].first->valid = false;
                        adjacency_list[p->x2][p->y2][p->z2].first = n;
                        nodes_to_copy.push_back(n);
                        adjacency_list[p->x2][p->y2][p->z2].second.erase(p);
                    } else {
                        Node *n = new Node(adjacency_list[p->x1][p->y1][p->z1].first->strength + p->strength, p->x1, p->y1, p->z1);
                        adjacency_list[p->x1][p->y1][p->z1].first->valid = false;
                        adjacency_list[p->x1][p->y1][p->z1].first = n;
                        nodes_to_copy.push_back(n);
                        adjacency_list[p->x1][p->y1][p->z1].second.erase(p);
                    }

                    // remove connected edges from adjacency list
                    p->valid = false;
                }

                // remove node and associated from adjacency list, will be popped immediately from priority queue so no need to invalidate
                adjacency_list[parameters.top()->x1][parameters.top()->y1][parameters.top()->z1].second.clear();adjacency_list[parameters.top()->x1][parameters.top()->y1][parameters.top()->z1].first = nullptr;

                // clean up memory and pop
                delete parameters.top();
                parameters.pop();

                // update priority queue after popping
                for (Node* n : nodes_to_copy)
                    parameters.push(n);
            } else {
                // edge decimation, print to file
                Edge *e = (Edge*) parameters.top();
                fprintf(output_file, "Edge strength %f: (%d, %d, %d) <-> (%d, %d, %d)\n", e->strength, e->x1, e->y1, e->z1, e->x2, e->y2, e->z2);

                // new node with combined strength
                Node *n = new Node(adjacency_list[e->x1][e->y1][e->z1].first->strength + adjacency_list[e->x2][e->y2][e->z2].first->strength, e->x1, e->y1, e->z1);

                // invalidate 1st and 2nd node
                adjacency_list[e->x1][e->y1][e->z1].first->valid = false;
                adjacency_list[e->x2][e->y2][e->z2].first->valid = false;

                // update 1 adjacency list to point to new node and other to null
                adjacency_list[e->x1][e->y1][e->z1].first = n;
                adjacency_list[e->x2][e->y2][e->z2].first = nullptr;

                // remove connecting edge from adjacency list, will be popped immediately from priority queue so no need to invalidate
                adjacency_list[e->x1][e->y1][e->z1].second.erase(e);
                adjacency_list[e->x2][e->y2][e->z2].second.erase(e);

                // combine the edges from 2nd node into 1st node
                vector<Edge*> edges_to_copy;
                for (Edge* p : adjacency_list[e->x2][e->y2][e->z2].second) {
                    bool overlap = false;
                    for (Edge* q : adjacency_list[e->x1][e->y1][e->z1].second) {
                        if (p->x1 == q->x1 && p->y1 == q->y1 && p->z1 == q->z1 || p->x1 == q->x2 && p->y1 == q->y2 && p->z1 == q->z2 || p->x2 == q->x1 && p->y2 == q->y1 && p->z2 == q->z1 || p->x2 == q->x2 && p->y2 == q->y2 && p->z2 == q->z2) {
                            // maximum rule for overlapping edges (doesn't matter because both are 1 for now)
                            overlap = true;

                            // get rid of one of the edges
                            p->valid = false;
                            if (p->x1 == e->x2 && p->y1 == e->y2 && p->z1 == e->z2)
                                adjacency_list[p->x2][p->y2][p->z2].second.erase(p);
                            else 
                                adjacency_list[p->x1][p->y1][p->z1].second.erase(p);
                            break;
                        }
                    }

                    // skip moving if overlapping
                    if (overlap) 
                        continue;

                    // move over non-overlapping edges
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

                // remove all edges from 2nd node
                adjacency_list[e->x2][e->y2][e->z2].second.clear();

                // update 1st adjacency list after all done
                for (Edge* p : edges_to_copy)
                    adjacency_list[e->x1][e->y1][e->z1].second.insert(p);

                // clean up memory and pop
                delete parameters.top();
                parameters.pop();

                // update priority queue after popping
                parameters.push(n);
            }
        } else {
            // clean up memory and pop
            delete parameters.top();
            parameters.pop();
        }
    }
    fclose(output_file);
    return 0;
}