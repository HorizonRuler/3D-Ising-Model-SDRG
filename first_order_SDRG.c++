#include <fstream>
#include <random>
#include <chrono>
#include <queue>
#include <unordered_set>
#include <utility>
#include <cstdio>
#include <math.h>
// #include <omp.h>
using namespace std;
const int LATTICE_SIDE_LENGTH = 64;
const int COUPLING_STRENGTH = 1;
const double COUPLING_NOISE_STANDARD_DEVIATION = 0.0000000001;
const double LONGITUDINAL_FIELD_STARTING_STANDARD_DEVIATION = 2;
const int LONGITUDINAL_FIELD_MEAN = 0;
const double STEP_SIZE = 0.01;
const int REPETITIONS_PER_STEP = 1000;
const int NUM_STEPS = 100;

// inheritance structure for priority queue
struct Parameter {
    double strength;
    int x1;
    int y1;
    int z1;
    string type;
    bool valid = true;
    Parameter (double strength, string type, int x1, int y1, int z1) : strength(strength), x1(x1), y1(y1), z1(z1), type(type) {}
};

struct Node : Parameter {
    vector<Node> domain;
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
        static const size_t shift = (size_t) log2(1 + sizeof(Edge));
        return (size_t)(e) >> shift;
    }
};

auto parameter_compare = [] (Parameter *a, Parameter *b) {return fabs(a->strength) < fabs(b->strength);};

double run_SDRG(double LONGITUDINAL_FIELD_STANDARD_DEVIATION = LONGITUDINAL_FIELD_STARTING_STANDARD_DEVIATION) {
    // 3D adjacency list for each node pointing to nodes and edges in priority queue
    pair<Node*, unordered_set<Edge*, EdgeHash>> adjacency_list[LATTICE_SIDE_LENGTH][LATTICE_SIDE_LENGTH][LATTICE_SIDE_LENGTH];
    priority_queue<Parameter*, vector<Parameter*>, decltype(parameter_compare)> parameters(parameter_compare);
    double final_domains[LATTICE_SIDE_LENGTH][LATTICE_SIDE_LENGTH][LATTICE_SIDE_LENGTH];

    // file output
    // FILE *decimation_steps_file = fopen("first_order_SDRG_steps.txt", "w");
    // FILE *final_domains_file = fopen("first_order_SDRG_domains.csv", "w");

	// generate the network randomly using gaussian distribution for fields centered at 0
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<double> fields(LONGITUDINAL_FIELD_MEAN, LONGITUDINAL_FIELD_STANDARD_DEVIATION);
    normal_distribution<double> couplings(COUPLING_STRENGTH, COUPLING_NOISE_STANDARD_DEVIATION);

    // initialize priority queue and adjacency list
    // omp_set_num_threads(4);
    // #pragma omp parallel for 
    for (int i = 0; i < LATTICE_SIDE_LENGTH; i++) {
        for (int j = 0; j < LATTICE_SIDE_LENGTH; j++) {
            for (int k = 0; k < LATTICE_SIDE_LENGTH; k++) {
                final_domains[i][j][k] = 0;

                // nodes first
                Node *n = new Node(fields(generator), i, j, k);
                parameters.push(n);
                adjacency_list[i][j][k].first = n;
                n->domain.push_back(*n);

                // mod to accound for periodic boundary conditions
                Edge *e1 = new Edge(couplings(generator), i, j, k, i, j, (k + 1) % LATTICE_SIDE_LENGTH);
                Edge *e2 = new Edge(couplings(generator), i, j, k, i, (j + 1) % LATTICE_SIDE_LENGTH, k);
                Edge *e3 = new Edge(couplings(generator), i, j, k, (i + 1) % LATTICE_SIDE_LENGTH, j, k);
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
                Node *n = (Node*) parameters.top();
                // fprintf(decimation_steps_file, "Node strength %f: (%d, %d, %d)\n", n->strength, n->x1, n->y1, n->z1);
                // fprintf(decimation_steps_file, "Original Nodes in domain:\n");
                for (Node m : n->domain) {
                    // fprintf(decimation_steps_file, "Original strength %f: (%d, %d, %d)\n", n.strength, n.x1, n.y1, n.z1);
                    final_domains[m.x1][m.y1][m.z1] = n->strength;
                }
                // fprintf(decimation_steps_file, "\n");

                vector<Node*> nodes_to_copy;
                int node_sign = (adjacency_list[n->x1][n->y1][n->z1].first->strength > 0) - (adjacency_list[n->x1][n->y1][n->z1].first->strength < 0);
                for (Edge* p : adjacency_list[n->x1][n->y1][n->z1].second) {
                    // increment or decrement neighboring node strengths by connected edge strengths based on sign of removed node and maintain heap by creating new nodes and invalidating old ones
                    int x, y, z;
                    if (p->x1 == n->x1 && p->y1 == n->y1 && p->z1 == n->z1) {
                        x = p->x2;
                        y = p->y2;
                        z = p->z2;
                    } else {
                        x = p->x1;
                        y = p->y1;
                        z = p->z1;
                    }
                    Node *m = new Node(adjacency_list[x][y][z].first->strength + p->strength * node_sign, x, y, z);
                    m->domain = adjacency_list[x][y][z].first->domain;
                    adjacency_list[x][y][z].first->valid = false;
                    adjacency_list[x][y][z].first = m;
                    nodes_to_copy.push_back(m);
                    adjacency_list[x][y][z].second.erase(p);

                    // remove connected edges from adjacency list
                    p->valid = false;
                }

                // remove node and associated from adjacency list, will be popped immediately from priority queue so no need to invalidate
                adjacency_list[n->x1][n->y1][n->z1].second.clear();
                adjacency_list[n->x1][n->y1][n->z1].first = nullptr;

                // clean up memory and pop
                delete n;
                parameters.pop();

                // update priority queue after popping
                for (Node* m : nodes_to_copy)
                    parameters.push(m);
            } else {
                // edge decimation, print to file
                Edge *e = (Edge*) parameters.top();
                // fprintf(decimation_steps_file, "Edge strength %f: (%d, %d, %d) <-> (%d, %d, %d)\n", e->strength, e->x1, e->y1, e->z1, e->x2, e->y2, e->z2);

                // new node with combined strength
                Node *n = new Node(adjacency_list[e->x1][e->y1][e->z1].first->strength + adjacency_list[e->x2][e->y2][e->z2].first->strength, e->x1, e->y1, e->z1);

                // combine domains for output
                n->domain = adjacency_list[e->x1][e->y1][e->z1].first->domain;
                n->domain.insert(n->domain.end(), adjacency_list[e->x2][e->y2][e->z2].first->domain.begin(), adjacency_list[e->x2][e->y2][e->z2].first->domain.end());

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
                        if ((p->x1 == q->x1 && p->y1 == q->y1 && p->z1 == q->z1) || (p->x1 == q->x2 && p->y1 == q->y2 && p->z1 == q->z2) || (p->x2 == q->x1 && p->y2 == q->y1 && p->z2 == q->z1) || (p->x2 == q->x2 && p->y2 == q->y2 && p->z2 == q->z2)) {
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
                delete e;
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

    // Calculate magnetic moment for this run
    int imbalanced_field_count = 0;
    for (int i = 0; i < LATTICE_SIDE_LENGTH; i++) {
        for (int j = 0; j < LATTICE_SIDE_LENGTH; j++) {
            // fprintf(final_domains_file, "%f", final_domains[i][j][0]);
            for (int k = 1; k < LATTICE_SIDE_LENGTH; k++) {
                // print final domains to file
                // fprintf(final_domains_file, ", %f", final_domains[i][j][k]);
                if (final_domains[i][j][k] > 0)
                    imbalanced_field_count++;
                else
                    imbalanced_field_count--;
            }
            // fprintf(final_domains_file, "\n");
        }
    }

    // close files
    // fclose(final_domains_file);
    // fclose(decimation_steps_file);

    return abs((double) imbalanced_field_count / (LATTICE_SIDE_LENGTH * LATTICE_SIDE_LENGTH * LATTICE_SIDE_LENGTH));
}

int main() {
    // open file for output
    FILE *magnetization_moment_file = fopen("first_order_SDRG_magnetization_moment.txt", "w");

    // run SDRG for each step
    for (int i = 0; i < NUM_STEPS; i++) {
        double total_magnetization_moment = 0;
        for (int j = 0; j < REPETITIONS_PER_STEP; j++) {
            total_magnetization_moment += run_SDRG(LONGITUDINAL_FIELD_STARTING_STANDARD_DEVIATION + i * STEP_SIZE);
        }
        fprintf(magnetization_moment_file, "h st. dev. %f : average moment %f\n", LONGITUDINAL_FIELD_STARTING_STANDARD_DEVIATION + i * STEP_SIZE, total_magnetization_moment / REPETITIONS_PER_STEP);
    }

    // close file
    fclose(magnetization_moment_file);
    return 0;
}