//
// Created by nick434434 on 14-10-19.
//
#include "prob_stuff.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

typedef boost::adjacency_list<> Graph;

using std::vector;


void s() {
    Graph g(10);
}


vector<vector<int>> generate_GRG(long n, const std::function<long()>& gen) {
    vector<vector<int>> G(n);

    Graph g(n);

    vector<long> w;
    std::generate_n(std::back_inserter(w), n, gen);

    long ln = 0;
    ln = std::accumulate(w.begin(), w.end(), ln);

    long tmp;
    double prob;
    for (long i = 0; i < n; ++i) {
        G[i] = vector<int>(n);
        for (long j = i + 1; j < n; ++j) {
            tmp = w[i] * w[j];
            prob = (double)tmp / (double)(ln + tmp);
            G[i][j] = edge_decider(prob);
            if (G[i][j])
                boost::add_edge(i, j, g);
        }
    }

    std::ofstream f("graph.dot");
    boost::write_graphviz(f, g);
    f.close();

    return G;
}