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
using std::pair;


void s() {
    Graph g(10);
}


vector<vector<int>> generate_GRG(long n, long max_weight, double tau_subs) {
    vector<vector<int>> G(n);

    Graph g(n);

    if (tau_subs > 1) {
        tau = tau_subs;
        dist = boost::math::pareto_distribution<>(tau, tau);
    }

    vector<long> w = generate_n_pareto(1, max_weight, n);

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


vector<pair<long, long>> generate_GRG_edge_pairs(long n, const vector<long>& w) {
    vector<pair<long, long>> G;

#ifdef VISUALIZE_GRAPH
    Graph g(n);
#endif

    long ln = 0;
    ln = std::accumulate(w.begin(), w.end(), ln);

    long tmp;
    double prob;
    for (long i = 0; i < n; ++i) {
        for (long j = i + 1; j < n; ++j) {
            tmp = w[i] * w[j];
            prob = (double)tmp / (double)(ln + tmp);
            if (edge_decider(prob)) {
#ifdef VISUALIZE_GRAPH
                boost::add_edge(i, j, g);
#endif
                G.push_back({i, j});
            }
        }
    }

#ifdef VISUALIZE_GRAPH
    std::ofstream f("graph.dot");
    boost::write_graphviz(f, g);
    f.close();
#endif

    return G;
}


vector<pair<long, long>> generate_GRG_edge_pairs(long n, long max_weight, double tau_subs = 0) {
    if (tau_subs > 1) {
        tau = tau_subs;
        dist = boost::math::pareto_distribution<>(tau, tau);
    }

    vector<long> w = generate_n_pareto(1, max_weight, n);

    return generate_GRG_edge_pairs(n, w);
}


vector<vector<pair<long, long>>> generate_m_GRG_edge_pairs(long n, long m, long max_weight, double tau_subs = 0) {
    vector<vector<pair<long, long>>> res;

    for (long i = 0; i < m; ++i)
        res.push_back(generate_GRG_edge_pairs(n, max_weight, tau_subs));

    return res;
}


vector<vector<pair<long, long>>> generate_m_GRG_edge_pairs_same_weights(long n, long m, long max_weight, vector<long>& weights, double tau_subs = 0) {
    vector<vector<pair<long, long>>> res;

    weights = generate_n_pareto(1, max_weight, n);

    for (long i = 0; i < m; ++i) {
        res.push_back(generate_GRG_edge_pairs(n, weights));
        std::cout << i + 1 << " constructed\n";
    }

    return res;
}


vector<long> count_degrees(const vector<pair<long, long>>& G) {
    vector<long> degrees;

    long max_v;
    for (auto e: G) {
        max_v = std::max(e.first, e.second);
        if (max_v + 1 > degrees.size())
            for (long i = degrees.size(); i < max_v + 1; ++i)
                degrees.push_back(0);
        degrees[e.first]++;
        degrees[e.second]++;
    }

    return degrees;
}
