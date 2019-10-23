//
// Created by nick434434 on 19-10-19.
//
#pragma once

#ifndef NUMERICALASSIGNMENTPS2_CONF_MODEL_H
#define NUMERICALASSIGNMENTPS2_CONF_MODEL_H

#include "prob_stuff.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

using std::pair;

// multisetS for OutEdgeList is used to ensure multiple edges are possible
typedef boost::adjacency_list<boost::multisetS, boost::vecS, boost::directedS> Graph;
// typedef boost:: Matrix;

class ConfigurationModel {
private:
    // Base variable for storing generated realization
    Graph g;
    // Base variable for storing degrees of vertices
    vector<long> degrees;
    // Number of vertices
    long n = 0;
    // Temporary variable for constructing realization
    vector<pair<long, long>> he = vector<pair<long, long>>();

public:
    // Graph distance
    long distance = 0;

    ConfigurationModel() = default;

    explicit ConfigurationModel(long n): n(n), degrees(vector<long>(n, 0)) {}

    // template<class Container>
    explicit ConfigurationModel(const vector<long>& new_degrees);

    ~ConfigurationModel();

    /*
     * Each half-edge has vertex attachment as 2nd elem of pair and reserved first elem for number of paired half-edge
     */
    void make_half_edges();

    /*
     * Connecting half-edges by creating edges in @field g and filling in numbers of paired half-edges
     */
    void connect_half_edges();

    void clear_realization();

    void compute_distance();

    void get_graphviz(const std::string& dot_fname);

};


#endif //NUMERICALASSIGNMENTPS2_CONF_MODEL_H
