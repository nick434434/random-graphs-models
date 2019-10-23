//
// Created by nick434434 on 19-10-19.
//


#include "conf_model.h"
#include "timing.h"
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <limits>

typedef boost::exterior_vertex_property<Graph, double> DistanceProperty;
typedef DistanceProperty::matrix_type DistanceMatrix;
typedef DistanceProperty::matrix_map_type DistanceMatrixMap;

typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;

using std::unordered_set;
using std::vector;
using std::pair;


// template<class Container>
ConfigurationModel::ConfigurationModel(const vector<long>& new_degrees) {
    // static_assert(has_const_iterator<Container>::value, "Used type does not have const iterator\n");
    degrees = vector<long>(new_degrees.begin(), new_degrees.end());
    n = degrees.size();
}


ConfigurationModel::~ConfigurationModel() {
    // Clearing the Graph object g
    clear_realization();
    degrees.erase(degrees.begin(), degrees.end());
}


void ConfigurationModel::make_half_edges() {
    he = vector<pair<long, long>>();

    long half_edge_counter = 0;
    for (size_t i = 0; i < n; ++i)
        for (long j = 0; j < degrees[i]; ++j) {
            // First elem will be number of connected half-edge, second is vertex for which half-edge is created
            he.push_back({half_edge_counter++, i});
        }
}


void ConfigurationModel::connect_half_edges() {

    std::uniform_int_distribution<unsigned long> picker(0, he.size() - 1);

    // Optimized container for search - just what we need
    // unordered_set<long> used;
    unordered_set<long> available;
    for (size_t i = 0; i < he.size(); ++i)
        available.insert(i);
    long connection;

    for (size_t i = 0; i < he.size(); ++i) {
        while (i < he.size() && available.find(i) == available.end())
            ++i;
        if (i >= he.size())
            break;
        // Marking current half-edge as used
        available.erase(available.find(i));

        // Picking not-used half-edge
        connection = (long)(std::trunc(uniform_distribution(generator_uniform) * (available.size() - 1)));
        auto it = available.begin();
        for (long j = 0; j < connection; ++j)
            it++;
        connection = *it;
        // Marking picked half-edge as used
        available.erase(available.find(connection));

        // Creating an edge in realization
        boost::add_edge(he[i].second, he[connection].second, 1,g);

        // Remembering which half-edges are paired in case we need this info
        he[connection].first = i;
        he[i].first = connection;
    }

}


void ConfigurationModel::clear_realization() {
    boost::graph_traits<Graph>::vertex_iterator vi, vi_end, next;
    std::tie(vi, vi_end) = boost::vertices(g);
    for (next = vi; next != vi_end; vi = next) {
        ++next;
        boost::remove_vertex(*vi, g);
    }

    he.erase(he.begin(), he.end());
}


void ConfigurationModel::compute_distance(bool use_johnson) {
    DistanceMatrix distances(n);
    DistanceMatrixMap dm(distances, g);

    WeightMap wm = boost::get(boost::edge_weight, g);

    // std::string message = "Distances computed";
    // timing::start_local_clock();
    if (use_johnson)
        boost::johnson_all_pairs_shortest_paths(g, dm, boost::weight_map(wm));
    else
        boost::floyd_warshall_all_pairs_shortest_paths(g, dm, boost::weight_map(wm));
    // timing::reset_local_clock(message);

    // message = "Maximum distance found";
    distance = 0;
    bool f;
    long count = 0;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j) {
            f = std::abs(distances[i][j] - std::numeric_limits<typeof(distance)>::max()) < 1e-5;
            distance += f ? 0 : distances[i][j];
            count += 1 - (int)f;
        }
            // distance = std::max(distance, std::abs(distances[i][j] - std::numeric_limits<typeof(distance)>::max()) < 1e-5 ? 0 : distances[i][j]);
    distance /= count;
    // timing::reset_local_clock(message);
    // distance = *std::max_element(d_max.begin(), d_max.end());
    // for (long i = 0; i < n; ++i)
    //     max_d = std::max(*std::max_element(d[i].begin(), d[i].end()), max_d);
}


void ConfigurationModel::get_graphviz(const std::string& dot_fname) {
    std::ofstream f(dot_fname);
    boost::write_graphviz(f, g);
    f.close();
}
