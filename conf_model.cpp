//
// Created by sv2019 on 19-10-19.
//


#include "conf_model.h"
#include <vector>
#include <unordered_set>
#include <algorithm>

using std::unordered_set;
using std::vector;
using std::pair;


// Helper to determine whether there's a const_iterator for T.
template<typename T>
struct has_const_iterator
{
private:
    template<typename C> static char test(typename C::const_iterator*);
    template<typename C> static int test(...);
public:
    enum { value = sizeof(test<T>(0)) == sizeof(char) };
};


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
        boost::add_edge(he[i].second, he[connection].second, g);

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


void ConfigurationModel::get_graphviz(const std::string& dot_fname) {
    std::ofstream f(dot_fname);
    boost::write_graphviz(f, g);
    f.close();
}
