#include "erdos_renyi.h"
#include "grg.h"
#include "matplotlibcpp.h"
#include "prob_stuff.h"
#include <iostream>
#include <map>
#include <algorithm>

// Use this for saving .dot Graphviz representation
// #define VISUALIZE_GRAPH

namespace plt = matplotlibcpp;

using std::map;
using std::cout;
using std::endl;


void example() {
    vector<vector<double>> x, y, z;
    for (double i = -5; i <= 5 + 1e-10;  i += 0.25) {
        vector<double> x_row, y_row, z_row;
        for (double j = -5; j <= 5 + 1e-10; j += 0.25) {
            x_row.push_back(i);
            y_row.push_back(j);
            z_row.push_back(std::sin(std::hypot(i, j)));
        }
        x.push_back(x_row);
        y.push_back(y_row);
        z.push_back(z_row);
    }

    plt::plot_surface(x, y, z);
    plt::show();
}


void construct_plot(int N_trials, double lmbd, int N_start, int N_step, int num_iters, std::string folder) {
    vector<pair<vector<int>, vector<long>>> results(N_trials);

    for (int i = 0; i < N_trials; ++i)
        results[i] = generate_for_single_lambda(lmbd, N_start, N_step, num_iters), cout << "Trial " << i + 1 << "!\n";

    vector<map<long, long>> counts;

    for (int i = 0; i < num_iters; ++i) {
        map<long, long> values_and_counts;
        for (int j = 0; j < N_trials; ++j) {
            auto smth = results[j].second[i];
            if (values_and_counts.find(smth) == values_and_counts.end())
                values_and_counts[smth] = 1;
            else
                values_and_counts[smth]++;
        }
        counts.push_back(values_and_counts);
    }

    vector<long> x, y;
    vector<int> Ns = results[0].first;
    if (folder != "/")
        folder = "/" + folder;
    for (int i = 0; i < num_iters; ++i) {
        long minx = INT_MAX, maxx = 0;
        for (auto map_elem: counts[i]) {
            x.push_back(map_elem.first);
            y.push_back(map_elem.second);
            if (map_elem.first < minx)
                minx = map_elem.first;
            else if (map_elem.first > maxx)
                maxx = map_elem.first;
        }
        vector<double> poisson_y = generate_poisson_function(lmbd * lmbd * lmbd / 6, minx, maxx, N_trials);
        vector<double> poisson_x(maxx - minx + 1);
        for (long j = 0; j < maxx - minx + 1; ++j)
            poisson_x[j] = minx + j;
        plt::named_plot("Experiments", x, y);
        plt::named_plot("Actual Poisson", poisson_x, poisson_y);
        plt::title("Graphs of size " + std::to_string(Ns[i]));
        plt::xlabel("# of triangles");
        plt::ylabel("# of occurences in experiments");
        cout << "../plots" + folder + "/size_" + std::to_string(Ns[i]) + "_lambda_" + std::to_string((int)(lmbd * 10)/10) + "_trials_" +
                std::to_string(N_trials) + ".png" << endl;
        plt::save("../plots" + folder + "/size_" + std::to_string(Ns[i]) + "_lambda_" + std::to_string((int)(lmbd * 10)/10) + "_trials_" +
                  std::to_string(N_trials) + ".png");
        x.clear();
        y.clear();
        plt::clf();
    }

}


void construct_plot(int N_trials, double lmbd, int N_start, int N_step, int num_iters) {
    construct_plot(N_trials, lmbd, N_start, N_step, num_iters, "/");
}


void test1() {
    int N = 10;

    cout << "Lambda | # of triangles | # of edges | expected # of edges" << endl;

    vector<vector<int>> G;
    for (int i = 0; i < 20; ++i) {
        G = generate_ER(N, 1. + 0.2 * i);
        cout << 1. + 0.2 * i << " | " << count_triangles(G) << " | " << count_edges(G) << " | " << (1. + 0.2 * i) * (N - 1) / 2 << endl;
    }
}


void test2() {
    auto results = generate_for_single_n(100, 1., 0.2, 10);

    plt::plot(results.first, results.second);
    // TODO: make legend function work. plt::legend({"delta as function of lambda"});
    plt::xlabel("lambda");
    plt::ylabel("# of triangles");
    plt::save("plots/1.png");
}


void plot_pareto(long n) {
    vector<long> values = pareto_vec(n, 1, (long)std::trunc(sqrt(n)));
    long min = *std::min_element(values.begin(), values.end());
    long max = *std::max_element(values.begin(), values.end());
    cout << min << " " << max << " " << (long)std::trunc(sqrt(n)) << endl;

    vector<long> counts(max - min + 1, 0);

    for (long i = 0; i < n; ++i)
        try {
            counts[values[i] - min]++;
        } catch(std::out_of_range& e) {
            cout << e.what();
        }

    vector<double> x(counts.size()), y(counts.size());
    for (long i = 0; i < counts.size(); ++i) {
        x[i] = std::log(i), y[i] = std::log(counts[i]);
        // cout << i << (i < 10 ? "  | " : " | ") << counts[i] << endl;
    }
    plt::plot(counts);
    // plt::plot(x, y);
    plt::show();
}


void plot_pareto_old(long n) {
    vector<double> values = pareto_vec(n);  // , std::trunc(sqrt(n)));
    cout << *std::min_element(values.begin(), values.end()) << endl;

    vector<long> counts(1 / 0.01, 0);

    for (long i = 0; i < n; ++i)
        try {
            counts[std::trunc(values[i] * 100)]++;
        } catch(std::out_of_range& e) {
            cout << e.what();
        }

    vector<double> x(counts.size()), y(counts.size());
    for (long i = 0; i < counts.size(); ++i)
        x[i] = std::log(i), y[i] = std::log(counts[i]);
    // plt::plot(counts);
    plt::plot(x, y);
    plt::show();
}


void generate_GRG_degrees_plot(long n, long m, bool plot_degrees = true) {
    vector<long> w;
    vector<vector<pair<long, long>>> graphs = generate_m_GRG_edge_pairs_same_weights(n, m, n/2, w);

    long w_max = *std::max_element(w.begin(), w.end());
    int n_vertices = 10;
    vector<long> v(n_vertices);
    v[0] = std::trunc(uniform_distribution(generator) * w.size());
    for (int i = 1; i < n_vertices; ++i)
        while ((v[i] = std::trunc(uniform_distribution(generator) * w.size())) == v[i - 1]);
    
    vector<vector<double>> poisson(n_vertices);
    for (int i = 0; i < n_vertices; ++i)
        poisson[i] = generate_poisson_function(w[v[i]], 1, 2000, 1);

    plt::clf();
    std::sort(w.begin(), w.end());
    plt::plot(w);

    vector<vector<long>> degrees(m);
    vector<vector<long>> v_distribution(n_vertices);
    long i = 0;
    for (const auto& g: graphs) {
        degrees[i] = count_degrees(g);
        for (int j = 0; j < n_vertices; ++j)
            v_distribution[j].push_back(degrees[i][v[j]]);
        if (i < 3) {
            std::sort(degrees[i].begin(), degrees[i].end());
            plt::plot(degrees[i]);
            i++;
        }
    }

    // This is a plot of similarity of ALL degrees with weights
    // plt::legend({"Initial weights", "1st graph's degrees", "2nd graph's degrees", "3rd graph's degrees"}, "upper left");
    plt::save(std::string("../plots/PS3/PS3_4_size_") + std::to_string(n) + std::string(".png"));

    if (plot_degrees) {
        for (i = 0; i < n_vertices; ++i) {
            long max_d = *std::max_element(v_distribution[i].begin(), v_distribution[i].end());
            vector<double> counts(max_d + 1, 0);
            double sum = 0;

            for (auto d: v_distribution[i])
                counts[d]++;
            for (long j = 0; j < max_d + 1; ++j)
                counts[j] /= v_distribution[i].size();

            plt::clf();
            plt::plot(counts);
            plt::plot(poisson[i]);
            plt::xlim((long)1, max_d + 1);
            // plt::legend({"Vertex degree distribution", "Corresponding Poisson(w)"});
            plt::save(std::string("../plots/PS3/PS3_3_size_") + std::to_string(n) + std::string("_") + std::to_string(i+1) + std::string("_vertex.png"));
        }
    }

}


int main() {

    // test1();

    // test2();

    // ==========================
    // Solution for numerical assigment from Problem Set 2, Option 1:
    // construct_plot(1000, 6., 1000, 200, 1, "big_graphs");
    // ==========================

    // ==========================
    // Solution for numerical assigment from Problem Set 3:
    long n = 100;
    // auto grg = generate_GRG(n, get_pareto_generator(1, 50));

    /*
     *  for (long i = 0; i < n; ++i) {
     *  for (long j = 0; j < n; ++j)
     *  cout << grg[i][j] << " ";
     *  cout << endl;
     *  }
     */

    // plot_pareto(40000);
    generate_GRG_degrees_plot(500, 500);
    generate_GRG_degrees_plot(1000, 500);
    generate_GRG_degrees_plot(1500, 500);
    generate_GRG_degrees_plot(2000, 500);

//    generate_GRG_degrees_plot(5000, 5, false);
//    generate_GRG_degrees_plot(7500, 5, false);
//    generate_GRG_degrees_plot(10000, 5, false);
//    generate_GRG_degrees_plot(15000, 5, false);
    // ==========================

    return 0;
}
