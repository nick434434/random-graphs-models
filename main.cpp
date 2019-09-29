#include<vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <random>
#include "matplotlibcpp.h"


namespace plt = matplotlibcpp;
using std::vector;
using std::find;
using std::map;
using std::pair;
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


std::random_device generator;
std::uniform_real_distribution<double> uniform_distribution(0, 1);
auto edge_decider = [](double prob) {
    assert(prob > 0 && prob < 1);
    return uniform_distribution(generator) < prob;
};


void print_graph(vector<vector<int>> G) {
    int N = G.size();
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            cout << G[i][j] << " ";
        cout << endl;
    }
}


vector<vector<int>> generate_ER(int n, double lambd) {
    vector<vector<int>> graph;
    double p = lambd / n;

    vector<int> line(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j)
            line[j] = (int)edge_decider(p);
        graph.push_back(line);
        line.assign(n, 0);
    }

    return graph;
}


long count_edges(vector<vector<int>> graph) {
    long total_num = 0;
    for (int i = 0; i < graph.size(); ++i)
        for (int j = i; j < graph.size(); ++j)
            total_num += graph[i][j];
    return total_num;
}


long count_triangles(vector<vector<int>> graph) {
    long total_num = 0;
    int n = graph.size();
    for (int i = 0; i < n; ++i) {
        vector<int> g_i = graph[i];
        for (int j = i + 1; j < n; ++j) {
            vector<int> g_j = graph[j];
            int g_i_j = g_i[j];
            for (int k = j + 1; k < n; ++k) {
                // cout << "(" << i << ", " << j << ", " << k << "): " << graph[i][j] * graph[j][k] * graph[k][i] << endl;
                total_num += g_i_j * g_j[k] * graph[i][k];
            }
        }
    }
    return total_num;
}


pair<vector<double>, vector<long>> generate_for_single_n(int N, double lambda_start, double lambda_step, int num_iters) {
    vector<long> results(num_iters, 0);
    vector<double> lambdas(num_iters, 0);
    double lmbd = lambda_start;

    for (int i = 0; i < num_iters; ++i) {
        lambdas[i] = lmbd;
        auto G = generate_ER(N, lmbd);
        results[i] = count_triangles(G);
        lmbd += lambda_step;
    }

    return {lambdas, results};
}


pair<vector<int>, vector<long>> generate_for_single_lambda(double lmbd, int N_start, int N_step, int num_iters) {
    vector<long> results(num_iters, 0);
    vector<int> Ns(num_iters, 0);
    int N = N_start;

    #pragma omp parallel for
    for (int i = 0; i < num_iters; ++i) {
        Ns[i] = N;
        auto G = generate_ER(N, lmbd);
        results[i] = count_triangles(G);
        N += N_step;
    }

    return {Ns, results};
}


vector<double> generate_poisson_function(double lmbd, long start, long end, long scale) {
    vector<double> res(end - start + 1);
    double fact = 1;
    for (long i = 2; i < start; ++i)
        fact *= i;
    for (long i = start; i <= end; ++i) {
        if (i > 0)
            fact *= i;
        res[i - start] = scale * std::pow(lmbd, i) / fact * exp(-lmbd);
        cout << res[i - start] << " ";
    }
    cout << endl;
    return res;
}


void construct_plot(int N_trials, double lmbd, int N_start, int N_step, int num_iters) {
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
        cout << "../plots/size_" + std::to_string(Ns[i]) + "_lambda_" + std::to_string((int)(lmbd * 10)/10) + "_trials_" +
                std::to_string(N_trials) + ".png" << endl;
        plt::save("../plots/size_" + std::to_string(Ns[i]) + "_lambda_" + std::to_string((int)(lmbd * 10)/10) + "_trials_" +
                  std::to_string(N_trials) + ".png");
        x.clear();
        y.clear();
        plt::clf();
    }

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
    plt::xlabel("lambda");
    plt::ylabel("# of triangles");
    plt::save("plots/1.png");
}


int main() {

    // test1();

    // test2();

    construct_plot(1000, 6., 600, 100, 2);

    return 0;
}