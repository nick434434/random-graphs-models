//
// Created by nick434434 on 14-10-19.
//
#include "prob_stuff.h"
#include <vector>
#include <cmath>
#include <iostream>

using std::vector;
using std::pair;
using std::cout;
using std::endl;


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
