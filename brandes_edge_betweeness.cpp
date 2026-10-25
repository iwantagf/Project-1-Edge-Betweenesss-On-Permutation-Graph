#ifndef BRANDES_EDGE_BETWEENESS_CPP
#define BRANDES_EDGE_BETWEENESS_CPP

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <utility>
#include <iomanip>

/**
 * Thuật toán Brandes tiêu chuẩn (2001) tính Độ trung gian cạnh (Edge Betweenness Centrality)
 * trên Đồ thị Hoán vị (Permutation Graph) G = (V, E).
 * Độ phức tạp: O(V * E), tương đương O(V^3) trên đồ thị hoán vị dày.
 */
class BrandesEdgeBetweenness {
public:
    int V;
    std::vector<int> pi;
    std::vector<std::pair<int, int>> edges; // Danh sách cạnh (v, u) với 1 <= v < u <= V và pi[v] > pi[u]
    std::vector<std::vector<std::pair<int, int>>> adj; // adj[u] chứa danh sách {v, edge_id}

    // Lưu trữ thông tin theo từng nguồn s để kiểm thử chuyên sâu (nếu cần)
    std::vector<std::vector<int>> dist_all;
    std::vector<std::vector<long long>> sigma_all;
    std::vector<std::vector<double>> delta_all;

    BrandesEdgeBetweenness(const std::vector<int>& permutation) {
        V = static_cast<int>(permutation.size());
        pi.assign(V + 1, 0);
        for (int i = 1; i <= V; ++i) {
            pi[i] = permutation[i - 1];
        }
        adj.resize(V + 1);
        build_graph();
    }

    // Dựng đồ thị hoán vị từ các cặp nghịch thế: (u - v) * (pi[u] - pi[v]) < 0
    void build_graph() {
        edges.clear();
        for (int i = 1; i <= V; ++i) {
            adj[i].clear();
        }
        for (int v = 1; v <= V; ++v) {
            for (int u = v + 1; u <= V; ++u) {
                if (pi[v] > pi[u]) {
                    int edge_id = static_cast<int>(edges.size());
                    edges.push_back({v, u});
                    adj[v].push_back({u, edge_id});
                    adj[u].push_back({v, edge_id});
                }
            }
        }
    }

    // Thực thi thuật toán Brandes tiêu chuẩn trong O(V * E)
    std::vector<double> compute_edge_betweenness(bool store_states = false) {
        int E = static_cast<int>(edges.size());
        std::vector<double> EB(E, 0.0);

        if (store_states) {
            dist_all.assign(V + 1, std::vector<int>(V + 1, -1));
            sigma_all.assign(V + 1, std::vector<long long>(V + 1, 0));
            delta_all.assign(V + 1, std::vector<double>(V + 1, 0.0));
        }

        for (int s = 1; s <= V; ++s) {
            std::stack<int> S;
            std::vector<std::vector<std::pair<int, int>>> P(V + 1); // P[w] lưu danh sách {v, edge_id}
            std::vector<long long> sigma(V + 1, 0);
            std::vector<int> d(V + 1, -1);
            std::vector<double> delta(V + 1, 0.0);

            sigma[s] = 1;
            d[s] = 0;
            std::queue<int> Q;
            Q.push(s);

            // Pha Xuôi (BFS): Tìm đường đi ngắn nhất và đếm số đường đi ngắn nhất sigma
            while (!Q.empty()) {
                int v = Q.front();
                Q.pop();
                S.push(v);

                for (const auto& edge : adj[v]) {
                    int w = edge.first;
                    int edge_id = edge.second;

                    if (d[w] < 0) {
                        Q.push(w);
                        d[w] = d[v] + 1;
                    }
                    if (d[w] == d[v] + 1) {
                        sigma[w] += sigma[v];
                        P[w].push_back({v, edge_id});
                    }
                }
            }

            // Pha Ngược: Tích lũy độ phụ thuộc delta và độ trung gian cạnh EB(e)
            while (!S.empty()) {
                int w = S.top();
                S.pop();
                for (const auto& pred : P[w]) {
                    int v = pred.first;
                    int edge_id = pred.second;
                    if (sigma[w] > 0) {
                        double c = (static_cast<double>(sigma[v]) / static_cast<double>(sigma[w])) * (1.0 + delta[w]);
                        EB[edge_id] += c;
                        delta[v] += c;
                    }
                }
            }

            if (store_states) {
                dist_all[s] = d;
                sigma_all[s] = sigma;
                delta_all[s] = delta;
            }
        }

        return EB;
    }
};

#ifdef STANDALONE_BRANDES
int main() {
    std::vector<int> perm = {4, 3, 1, 5, 2};
    BrandesEdgeBetweenness solver(perm);
    std::vector<double> eb = solver.compute_edge_betweenness(true);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== KẾT QUẢ BRANDES TIÊU CHUẨN ===\n";
    for (size_t i = 0; i < solver.edges.size(); ++i) {
        std::cout << "Cạnh (" << solver.edges[i].first << ", " << solver.edges[i].second << "): "
                  << "EB = " << eb[i] << "\n";
    }
    return 0;
}
#endif

#endif // BRANDES_EDGE_BETWEENESS_CPP