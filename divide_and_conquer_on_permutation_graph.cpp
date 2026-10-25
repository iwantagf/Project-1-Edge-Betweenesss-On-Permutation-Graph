#ifndef DIVIDE_AND_CONQUER_ON_PERMUTATION_GRAPH_CPP
#define DIVIDE_AND_CONQUER_ON_PERMUTATION_GRAPH_CPP

#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>
#include <iomanip>

/**
 * Thuật toán Đề xuất:
 * Tính Độ Trung Gian Cạnh (Edge Betweenness Centrality) trên Đồ thị Hoán vị
 * thông qua Đồ thị Nén Chia để Trị (CDQ) kết hợp Chuỗi đỉnh giả (Prefix Dummy Sparsification).
 *
 * - Độ phức tạp không gian và dựng đồ thị G_CDQ: O(V log V)
 * - Độ phức tạp chạy 0-1 BFS và tích lũy ngược trên G_CDQ: O(V^2 log V)
 * - Độ phức tạp hợp nhất toàn cục: O(min(V * E, V^2 log V + E)) = O(E^1.25 * sqrt(log E))
 */
class CDQPermutationGraphEBC {
public:
    struct DirectedEdge {
        int to;
        int weight; // Trọng số thuộc {0, 1}
    };

    int V;
    int num_nodes; // Tổng số đỉnh trong G_CDQ (1..V là đỉnh thật, > V là đỉnh giả)
    std::vector<int> pi;
    std::vector<std::pair<int, int>> original_edges; // Các cạnh gốc (v, u) với v < u và pi[v] > pi[u]
    std::vector<std::vector<DirectedEdge>> adj_cdq;

    // Bảng kết quả toàn cục từ G_CDQ cho mọi nguồn s in [1..V] và đỉnh thật x in [1..V]
    std::vector<std::vector<int>> D_real;
    std::vector<std::vector<long long>> sigma_real;
    std::vector<std::vector<double>> delta_real;

    CDQPermutationGraphEBC(const std::vector<int>& permutation) {
        V = static_cast<int>(permutation.size());
        pi.assign(V + 1, 0);
        for (int i = 1; i <= V; ++i) {
            pi[i] = permutation[i - 1];
        }
        num_nodes = V;
        adj_cdq.resize(V + 1);

        // Liệt kê danh sách cạnh gốc E của đồ thị hoán vị
        for (int v = 1; v <= V; ++v) {
            for (int u = v + 1; u <= V; ++u) {
                if (pi[v] > pi[u]) {
                    original_edges.push_back({v, u});
                }
            }
        }

        // Xây dựng Đồ thị Nén Chia để Trị G_CDQ theo Algorithm 1
        if (V >= 1) {
            cdq_build(1, V);
        }
    }

    int new_dummy_node() {
        ++num_nodes;
        adj_cdq.emplace_back();
        return num_nodes;
    }

    void add_edge(int u, int v, int w) {
        adj_cdq[u].push_back({v, w});
    }

    /**
     * Algorithm 1: Xây dựng Đồ thị Nén Chia để Trị G_CDQ (Build-CDQ-Dummy-Graph)
     */
    void cdq_build(int l, int r) {
        if (l >= r) return;
        int mid = l + (r - l) / 2;
        cdq_build(l, mid);
        cdq_build(mid + 1, r);

        // 1. Chiều từ phải [mid + 1, r] sang trái [l, mid] (Sort giảm dần theo pi)
        std::vector<int> L(mid - l + 1);
        std::iota(L.begin(), L.end(), l);
        std::sort(L.begin(), L.end(), [&](int a, int b) {
            return pi[a] > pi[b];
        });

        std::vector<int> R(r - mid);
        std::iota(R.begin(), R.end(), mid + 1);
        std::sort(R.begin(), R.end(), [&](int a, int b) {
            return pi[a] > pi[b];
        });

        std::vector<int> dummy_L(L.size());
        for (size_t i = 0; i < L.size(); ++i) {
            dummy_L[i] = new_dummy_node();
            add_edge(dummy_L[i], L[i], 0);
            if (i > 0) {
                add_edge(dummy_L[i], dummy_L[i - 1], 0);
            }
        }

        size_t ptr = 0;
        for (int u : R) {
            while (ptr < L.size() && pi[L[ptr]] > pi[u]) {
                ++ptr;
            }
            if (ptr >= 1) {
                add_edge(u, dummy_L[ptr - 1], 1);
            }
        }

        // 2. Chiều từ trái [l, mid] sang phải [mid + 1, r] (Sort tăng dần theo pi)
        std::vector<int> L_prime(L.rbegin(), L.rend());
        std::vector<int> R_prime(R.rbegin(), R.rend());

        std::vector<int> dummy_R(R_prime.size());
        for (size_t j = 0; j < R_prime.size(); ++j) {
            dummy_R[j] = new_dummy_node();
            add_edge(dummy_R[j], R_prime[j], 0);
            if (j > 0) {
                add_edge(dummy_R[j], dummy_R[j - 1], 0);
            }
        }

        ptr = 0;
        for (int v : L_prime) {
            while (ptr < R_prime.size() && pi[R_prime[ptr]] < pi[v]) {
                ++ptr;
            }
            if (ptr >= 1) {
                add_edge(v, dummy_R[ptr - 1], 1);
            }
        }
    }

    /**
     * Algorithm 2: Pha Xuôi 0-1 BFS và Pha Ngược Brandes trên G_CDQ từ một đỉnh nguồn s
     */
    void brandes_single_source_cdq(
        int s,
        std::vector<int>& D,
        std::vector<long long>& sigma,
        std::vector<double>& delta,
        std::vector<double>& acc
    ) {
        const int INF = 1e9;
        std::fill(D.begin(), D.end(), INF);
        std::fill(sigma.begin(), sigma.end(), 0);
        std::fill(delta.begin(), delta.end(), 0.0);
        std::fill(acc.begin(), acc.end(), 0.0);

        D[s] = 0;
        sigma[s] = 1;
        std::deque<int> Q;
        Q.push_back(s);

        std::vector<int> raw_order;
        raw_order.reserve(num_nodes * 2);

        // 0-1 BFS tìm khoảng cách ngắn nhất D[x] trên G_CDQ
        while (!Q.empty()) {
            int u = Q.front();
            Q.pop_front();
            raw_order.push_back(u);

            for (const auto& edge : adj_cdq[u]) {
                int y = edge.to;
                int w = edge.weight;
                if (D[u] + w < D[y]) {
                    D[y] = D[u] + w;
                    if (w == 0) {
                        Q.push_front(y);
                    } else {
                        Q.push_back(y);
                    }
                }
            }
        }

        // Chuẩn hóa thứ tự topo trên DAG cạnh chặt:
        // Trong cùng một lớp khoảng cách d: dummy[i] -> dummy[i-1] -> v_i
        // Do các đỉnh giả có chỉ số dummy[i] > dummy[i-1] > V >= v_i, sắp xếp giảm dần theo chỉ số đỉnh
        // trong cùng lớp khoảng cách d đảm bảo thứ tự topo tuyệt đối cho mọi cạnh trọng số 0.
        std::vector<int> order;
        order.reserve(num_nodes);
        std::vector<char> visited(num_nodes + 1, 0);
        for (int u : raw_order) {
            if (!visited[u]) {
                visited[u] = 1;
                order.push_back(u);
            }
        }
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            if (D[a] != D[b]) return D[a] < D[b];
            return a > b;
        });

        // Pha Xuôi: Cộng dồn sigma theo thứ tự topo của các cạnh chặt
        for (int u : order) {
            for (const auto& edge : adj_cdq[u]) {
                int y = edge.to;
                int w = edge.weight;
                if (D[y] == D[u] + w) {
                    sigma[y] += sigma[u];
                }
            }
        }

        // Pha Ngược: Lan truyền ngược độ phụ thuộc qua chuỗi đỉnh giả
        for (int idx = static_cast<int>(order.size()) - 1; idx >= 0; --idx) {
            int u = order[idx];
            if (u <= V) { // u là đỉnh thật thuộc V
                for (const auto& edge : adj_cdq[u]) {
                    int dummy = edge.to;
                    if (D[dummy] == D[u] + 1) {
                        delta[u] += static_cast<double>(sigma[u]) * acc[dummy];
                    }
                }
            } else { // u là đỉnh giả dummy
                for (const auto& edge : adj_cdq[u]) {
                    int y = edge.to;
                    if (D[y] == D[u]) { // Chỉ xét cạnh trọng số 0 chặt
                        if (y <= V) {
                            if (sigma[y] > 0) {
                                acc[u] += (1.0 + delta[y]) / static_cast<double>(sigma[y]);
                            }
                        } else {
                            acc[u] += acc[y];
                        }
                    }
                }
            }
        }
    }

    // Chạy Algorithm 2 cho toàn bộ V đỉnh nguồn trên G_CDQ trong O(V^2 log V)
    void run_all_pairs_cdq() {
        D_real.assign(V + 1, std::vector<int>(V + 1, -1));
        sigma_real.assign(V + 1, std::vector<long long>(V + 1, 0));
        delta_real.assign(V + 1, std::vector<double>(V + 1, 0.0));

        std::vector<int> D(num_nodes + 1);
        std::vector<long long> sigma(num_nodes + 1);
        std::vector<double> delta(num_nodes + 1);
        std::vector<double> acc(num_nodes + 1);

        const int INF = 1e9;
        for (int s = 1; s <= V; ++s) {
            brandes_single_source_cdq(s, D, sigma, delta, acc);
            for (int x = 1; x <= V; ++x) {
                D_real[s][x] = (D[x] >= INF) ? -1 : D[x];
                sigma_real[s][x] = sigma[x];
                delta_real[s][x] = delta[x];
            }
        }
    }


    std::vector<double> compute_edge_betweenness() {
        run_all_pairs_cdq();
        int E = static_cast<int>(original_edges.size());
        std::vector<double> EB(E, 0.0);

        for (int e_id = 0; e_id < E; ++e_id) {
            int v = original_edges[e_id].first;
            int u = original_edges[e_id].second;
            double total_eb = 0.0;

            for (int s = 1; s <= V; ++s) {
                int dv = D_real[s][v];
                int du = D_real[s][u];
                if (dv < 0 || du < 0) continue;

                if (du == dv + 1 && sigma_real[s][u] > 0) {
                    total_eb += (static_cast<double>(sigma_real[s][v]) / static_cast<double>(sigma_real[s][u])) * (1.0 + delta_real[s][u]);
                } else if (dv == du + 1 && sigma_real[s][v] > 0) {
                    total_eb += (static_cast<double>(sigma_real[s][u]) / static_cast<double>(sigma_real[s][v])) * (1.0 + delta_real[s][v]);
                }
            }
            EB[e_id] = total_eb;
        }
        return EB;
    }

    // Tổng số cạnh trong đồ thị nén G_CDQ
    size_t count_cdq_edges() const {
        size_t total = 0;
        for (const auto& list : adj_cdq) {
            total += list.size();
        }
        return total;
    }
};

#ifdef STANDALONE_CDQ
int main() {
    std::vector<int> perm = {4, 3, 1, 5, 2};
    CDQPermutationGraphEBC solver(perm);
    std::vector<double> eb = solver.compute_edge_betweenness();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== KẾT QUẢ THUẬT TOÁN ĐỀ XUẤT (ĐỒ THỊ NÉN CDQ) ===\n";
    std::cout << "Số đỉnh G_CDQ: " << solver.num_nodes << ", Số cạnh G_CDQ: " << solver.count_cdq_edges() << "\n";
    for (size_t i = 0; i < solver.original_edges.size(); ++i) {
        std::cout << "Cạnh (" << solver.original_edges[i].first << ", " << solver.original_edges[i].second << "): "
                  << "EB = " << eb[i] << "\n";
    }
    return 0;
}
#endif

#endif // DIVIDE_AND_CONQUER_ON_PERMUTATION_GRAPH_CPP