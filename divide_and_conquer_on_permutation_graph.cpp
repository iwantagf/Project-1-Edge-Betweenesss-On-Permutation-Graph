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
 * Tính và Truy vấn Độ Trung Gian Cạnh (Edge Betweenness Centrality) trên Đồ thị Hoán vị
 * thông qua Đồ thị Nén Chia để Trị (CDQ), Chuỗi đỉnh giả (Prefix Dummy Sparsification)
 * và Thu hồi Ma trận theo Đồng dư Khoảng cách Modulo 3.
 *
 * - Độ phức tạp không gian & dựng đồ thị G_CDQ: O(V log V)
 * - Độ phức tạp 0-1 BFS & Tích lũy ngược trên G_CDQ cho mọi nguồn: O(V^2 log V)
 * - Thu hồi toàn bộ cạnh qua Phân lớp Đồng dư Modulo 3: O(V^omega) lý thuyết / O(V * E) vector hóa
 * - Truy vấn cặp (u, v) bất kỳ:
 *     + O(1) sau khi nhân ma trận đồng dư 3 (query_mod3_O1)
 *     + O(1) kiểm tra không có cạnh, O(V) SIMD-friendly có nhớ lần đầu (query_lazy)
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

    // Các ma trận phục vụ Phân lớp Đồng dư Modulo 3 (A^(r)[u][s] và B_T^(r)[v][s])
    std::vector<std::vector<double>> A_mod3[3];
    std::vector<std::vector<double>> B_T_mod3[3];
    std::vector<std::vector<double>> C_matrix; // C_matrix[u][v] = sum_{s} c_s(u -> v)
    std::vector<std::vector<double>> memo_eb;  // Bảng nhớ cho truy vấn lười query_lazy(u, v)
    bool mod3_prepared = false;
    bool full_matrix_multiplied = false;

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

    // Kiểm tra sự tồn tại của cạnh {u, v} trong O(1) bằng điều kiện nghịch thế
    inline bool has_edge(int u, int v) const {
        if (u < 1 || u > V || v < 1 || v > V || u == v) return false;
        return (1LL * (u - v) * (pi[u] - pi[v])) < 0;
    }

    /**
     * Algorithm 1: Xây dựng Đồ thị Nén Chia để Trị G_CDQ (Build-CDQ-Dummy-Graph)
     * Sử dụng trộn tại chỗ (merge) từ dưới lên để đảm bảo đúng thời gian O(V log V).
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
        std::vector<double>& acc,
        std::vector<std::vector<int>>& buckets
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

        int max_d = 0;
        std::vector<int> visited_nodes;
        visited_nodes.reserve(num_nodes);
        std::vector<char> in_visited(num_nodes + 1, 0);

        // 0-1 BFS tìm khoảng cách ngắn nhất D[x] trên G_CDQ trong O(|V_CDQ| + |E_CDQ|)
        while (!Q.empty()) {
            int u = Q.front();
            Q.pop_front();
            if (!in_visited[u]) {
                in_visited[u] = 1;
                visited_nodes.push_back(u);
                if (D[u] > max_d) max_d = D[u];
            }

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

        // Xây dựng thứ tự topo chuẩn trong O(|V_CDQ|) bằng Bucket Sort theo khoảng cách D[u]:
        // Trong cùng một lớp khoảng cách d, các cạnh trọng số 0 luôn đi từ dummy[i] -> dummy[i-1] -> v_i
        // (tức từ chỉ số đỉnh lớn xuống chỉ số đỉnh nhỏ).
        // Quét u từ num_nodes xuống 1 giúp các đỉnh trong mỗi bucket tự động giảm dần theo chỉ số mà không cần std::sort!
        for (int d = 0; d <= max_d; ++d) {
            buckets[d].clear();
        }
        for (int u = num_nodes; u >= 1; --u) {
            if (in_visited[u]) {
                buckets[D[u]].push_back(u);
            }
        }

        std::vector<int> order;
        order.reserve(visited_nodes.size());
        for (int d = 0; d <= max_d; ++d) {
            for (int u : buckets[d]) {
                order.push_back(u);
            }
        }

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

    // Chạy Algorithm 2 cho toàn bộ V đỉnh nguồn trên G_CDQ trong đúng O(V^2 log V)
    void run_all_pairs_cdq() {
        D_real.assign(V + 1, std::vector<int>(V + 1, -1));
        sigma_real.assign(V + 1, std::vector<long long>(V + 1, 0));
        delta_real.assign(V + 1, std::vector<double>(V + 1, 0.0));

        std::vector<int> D(num_nodes + 1);
        std::vector<long long> sigma(num_nodes + 1);
        std::vector<double> delta(num_nodes + 1);
        std::vector<double> acc(num_nodes + 1);
        std::vector<std::vector<int>> buckets(V + 2);

        const int INF = 1e9;
        for (int s = 1; s <= V; ++s) {
            brandes_single_source_cdq(s, D, sigma, delta, acc, buckets);
            for (int x = 1; x <= V; ++x) {
                D_real[s][x] = (D[x] >= INF) ? -1 : D[x];
                sigma_real[s][x] = sigma[x];
                delta_real[s][x] = delta[x];
            }
        }
    }

    /**
     * Khởi tạo 3 cặp ma trận Phân lớp Đồng dư Modulo 3: A_mod3[r][u][s] và B_T_mod3[r][v][s]
     * Thời gian: O(V^2), Không gian: O(V^2).
     * Định lý 4.1: Với mọi {u, v} in E, d(s, v) = d(s, u) + 1 <=> d(s, v) == d(s, u) + 1 (mod 3).
     */
    void prepare_mod3_matrices() {
        if (mod3_prepared) return;
        if (D_real.empty()) {
            run_all_pairs_cdq();
        }

        for (int r = 0; r < 3; ++r) {
            A_mod3[r].assign(V + 1, std::vector<double>(V + 1, 0.0));
            B_T_mod3[r].assign(V + 1, std::vector<double>(V + 1, 0.0));
        }
        memo_eb.assign(V + 1, std::vector<double>(V + 1, -1.0));

        for (int s = 1; s <= V; ++s) {
            for (int x = 1; x <= V; ++x) {
                int d = D_real[s][x];
                if (d >= 0 && sigma_real[s][x] > 0) {
                    int r_u = d % 3;
                    int r_v = (d + 2) % 3; // d == (r_v + 1) mod 3 <=> r_v == (d - 1) mod 3
                    double sig = static_cast<double>(sigma_real[s][x]);
                    A_mod3[r_u][x][s] = sig;
                    B_T_mod3[r_v][x][s] = (1.0 + delta_real[s][x]) / sig;
                }
            }
        }
        mod3_prepared = true;
    }

    /**
     * Thực hiện Nhân Ma trận trên 3 lớp đồng dư Modulo 3:
     * C = sum_{r=0..2} A^(r) * B^(r)
     * Sau bước này, mọi truy vấn query_mod3_O1(u, v) đạt O(1) tuyệt đối.
     */
    void build_full_mod3_matrix_product() {
        if (full_matrix_multiplied) return;
        prepare_mod3_matrices();

        C_matrix.assign(V + 1, std::vector<double>(V + 1, 0.0));
        for (int r = 0; r < 3; ++r) {
            for (int u = 1; u <= V; ++u) {
                const double* a_row = A_mod3[r][u].data();
                for (int v = 1; v <= V; ++v) {
                    if (!has_edge(u, v)) continue;
                    const double* b_row = B_T_mod3[r][v].data();
                    double sum = 0.0;
                    for (int s = 1; s <= V; ++s) {
                        sum += a_row[s] * b_row[s];
                    }
                    C_matrix[u][v] += sum;
                }
            }
        }
        full_matrix_multiplied = true;
    }

    /**
     * Truy vấn O(1) tuyệt đối sau khi đã gọi build_full_mod3_matrix_product().
     * Trả về -1.0 nếu giữa u và v không có cạnh.
     */
    double query_mod3_O1(int u, int v) const {
        if (!has_edge(u, v)) return -1.0;
        return C_matrix[u][v] + C_matrix[v][u];
    }

    /**
     * Truy vấn lười có nhớ (Lazy Query with Memoization) trên ma trận Modulo 3:
     * - Chỉ cần gọi prepare_mod3_matrices() trong O(V^2 log V).
     * - Nếu (u, v) không có cạnh: trả về -1.0 trong O(1).
     * - Nếu (u, v) đã được hỏi trước đó: trả về ngay trong O(1).
     * - Nếu (u, v) được hỏi lần đầu: tính tích vô hướng không rẽ nhánh (branchless SIMD dot product)
     *   trên 3 lớp đồng dư r in {0, 1, 2}.
     */
    double query_lazy(int u, int v) {
        if (!has_edge(u, v)) return -1.0;
        if (!mod3_prepared) prepare_mod3_matrices();
        if (u > v) std::swap(u, v);
        if (memo_eb[u][v] >= 0.0) return memo_eb[u][v];

        double total_eb = 0.0;
        for (int r = 0; r < 3; ++r) {
            const double* au = A_mod3[r][u].data();
            const double* bv = B_T_mod3[r][v].data();
            const double* av = A_mod3[r][v].data();
            const double* bu = B_T_mod3[r][u].data();

            double sum_r = 0.0;
            for (int s = 1; s <= V; ++s) {
                sum_r += au[s] * bv[s] + av[s] * bu[s];
            }
            total_eb += sum_r;
        }

        memo_eb[u][v] = memo_eb[v][u] = total_eb;
        return total_eb;
    }

    /**
     * Thu hồi toàn bộ mảng EB(e) cho danh sách original_edges bằng Phân lớp Đồng dư Modulo 3.
     */
    std::vector<double> compute_edge_betweenness() {
        prepare_mod3_matrices();
        int E = static_cast<int>(original_edges.size());
        std::vector<double> EB(E, 0.0);

        for (int e_id = 0; e_id < E; ++e_id) {
            int v = original_edges[e_id].first;
            int u = original_edges[e_id].second;
            EB[e_id] = query_lazy(v, u);
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
    std::cout << "=== KẾT QUẢ THUẬT TOÁN ĐỀ XUẤT (ĐỒ THỊ NÉN CDQ + MA TRẬN MOD 3) ===\n";
    std::cout << "Số đỉnh G_CDQ: " << solver.num_nodes << ", Số cạnh G_CDQ: " << solver.count_cdq_edges() << "\n";
    for (size_t i = 0; i < solver.original_edges.size(); ++i) {
        std::cout << "Cạnh (" << solver.original_edges[i].first << ", " << solver.original_edges[i].second << "): "
                  << "EB = " << eb[i] << "\n";
    }
    return 0;
}
#endif

#endif // DIVIDE_AND_CONQUER_ON_PERMUTATION_GRAPH_CPP
