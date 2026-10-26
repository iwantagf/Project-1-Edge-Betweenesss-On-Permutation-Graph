#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <random>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <chrono>

#include "brandes_edge_betweeness.cpp"
#include "divide_and_conquer_on_permutation_graph.cpp"

const double EPS = 1e-6;

bool almost_equal(double a, double b, double eps = EPS) {
    if (a < 0 && b < 0) return true;
    return std::fabs(a - b) <= eps * std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

// Kiểm tra toàn diện hoán vị: Vét cạn đủ V^2 cặp (u, v) và Q lần truy vấn ngẫu nhiên
void verify_permutation_and_queries(const std::vector<int>& perm, int Q_queries, const std::string& test_name) {
    int V = static_cast<int>(perm.size());
    BrandesEdgeBetweenness std_solver(perm);
    CDQPermutationGraphEBC cdq_solver(perm);

    // Chạy Brandes tiêu chuẩn
    std_solver.compute_edge_betweenness(true);

    // Chuẩn bị ma trận Modulo 3 cho thuật toán đề xuất
    cdq_solver.prepare_mod3_matrices();
    cdq_solver.build_full_mod3_matrix_product();

    // 1. Kiểm tra Định lý 2.1: Ràng buộc kích thước đồ thị nén G_CDQ
    if (V >= 2) {
        int ceil_log2_V = 0;
        while ((1 << ceil_log2_V) < V) ++ceil_log2_V;
        int max_allowed_nodes = V + V * ceil_log2_V;
        size_t max_allowed_edges = static_cast<size_t>(3 * V * ceil_log2_V);

        assert(cdq_solver.num_nodes <= max_allowed_nodes && "Vi phạm cận trên số đỉnh của Định lý 2.1!");
        assert(cdq_solver.count_cdq_edges() <= max_allowed_edges && "Vi phạm cận trên số cạnh của Định lý 2.1!");
    }

    // 2. Kiểm tra Bổ đề 3.2 & Định lý 3.3: Bảo toàn khoảng cách D(s, x), số đường đi ngắn nhất sigma và độ phụ thuộc delta
    for (int s = 1; s <= V; ++s) {
        for (int x = 1; x <= V; ++x) {
            assert(std_solver.dist_all[s][x] == cdq_solver.D_real[s][x] && "Sai lệch khoảng cách ngắn nhất D(s, x)!");
            assert(std_solver.sigma_all[s][x] == cdq_solver.sigma_real[s][x] && "Sai lệch số đường đi ngắn nhất sigma(s, x)!");
            assert(almost_equal(std_solver.delta_all[s][x], cdq_solver.delta_real[s][x], 1e-6) && "Sai lệch độ phụ thuộc delta(s, x)!");
        }
    }

    // 3. Kiểm tra vét cạn ĐỦ V^2 cặp (u, v) giữa Brandes và CDQ (query_lazy & query_mod3_O1)
    for (int u = 1; u <= V; ++u) {
        for (int v = 1; v <= V; ++v) {
            double q_std = std_solver.query_edge(u, v);
            double q_lazy = cdq_solver.query_lazy(u, v);
            double q_mod3 = cdq_solver.query_mod3_O1(u, v);

            if (!almost_equal(q_std, q_lazy, 1e-6)) {
                std::cerr << "[ERROR] V^2 Exhaustive Mismatch at u=" << u << ", v=" << v
                          << " | Brandes: " << q_std << " vs CDQ Lazy: " << q_lazy << "\n";
                assert(false);
            }
            if (!almost_equal(q_std, q_mod3, 1e-6)) {
                std::cerr << "[ERROR] V^2 Exhaustive Mismatch at u=" << u << ", v=" << v
                          << " | Brandes: " << q_std << " vs CDQ Mod3 O(1): " << q_mod3 << "\n";
                assert(false);
            }
        }
    }

    // 4. Kiểm tra Q lần truy vấn ngẫu nhiên từng cặp (u, v)
    std::mt19937 rng(1337 + V);
    std::uniform_int_distribution<int> dist(1, V);
    for (int q = 0; q < Q_queries; ++q) {
        int u = dist(rng);
        int v = dist(rng);
        double q_std = std_solver.query_edge(u, v);
        double q_res = cdq_solver.query_mod3_O1(u, v);
        if (!almost_equal(q_std, q_res, 1e-6)) {
            std::cerr << "[ERROR] Random Query Mismatch at u=" << u << ", v=" << v
                      << " | Brandes: " << q_std << " vs CDQ: " << q_res << "\n";
            assert(false);
        }
    }

    std::cout << "[PASS] " << std::left << std::setw(42) << test_name
              << " | V = " << std::setw(3) << V
              << " | E = " << std::setw(5) << std_solver.count_edges()
              << " | V^2 pairs checked"
              << " | Random Queries Q = " << std::setw(6) << Q_queries
              << " | V_CDQ = " << std::setw(4) << cdq_solver.num_nodes
              << " | E_CDQ = " << std::setw(5) << cdq_solver.count_cdq_edges() << "\n";
}

// Test 1: Các trường hợp biên và đồ thị đặc biệt
void test_edge_cases() {
    std::cout << "=== 1. KIỂM THỬ TRƯỜNG HỢP BIÊN & ĐỒ THỊ ĐẶC BIỆT ===\n";
    verify_permutation_and_queries({1}, 100, "Single Vertex (V=1)");
    verify_permutation_and_queries({1, 2, 3, 4, 5, 6, 7, 8}, 1000, "Empty Graph / Identity Permutation (V=8)");
    verify_permutation_and_queries({8, 7, 6, 5, 4, 3, 2, 1}, 1000, "Complete Graph K_8 / Reverse Perm (V=8)");
    verify_permutation_and_queries({2, 3, 4, 5, 6, 7, 1}, 1000, "Star Graph (V=7)");
    verify_permutation_and_queries({3, 2, 1, 6, 5, 4, 9, 8, 7}, 1000, "Disconnected Components (3 x K_3, V=9)");
}

// Test 2: Kiểm tra vét cạn mọi hoán vị nhỏ (V = 1..6)
void test_exhaustive_small_permutations() {
    std::cout << "\n=== 2. KIỂM THỬ VÉT CẠN MỌI HOÁN VỊ NHỎ (V = 1..10) ===\n";
    int total_tested = 0;
    for (int V = 1; V <= 10; ++V) {
        std::vector<int> perm(V);
        std::iota(perm.begin(), perm.end(), 1);
        do {
            verify_permutation_and_queries(perm, 50, "Perm (V=" + std::to_string(V) + ")");
            ++total_tested;
        } while (std::next_permutation(perm.begin(), perm.end()));
    }
    std::cout << "[PASS] Đã kiểm tra chính xác tuyệt đối toàn bộ " << total_tested << " hoán vị với V in [1..10] (V^2 pairs + Q queries).\n";
}

// Test 3: Kiểm thử ngẫu nhiên với kích thước (V = 10, 25, 50, 100) và Q = 5000 truy vấn ngẫu nhiên
void test_random_permutations() {
    std::cout << "\n=== 3. KIỂM THỬ HOÁN VỊ NGẪU NHIÊN VÀ TRUY VẤN Q LẦN ===\n";
    std::mt19937 rng(20261026);
    std::vector<int> sizes = {10, 25, 50, 100, 200, 400, 500, 600, 700, 800, 990, 1000};

    for (int V : sizes) {
        std::vector<int> perm(V);
        std::iota(perm.begin(), perm.end(), 1);
        std::shuffle(perm.begin(), perm.end(), rng);
        verify_permutation_and_queries(perm, 100000, "Random Permutation (V=" + std::to_string(V) + ")");
    }

    std::cout << "\n=== 4. ĐỒ THỊ ĐẦY ĐỦ THEO TRONG V BÊN TRÊN ===\n";

    for (int V : sizes) {
        std::vector<int> perm(V);
        std::iota(perm.rbegin(), perm.rend(), 1);
        verify_permutation_and_queries(perm, 100000, "Complete Graph (V=" + std::to_string(V) + ")");
    }
}

int main() {
    std::cout << "================================================================================\n";
    std::cout << "   BỘ KIỂM THỬ ĐƠN VỊ (UNITTEST) - V^2 PAIRS & Q RANDOM QUERIES\n";
    std::cout << "================================================================================\n\n";

    test_edge_cases();
    test_exhaustive_small_permutations();
    test_random_permutations();

    std::cout << "\n================================================================================\n";
    std::cout << "   TẤT CẢ CÁC BÀI KIỂM THỬ (UNITTESTS) ĐỀU THÀNH CÔNG TUYỆT ĐỐI (100% PASSED)!\n";
    std::cout << "================================================================================\n";
    return 0;
}
