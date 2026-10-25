#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <random>
#include <algorithm>
#include <numeric>
#include <iomanip>

#include "brandes_edge_betweeness.cpp"
#include "divide_and_conquer_on_permutation_graph.cpp"

const double EPS = 1e-9;

// Hàm so sánh hai số thực
bool almost_equal(double a, double b, double eps = EPS) {
    return std::fabs(a - b) <= eps * std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

// Kiểm tra toàn diện một hoán vị: So sánh D, sigma, delta, EB và Định lý 2.1
void verify_permutation(const std::vector<int>& perm, const std::string& test_name) {
    int V = static_cast<int>(perm.size());
    BrandesEdgeBetweenness std_solver(perm);
    CDQPermutationGraphEBC cdq_solver(perm);

    std::vector<double> eb_std = std_solver.compute_edge_betweenness(true);
    std::vector<double> eb_cdq = cdq_solver.compute_edge_betweenness();

    // 1. Kiểm tra: Ràng buộc kích thước đồ thị nén G_CDQ
    if (V >= 2) {
        int ceil_log2_V = 0;
        while ((1 << ceil_log2_V) < V) ++ceil_log2_V;
        int max_allowed_nodes = V + V * ceil_log2_V;
        size_t max_allowed_edges = static_cast<size_t>(3 * V * ceil_log2_V);

        assert(cdq_solver.num_nodes <= max_allowed_nodes && "Vi phạm cận trên số đỉnh của Định lý 2.1!");
        assert(cdq_solver.count_cdq_edges() <= max_allowed_edges && "Vi phạm cận trên số cạnh của Định lý 2.1!");
    }

    // 2. Kiểm tra: Bảo toàn khoảng cách D(s, x), số đường đi ngắn nhất sigma và độ phụ thuộc delta
    for (int s = 1; s <= V; ++s) {
        for (int x = 1; x <= V; ++x) {
            assert(std_solver.dist_all[s][x] == cdq_solver.D_real[s][x] && "Sai lệch khoảng cách ngắn nhất D(s, x)!");
            assert(std_solver.sigma_all[s][x] == cdq_solver.sigma_real[s][x] && "Sai lệch số đường đi ngắn nhất sigma(s, x)!");
            assert(almost_equal(std_solver.delta_all[s][x], cdq_solver.delta_real[s][x]) && "Sai lệch độ phụ thuộc delta(s, x)!");
        }
    }

    // 3. Kiểm tra Độ trung gian cạnh EB(e)
    assert(eb_std.size() == eb_cdq.size() && "Số lượng cạnh không khớp!");
    for (size_t i = 0; i < eb_std.size(); ++i) {
        assert(std_solver.edges[i] == cdq_solver.original_edges[i]);
        assert(almost_equal(eb_std[i], eb_cdq[i]) && "Sai lệch giá trị Edge Betweenness Centrality EB(e)!");
    }

    std::cout << "[PASS] " << std::left << std::setw(42) << test_name
              << " | V = " << std::setw(3) << V
              << " | E = " << std::setw(5) << eb_std.size()
              << " | V_CDQ = " << std::setw(4) << cdq_solver.num_nodes
              << " | E_CDQ = " << std::setw(5) << cdq_solver.count_cdq_edges() << "\n";
}

// Test 1: Các trường hợp biên và đồ thị đặc biệt
void test_edge_cases() {
    std::cout << "=== 1. KIỂM THỬ TRƯỜNG HỢP BIÊN & ĐỒ THỊ ĐẶC BIỆT ===\n";
    // V = 1 (Đồ thị 1 đỉnh, 0 cạnh)
    verify_permutation({1}, "Single Vertex (V=1)");

    // Hoán vị đồng nhất (Đồ thị rỗng, 0 cạnh)
    verify_permutation({1, 2, 3, 4, 5, 6, 7, 8}, "Empty Graph / Identity Permutation (V=8)");

    // Hoán vị đảo ngược hoàn toàn (Đồ thị đầy đủ K_V, E = V(V-1)/2)
    verify_permutation({8, 7, 6, 5, 4, 3, 2, 1}, "Complete Graph K_8 / Reverse Perm (V=8)");

    // Đồ thị hình sao (Star Graph): Đỉnh cuối nối với tất cả các đỉnh trước đó
    verify_permutation({2, 3, 4, 5, 6, 7, 1}, "Star Graph (V=7)");

    // Đồ thị nhiều thành phần liên thông
    verify_permutation({3, 2, 1, 6, 5, 4, 9, 8, 7}, "Disconnected Components (3 x K_3, V=9)");
}

// Test 2: Kiểm tra vét cạn mọi hoán vị nhỏ (V = 1..6)
void test_exhaustive_small_permutations() {
    std::cout << "\n=== 2. KIỂM THỬ VÉT CẠN MỌI HOÁN VỊ NHỎ (V = 1..6) ===\n";
    int total_tested = 0;
    for (int V = 1; V <= 10; ++V) {
        std::vector<int> perm(V);
        std::iota(perm.begin(), perm.end(), 1);
        do {
            BrandesEdgeBetweenness std_solver(perm);
            CDQPermutationGraphEBC cdq_solver(perm);
            auto eb_std = std_solver.compute_edge_betweenness(true);
            auto eb_cdq = cdq_solver.compute_edge_betweenness();

            for (int s = 1; s <= V; ++s) {
                for (int x = 1; x <= V; ++x) {
                    assert(std_solver.dist_all[s][x] == cdq_solver.D_real[s][x]);
                    assert(std_solver.sigma_all[s][x] == cdq_solver.sigma_real[s][x]);
                    assert(almost_equal(std_solver.delta_all[s][x], cdq_solver.delta_real[s][x]));
                }
            }
            for (size_t i = 0; i < eb_std.size(); ++i) {
                assert(almost_equal(eb_std[i], eb_cdq[i]));
            }
            ++total_tested;
        } while (std::next_permutation(perm.begin(), perm.end()));
    }
    std::cout << "[PASS] Đã kiểm tra chính xác tuyệt đối toàn bộ " << total_tested << " hoán vị với V in [1..10].\n";
}

// Test 3: Kiểm thử ngẫu nhiên với kích thước trung bình và lớn
void test_random_permutations() {
    std::cout << "\n=== 3. KIỂM THỬ HOÁN VỊ NGẪU NHIÊN (RANDOM STRESS TESTS) ===\n";
    std::mt19937 rng(20261024);
    std::vector<int> sizes = {15, 30, 50, 80, 120, 300, 400, 500};

    for (int V : sizes) {
        std::vector<int> perm(V);
        std::iota(perm.begin(), perm.end(), 1);
        std::shuffle(perm.begin(), perm.end(), rng);
        verify_permutation(perm, "Random Permutation (V=" + std::to_string(V) + ")");
    }
}

int main() {
    std::cout << "================================================================================\n";
    std::cout << "   BỘ KIỂM THỬ ĐƠN VỊ (UNITTEST) - EDGE BETWEENNESS TRÊN ĐỒ THỊ HOÁN VỊ\n";
    std::cout << "================================================================================\n\n";

    test_edge_cases();
    test_exhaustive_small_permutations();
    test_random_permutations();

    std::cout << "\n================================================================================\n";
    std::cout << "   TẤT CẢ CÁC BÀI KIỂM THỬ (UNITTESTS) ĐỀU THÀNH CÔNG TUYỆT ĐỐI (100% PASSED)!\n";
    std::cout << "================================================================================\n";
    return 0;
}