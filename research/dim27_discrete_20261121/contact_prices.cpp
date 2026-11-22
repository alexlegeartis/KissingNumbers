// Price a fresh midpoint of two deleted Leech roots after optimizing all 32 switches.
#include "switch_bounds.hpp"
#include <algorithm>
#include <array>
#include <filesystem>
#include <iostream>
#include <set>
#ifdef __AVX2__
#include <immintrin.h>
#endif

using Profile = std::array<int, 6>;
namespace fs = std::filesystem;
constexpr int n = 13248;

template<class T>
std::vector<T> read(const fs::path& path, size_t count) {
    std::vector<T> rows(count);
    std::ifstream in(path, std::ios::binary);
    in.read(reinterpret_cast<char*>(rows.data()), count * sizeof(T));
    if (!in) throw std::runtime_error("Incomplete input: " + path.string());
    return rows;
}

std::vector<Profile> minimal_profiles(const fs::path& path, const SwitchGroups& groups) {
    auto data = read<int32_t>(path, 1 + 6 * 85);
    if (data.front() != 85) throw std::runtime_error("Unexpected profile count");
    std::set<Profile> distinct;
    for (int i = 0; i < 85; ++i) {
        Profile p{};
        for (int axis = 0; axis < 3; ++axis) {
            auto colors = groups.colors[axis];
            int a = data[1 + 6 * i + colors[0]], b = data[1 + 6 * i + colors[1]];
            p[2 * axis] = std::min(a, b);
            p[2 * axis + 1] = std::max(a, b);
        }
        distinct.insert(p);
    }
    std::vector<Profile> profiles;
    for (auto p : distinct) {
        bool dominated = false;
        for (auto q : distinct) {
            if (q == p) continue;
            bool below = true;
            for (int k = 0; k < 6; ++k) below &= q[k] <= p[k];
            dominated |= below;
        }
        if (!dominated) profiles.push_back(p);
    }
    return profiles;
}

int price(const int8_t* a, const int8_t* b, const SwitchGroups& groups,
          const std::vector<Profile>& profiles) {
    int cost[64][5]{};
    auto count = [&](int m) {
        const int bound = 8 - a[m] - b[m];
        if (bound < 4) ++cost[groups.slot[m]][std::max(0, bound + 1)];
    };
#ifdef __AVX2__
    static_assert(n % 32 == 0);
    for (int k = 0; k < n; k += 32) {
        auto first = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + k));
        auto second = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b + k));
        auto eligible = _mm256_cmpgt_epi8(_mm256_add_epi8(first, second), _mm256_set1_epi8(4));
        unsigned mask = static_cast<unsigned>(_mm256_movemask_epi8(eligible));
        while (mask) {
            count(k + __builtin_ctz(mask));
            mask &= mask - 1;
        }
    }
#else
    for (int m = 0; m < n; ++m) count(m);
#endif
    for (auto& row : cost)
        for (int i = 1; i < 5; ++i) row[i] += row[i - 1];
    int plane_cost[3][5][5]{};
    for (int g = 0; g < 32; ++g) {
        auto& total = plane_cost[groups.plane[g]];
        auto& first = cost[2 * g];
        auto& second = cost[2 * g + 1];
        for (int low = 0; low < 5; ++low)
            for (int high = low; high < 5; ++high)
                total[low][high] += std::min(first[low] + second[high], first[high] + second[low]);
    }
    int best = n;
    for (auto p : profiles) {
        int sum = 0;
        for (int axis = 0; axis < 3; ++axis)
            sum += plane_cost[axis][p[2 * axis]][p[2 * axis + 1]];
        best = std::min(best, sum);
    }
    return best;
}

int main(int argc, char** argv) {
    if (argc != 4) return 1;
    fs::path faces = argv[1], contacts = argv[2];
    auto products = read<int8_t>(faces / "products", size_t(n) * n);
    auto offsets = read<int32_t>(faces / "offsets", n + 1);
    auto edges = read<int32_t>(faces / "edges", offsets.back());
    auto incident = read<int32_t>(faces / "incident", 2 * n);
    SwitchGroups groups;
    groups.load(argv[3], read<int32_t>(faces / "colors", n));
    auto profiles = minimal_profiles(contacts / "profiles", groups);
    int min_price = n;
    int64_t fresh = 0, by_price[8]{};
    #pragma omp parallel for schedule(dynamic,1) reduction(+:fresh,by_price[:8]) reduction(min:min_price)
    for (int a = 0; a < n; ++a) {
        for (int j = offsets[a]; j < offsets[a + 1]; ++j) {
            int b = edges[j];
            bool original = false;
            for (int i = 0; i < 2; ++i) for (int k = 0; k < 2; ++k)
                original |= incident[2*a+i] == incident[2*b+k];
            if (original) continue;
            int cuts = price(products.data() + size_t(a) * n, products.data() + size_t(b) * n, groups, profiles);
            ++fresh;
            ++by_price[std::min(7, cuts)];
            min_price = std::min(min_price, cuts);
        }
    }
    std::cout << "{\"fresh_contacts\":" << fresh << ",\"sorted_minimal_profiles\":" << profiles.size()
              << ",\"minimum_conflicting_midpoint_pairs\":" << min_price << ",\"attaining_minimum\":"
              << (min_price < 7 ? by_price[min_price] : -1) << "}\n";
}
