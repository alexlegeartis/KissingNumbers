// Exact profile feasibility for every midpoint of a Leech contact.
#include "switch_bounds.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace fs = std::filesystem;
using Profile = std::array<int32_t, 6>;
constexpr std::array<int, 6> powers = {1, 5, 25, 125, 625, 3125};
constexpr int profiles_count = 15625;

template<class T>
std::vector<T> read_array(const fs::path& path, size_t size) {
    std::vector<T> result(size);
    std::ifstream input(path, std::ios::binary);
    input.read(reinterpret_cast<char*>(result.data()), sizeof(T) * size);
    if (!input) throw std::runtime_error("Incomplete input: " + path.string());
    return result;
}

std::array<bool, profiles_count> feasibility(const fs::path& path) {
    int count = read_array<int32_t>(path, 1).front();
    if (count != 85) throw std::runtime_error("Unexpected minimal-profile count");
    auto table = read_array<int32_t>(path, 1 + 6 * count);
    std::array<bool, profiles_count> allowed{};
    for (int key = 0; key < profiles_count; ++key) {
        for (int i = 0; i < count && !allowed[key]; ++i) {
            bool covers = true;
            for (int c = 0; c < 6; ++c)
                covers &= key / powers[c] % 5 >= table[1 + 6 * i + c];
            allowed[key] = covers;
        }
    }
    return allowed;
}

int main(int argc, char** argv) {
    if (argc != 2 && argc != 3) return 1;
    const fs::path dir = argv[1];
    auto header = read_array<int32_t>(dir / "capmeta", 9);
    const int roots = header[0], midpoints = header[1];
    if (roots != 196560 || midpoints != 13248 || header[2] != 0 ||
        header[8] != midpoints || !std::is_sorted(header.begin() + 2, header.end()))
        throw std::runtime_error("Invalid matrix dimensions");
    auto products = read_array<int8_t>(dir / "cap.bin", size_t(roots) * midpoints);
    auto offset = read_array<int64_t>(dir / "live_offsets", roots + 1);
    if (offset.front() != 0 || offset.back() != 452088000 || !std::is_sorted(offset.begin(), offset.end()))
        throw std::runtime_error("Invalid Leech contact offsets");
    auto neighbors = read_array<int32_t>(dir / "live_neighbors", offset.back());
    auto allowed = feasibility(dir / "profiles");
    SwitchGroups groups;
    const bool switched = argc == 3;
    if (switched) {
        std::vector<int32_t> labels(midpoints);
        for (int c = 0; c < 6; ++c)
            std::fill(labels.begin() + header[2+c], labels.begin() + header[3+c], c);
        groups.load(argv[2], labels);
        groups.fold(allowed, 5);
    }
    int max_degree = 0;
    for (int r = 0; r < roots; ++r)
        max_degree = std::max(max_degree, int(offset[r + 1] - offset[r]));

    int64_t checked = 0, surviving = 0;
    constexpr int block = 32;
    #pragma omp parallel reduction(+:checked,surviving)
    {
        std::vector<int8_t> local(size_t(block) * midpoints);
        std::vector<int> live(max_degree), testing(max_degree), upper(max_degree);
        std::vector<std::array<int8_t, 6>> lower(max_degree);
        std::array<std::array<std::vector<int>, 8>, 6> grouped;
        #pragma omp for schedule(dynamic,1)
        for (int start = 0; start < roots; start += block) {
            int size = std::min(block, roots - start);
            for (int m = 0; m < midpoints; ++m)
                for (int i = 0; i < size; ++i)
                    local[size_t(i) * midpoints + m] = products[size_t(m) * roots + start + i];
            for (int i = 0; i < size; ++i) {
                int r = start + i;
                const auto* scores = local.data() + size_t(i) * midpoints;
                for (int c = 0; c < 6; ++c) {
                    for (auto& v : grouped[c]) v.clear();
                    for (int m = header[2 + c]; m < header[3 + c]; ++m)
                        if (scores[m] >= -1 && (!switched || groups.largest[m]))
                            grouped[c][scores[m] + 1].push_back(m);
                }
                const auto* partners = neighbors.data() + offset[r];
                int remaining = int(offset[r + 1] - offset[r]);
                checked += remaining;
                std::iota(live.begin(), live.begin() + remaining, 0);
                std::fill_n(upper.begin(), remaining, profiles_count - 1);
                for (int j = 0; j < remaining; ++j) lower[j].fill(-1);

                // A level l requires scores[m] + products[m,s] <= 8-l for the whole color.
                // Levels are visited in a pruning order. A failed level lowers a base-5 digit.
                for (int level : {3, 2, 4, 1, 0}) {
                    for (int c = 0; c < 6 && remaining; ++c) {
                        int active = 0;
                        for (int j = 0; j < remaining; ++j) {
                            int pos = live[j];
                            if (upper[pos] / powers[c] % 5 >= level && lower[pos][c] < level)
                                testing[active++] = pos;
                        }
                        if (!active) continue;
                        for (int score = 6; score >= -1 && active; --score) {
                            int threshold = 8 - level - score;
                            if (threshold >= 6) break;
                            for (int m : grouped[c][score + 1]) {
                                const auto* row = products.data() + size_t(m) * roots;
                                int retained = 0;
                                for (int j = 0; j < active; ++j) {
                                    int pos = testing[j];
                                    if (row[partners[pos]] <= threshold) testing[retained++] = pos;
                                    else if (level == 0) upper[pos] = -1;
                                    else {
                                        int digit = upper[pos] / powers[c] % 5;
                                        upper[pos] -= (digit - level + 1) * powers[c];
                                    }
                                }
                                active = retained;
                                if (!active) break;
                            }
                        }
                        for (int j = 0; j < active; ++j) lower[testing[j]][c] = level;
                        int retained = 0;
                        for (int j = 0; j < remaining; ++j) {
                            int pos = live[j];
                            if (upper[pos] >= 0 && allowed[upper[pos]]) live[retained++] = pos;
                        }
                        remaining = retained;
                    }
                }
                if (!switched) {
                    surviving += remaining;
                    continue;
                }
                for (int j = 0; j < remaining; ++j) {
                    int s = partners[live[j]], key = profiles_count - 1;
                    GroupBounds bounds(4);
                    bool possible = true;
                    for (int m = 0; m < midpoints && possible; ++m) {
                        int bound = 8 - scores[m] - products[size_t(m) * roots + s];
                        possible = bounds.update(groups, allowed, powers, 5, m, bound, key);
                    }
                    surviving += possible;
                }
            }
        }
    }
    std::cout << "{\"contacts\":" << checked << ",\"minimal_profiles\":85,\"feasible_bodies\":"
              << surviving << "}\n";
    return checked == offset.back() && surviving == 0 ? 0 : 2;
}
