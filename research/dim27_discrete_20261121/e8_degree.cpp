// Bound the number of E8 root lines usable by two compatible midpoint orbits.
#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

constexpr int roots = 120;

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    std::ifstream file(argv[1], std::ios::binary);
    int points = 0;
    std::array<int, 3> sizes{};
    auto read = [&](int& value) { file.read(reinterpret_cast<char*>(&value), sizeof(value)); };
    read(points);
    for (int& size : sizes) read(size);
    if (!file || points != 5696 || sizes != std::array<int, 3>{144, 275, 144}) return 2;

    std::vector<std::array<int, 2>> parents(points);
    for (auto& parent : parents) {
        for (int& root : parent) {
            read(root);
            if (root < 0 || root >= roots) return 2;
        }
        if (parent[0] == parent[1]) return 2;
    }
    std::vector<std::array<int, 11>> coordinates(points);
    for (auto& row : coordinates) for (int& coordinate : row) read(coordinate);
    std::array<std::vector<std::vector<unsigned char>>, 3> frames;
    for (int axis = 0; axis < 3; ++axis) {
        frames[axis].resize(sizes[axis], std::vector<unsigned char>(points));
        for (auto& mask : frames[axis]) file.read(reinterpret_cast<char*>(mask.data()), points);
    }
    if (!file) return 2;

    std::array<std::vector<int>, roots> incident;
    for (int i = 0; i < points; ++i) {
        for (int root : parents[i]) incident[root].push_back(i);
    }
    std::array<std::vector<std::array<int, 2>>, roots> compatible;
    for (int root = 0; root < roots; ++root) {
        for (int i : incident[root]) {
            for (int j : incident[root]) {
                if (i >= j) continue;
                int body = 0, tail = 0;
                for (int k = 0; k < 8; ++k) body += coordinates[i][k] * coordinates[j][k];
                for (int k = 8; k < 11; ++k) tail += coordinates[i][k] * coordinates[j][k];
                if (std::abs(body) + std::abs(tail) <= 64) compatible[root].push_back({i, j});
            }
        }
    }

    std::array<long long, roots + 1> histogram{};
    #pragma omp parallel
    {
        std::array<long long, roots + 1> local{};
        std::vector<int> active;
        std::vector<bool> allowed(points);
        #pragma omp for schedule(dynamic, 4)
        for (int pair = 0; pair < sizes[0] * sizes[2]; ++pair) {
            const int a = pair / sizes[2], c = pair % sizes[2];
            active.clear();
            for (int i = 0; i < points; ++i) {
                if (frames[0][a][i] && frames[2][c][i]) active.push_back(i);
            }
            for (const auto& middle : frames[1]) {
                std::array<int, roots> degree{};
                int bound = 0;
                for (int i : active) {
                    if (!middle[i]) continue;
                    for (int root : parents[i]) bound += (++degree[root] == 2);
                }
                // A root incident to two individually allowed choices only
                // counts when at least one such pair is jointly compatible.
                if (bound >= 65) {
                    std::fill(allowed.begin(), allowed.end(), false);
                    for (int i : active) if (middle[i]) allowed[i] = true;
                    bound = 0;
                    for (int root = 0; root < roots; ++root) {
                        if (degree[root] < 2) continue;
                        for (auto [i, j] : compatible[root]) {
                            if (allowed[i] && allowed[j]) {
                                ++bound;
                                break;
                            }
                        }
                    }
                }
                ++local[bound];
            }
        }
        #pragma omp critical
        for (int i = 0; i <= roots; ++i) histogram[i] += local[i];
    }
    long long total = 0, surviving = 0;
    int maximum = 0;
    for (int i = 0; i <= roots; ++i) {
        total += histogram[i];
        if (histogram[i]) maximum = i;
        if (i >= 65) surviving += histogram[i];
    }
    std::cout << "{\n  \"frame_triples\": " << total
              << ",\n  \"maximum_degree_two_upper_bound\": " << maximum
              << ",\n  \"frame_triples_passing_65\": " << surviving << "\n}\n";
    return 0;
}
