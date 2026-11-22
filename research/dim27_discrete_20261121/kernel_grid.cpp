// Exhaust the even or full integer grid in the fixed 11D kernel.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

constexpr int dimension = 11;
constexpr int maximum_norm = 128;
using Row = std::array<int, dimension>;
using Counts = std::array<long long, 26>;
constexpr Row order = {2, 4, 0, 1, 3, 5, 6, 7, 8, 9, 10};

struct Pattern {
    Row absolute;
    std::vector<Row> signed_rows;
};

std::array<std::vector<Pattern>, dimension> constraints;
std::vector<Row> incumbent;
bool attainable[dimension + 1][maximum_norm + 1]{};
int limit;

int dot(const Row& a, const Row& b, int size) {
    int sum = 0;
    for (int i = 0; i < size; ++i) sum += a[i] * b[i];
    return sum;
}

bool allowed(const Row& x, int at) {
    for (const auto& pattern : constraints[at]) {
        int upper = 0;
        for (int j = 0; j <= at; ++j) upper += pattern.absolute[j] * std::abs(x[j]);
        if (upper <= limit) continue;
        // With all sign patterns present, the upper bound is attained.
        if (pattern.signed_rows.empty()) return false;
        for (const Row& row : pattern.signed_rows) {
            if (std::abs(dot(row, x, at + 1)) > limit) return false;
        }
    }
    return true;
}

void walk(int at, int remaining, Row& x, bool sign_chosen, Counts& counts) {
    if (at == dimension) {
        int blockers = 0;
        for (const Row& row : incumbent) {
            blockers += std::abs(dot(row, x, dimension)) > limit;
            if (blockers == 25) break;
        }
        ++counts[blockers];
        return;
    }
    const int cap = std::sqrt(remaining);
    for (int value = sign_chosen ? -cap : 0; value <= cap; ++value) {
        const int rest = remaining - value * value;
        if (!attainable[dimension - at - 1][rest]) continue;
        x[at] = value;
        if (allowed(x, at)) walk(at + 1, rest, x, sign_chosen || value != 0, counts);
    }
}

Row read_row(std::istream& file) {
    Row raw{}, sorted{};
    for (int& value : raw) file >> value;
    for (int i = 0; i < dimension; ++i) sorted[i] = raw[order[i]];
    return sorted;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: kernel_grid <1 for even grid or 2 for integer grid> <input>\n";
        return 2;
    }
    const int scale = std::atoi(argv[1]);
    if (scale != 1 && scale != 2) return 2;
    const int norm = 32 * scale * scale;
    limit = 32 * scale;
    std::ifstream file(argv[2]);
    int normal_count = 0, incumbent_count = 0;
    if (!(file >> normal_count >> incumbent_count) || normal_count != 2895 || incumbent_count != 238) {
        std::cerr << "Incorrect input header\n";
        return 2;
    }

    std::map<Row, std::vector<Row>> groups;
    for (int j = 0; j < normal_count; ++j) {
        const Row row = read_row(file);
        const auto first = std::find_if(row.begin(), row.end(), [](int x) { return x != 0; });
        if (first == row.end() || *first < 0) continue;
        Row absolute;
        std::transform(row.begin(), row.end(), absolute.begin(), [](int x) { return std::abs(x); });
        groups[absolute].push_back(row);
    }
    for (auto& [absolute, rows] : groups) {
        int last = dimension - 1;
        while (absolute[last] == 0) --last;
        const int support = std::count_if(absolute.begin(), absolute.end(), [](int x) { return x != 0; });
        if (rows.size() == (1u << (support - 1))) rows.clear();
        constraints[last].push_back({absolute, rows});
    }
    for (auto& group : constraints) {
        std::sort(group.begin(), group.end(), [](const Pattern& a, const Pattern& b) {
            return a.signed_rows.size() < b.signed_rows.size();
        });
    }
    for (int j = 0; j < incumbent_count; ++j) incumbent.push_back(read_row(file));
    if (!file) return 2;

    attainable[0][0] = true;
    for (int n = 1; n <= dimension; ++n) {
        for (int s = 0; s <= norm; ++s) {
            for (int v = 0; v * v <= s; ++v) attainable[n][s] |= attainable[n - 1][s - v * v];
        }
    }

    const int cap = std::sqrt(norm);
    std::vector<std::array<int, 3>> prefixes;
    for (int a = 0; a <= cap; ++a) {
        for (int b = a ? -cap : 0; b <= cap; ++b) {
            for (int c = a || b ? -cap : 0; c <= cap; ++c) {
                if (a * a + b * b + c * c <= norm) prefixes.push_back({a, b, c});
            }
        }
    }

    Counts total{};
    #pragma omp parallel
    {
        Counts local{};
        #pragma omp for schedule(dynamic, 4)
        for (size_t j = 0; j < prefixes.size(); ++j) {
            const auto [a, b, c] = prefixes[j];
            Row x{};
            x[0] = a;
            x[1] = b;
            x[2] = c;
            if (allowed(x, 0) && allowed(x, 1) && allowed(x, 2)) {
                walk(3, norm - a * a - b * b - c * c, x, a || b || c, local);
            }
        }
        #pragma omp critical
        for (size_t j = 0; j < total.size(); ++j) total[j] += local[j];
    }
    long long feasible = 0;
    for (long long count : total) feasible += count;
    std::cout << "{\n  \"scale\": " << scale
              << ",\n  \"lines_clearing_exterior\": " << feasible
              << ",\n  \"lines_clearing_incumbent_as_well\": " << total[0]
              << ",\n  \"blocker_histogram_capped_at_25\": [";
    for (size_t j = 0; j < total.size(); ++j) {
        if (j) std::cout << ", ";
        std::cout << total[j];
    }
    std::cout << "]\n}\n";
}
