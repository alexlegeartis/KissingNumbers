// Exhaust all contact triangles and tetrahedra on the deleted roots.
#include "switch_bounds.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace fs = std::filesystem;
constexpr int roots = 13248;
constexpr int per_root = 48 * 6;
using Signature = std::array<uint8_t, 18>;

template<class T>
std::vector<T> read(const fs::path& path, size_t size) {
    std::vector<T> data(size);
    std::ifstream input(path, std::ios::binary);
    input.read(reinterpret_cast<char*>(data.data()), size * sizeof(T));
    if (!input) throw std::runtime_error("Incomplete input: " + path.string());
    return data;
}

struct Geometry {
    std::vector<int8_t> products;
    std::vector<int32_t> incident, strongest, colors;
    std::vector<uint8_t> triangle, tetrahedron;

    SwitchGroups groups;

    template<size_t K, bool Switched = false>
    bool possible(const std::array<int, K>& face) const {
        constexpr int bits = K == 3 ? 3 : 4;
        constexpr int mask = (1 << bits) - 1;
        const auto& radius = K == 3 ? triangle : tetrahedron;
        int key = K == 3 ? 0x3ffff : 0xaaaaaa;
        constexpr int base = 1 << bits;
        constexpr std::array<int,6> powers = {1, base, 1 << (2*bits), 1 << (3*bits),
                                               1 << (4*bits), 1 << (5*bits)};
        GroupBounds bounds(K == 3 ? 7 : 10);
        std::array<const int8_t*, K> rows;
        for (size_t i = 0; i < K; ++i) rows[i] = products.data() + size_t(face[i]) * roots;
        auto clears = [&](int m) {
            int limit = 4 * K;
            for (auto row : rows) limit -= row[m];
            if constexpr (Switched)
                return bounds.update(groups, radius, powers, base, m, limit, key);
            if (limit < 0) return false;
            int shift = bits * colors[m];
            if (limit < ((key >> shift) & mask))
                key = (key & ~(mask << shift)) | (limit << shift);
            return bool(radius[key]);
        };
        for (int root : face)
            for (int side = 0; side < 2; ++side)
                if (!clears(incident[2 * root + side])) return false;
        for (int rank = 0; rank < 48; ++rank)
            for (int root : face)
                for (int color = 0; color < 6; ++color)
                    if (!clears(strongest[root * per_root + 6 * rank + color])) return false;
        for (int m = 0; m < roots; ++m)
            if (!clears(m)) return false;
        return true;
    }

    Signature chord(int a, int b) const {
        std::array<int, 6> best, low, high;
        best.fill(-100);
        low.fill(9);
        high.fill(-9);
        for (int r : {a, b}) {
            for (int j = 0; j < per_root; ++j) {
                int m = strongest[r * per_root + j], c = colors[m];
                int x = products[size_t(a) * roots + m];
                int sum = x + products[size_t(b) * roots + m];
                if (sum > best[c]) {
                    best[c] = sum;
                    low[c] = high[c] = x;
                } else if (sum == best[c]) {
                    low[c] = std::min(low[c], x);
                    high[c] = std::max(high[c], x);
                }
            }
        }
        Signature result{};
        for (int c = 0; c < 6; ++c) {
            result[3*c] = 16 + best[c];
            result[3*c+1] = 8 + low[c];
            result[3*c+2] = 8 + high[c];
        }
        return result;
    }
};

int main(int argc, char** argv) {
    if (argc != 2 && argc != 3) return 1;
    fs::path dir = argv[1];
    Geometry geometry{
        read<int8_t>(dir / "products", size_t(roots) * roots),
        read<int32_t>(dir / "incident", 2 * roots),
        read<int32_t>(dir / "strongest", per_root * roots),
        read<int32_t>(dir / "colors", roots),
        read<uint8_t>(dir / "radius3", 1 << 18),
        read<uint8_t>(dir / "radius4", 1 << 24)
    };
    bool switched = argc == 3;
    if (switched) {
        geometry.groups.load(argv[2], geometry.colors);
        geometry.groups.fold(geometry.triangle, 8);
        geometry.groups.fold(geometry.tetrahedron, 16);
    }
    auto offsets = read<int32_t>(dir / "offsets", roots + 1);
    if (offsets.back() != 2966566) return 2;
    auto edges = read<int32_t>(dir / "edges", offsets.back());
    std::vector<Signature> signatures(switched ? 0 : edges.size());
    int64_t triangles = 0, tetrahedra = 0, feasible_triangles = 0, feasible_tetrahedra = 0;
    #pragma omp parallel for schedule(dynamic,1) reduction(+:triangles,tetrahedra,feasible_triangles,feasible_tetrahedra)
    for (int a = 0; a < roots; ++a) {
        std::vector<int> common(offsets[a + 1] - offsets[a]);
        for (int i = offsets[a]; i < offsets[a + 1]; ++i) {
            int b = edges[i];
            if (!switched) signatures[i] = geometry.chord(a, b);
            auto end = std::set_intersection(
                edges.begin() + i + 1, edges.begin() + offsets[a + 1],
                edges.begin() + offsets[b], edges.begin() + offsets[b + 1], common.begin()
            );
            for (auto current = common.begin(); current != end; ++current) {
                int c = *current;
                ++triangles;
                feasible_triangles += switched ? geometry.possible<3, true>({a, b, c})
                                               : geometry.possible<3>({a, b, c});
                auto next = current + 1;
                int j = offsets[c];
                while (next != end && j < offsets[c + 1]) {
                    if (*next == edges[j]) {
                        ++tetrahedra;
                        feasible_tetrahedra += switched ? geometry.possible<4, true>({a, b, c, *next})
                                                       : geometry.possible<4>({a, b, c, *next});
                        ++next;
                        ++j;
                    } else if (*next < edges[j]) ++next;
                    else ++j;
                }
            }
        }
    }
    if (!switched) {
        std::ofstream file(dir / "chords", std::ios::binary);
        file.write(reinterpret_cast<const char*>(signatures.data()), signatures.size() * sizeof(Signature));
        if (!file) return 3;
    }
    std::cout << "{\"triangles\":" << triangles << ",\"tetrahedra\":" << tetrahedra
              << ",\"feasible_triangles\":" << feasible_triangles
              << ",\"feasible_tetrahedra\":" << feasible_tetrahedra << "}\n";
    return feasible_triangles || feasible_tetrahedra;
}
