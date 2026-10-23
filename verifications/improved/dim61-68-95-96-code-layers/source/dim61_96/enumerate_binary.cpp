// Exhaustive Gray-code enumeration of a binary generator with <=128 columns.
// Input file: n k, then k pairs of hexadecimal low/high 64-bit words.
// Optional output file receives the weight-16 words as consecutive lo,hi pairs.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    std::ifstream input(argv[1]);
    unsigned n, k;
    input >> n >> k;
    if (!input || n > 128 || k >= 63) return 3;
    std::vector<std::array<uint64_t, 2>> generator(k), selected;
    for (auto& row : generator) input >> std::hex >> row[0] >> row[1];
    if (!input) return 4;
    std::array<uint64_t, 129> histogram{};
    histogram[0] = 1;
    uint64_t lo = 0, hi = 0, total = uint64_t(1) << k;
    for (uint64_t i = 1; i < total; ++i) {
        const unsigned bit = __builtin_ctzll(i);
        lo ^= generator[bit][0];
        hi ^= generator[bit][1];
        const unsigned weight = __builtin_popcountll(lo) + __builtin_popcountll(hi);
        ++histogram[weight];
        if (weight == 16 && argc >= 3) selected.push_back({lo, hi});
    }
    if (argc >= 3) {
        std::ofstream out(argv[2], std::ios::binary);
        out.write(reinterpret_cast<const char*>(selected.data()), selected.size() * 16);
        if (!out) return 5;
    }
    std::cout << "{\"n\":" << n << ",\"k\":" << k << ",\"histogram\":{";
    bool first = true;
    for (unsigned i = 0; i <= n; ++i) if (histogram[i]) {
        if (!first) std::cout << ',';
        first = false;
        std::cout << '"' << i << "\":" << histogram[i];
    }
    std::cout << "}}\n";
}
