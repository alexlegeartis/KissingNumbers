#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

// Each group can swap its two perpendicular tail lines. For two bounds a,b,
// the smaller absolute projection must fit min(a,b) and the larger fits max(a,b).
// Taking the minimum of each sorted bound over groups is exact for one new tail.
struct SwitchGroups {
    std::array<std::array<int32_t, 2>, 3> colors{};
    std::array<int32_t, 32> plane{};
    std::vector<int> slot;
    std::vector<uint8_t> largest;

    void load(const std::filesystem::path& path, const std::vector<int32_t>& labels) {
        std::ifstream f(path, std::ios::binary);
        std::vector<int32_t> group(labels.size());
        f.read(reinterpret_cast<char*>(colors.data()), 6 * sizeof(int32_t));
        f.read(reinterpret_cast<char*>(plane.data()), 32 * sizeof(int32_t));
        f.read(reinterpret_cast<char*>(group.data()), group.size() * sizeof(int32_t));
        if (!f || f.peek() != std::ifstream::traits_type::eof())
            throw std::runtime_error("Invalid switch-group file");
        std::array<int, 6> frequencies{};
        for (auto pair : colors) for (int c : pair) ++frequencies.at(c);
        if (std::count(frequencies.begin(), frequencies.end(), 1) != 6)
            throw std::runtime_error("Invalid tail color partition");
        std::array<int, 32> sizes{};
        slot.resize(group.size());
        for (size_t i = 0; i < group.size(); ++i) {
            const int g = group[i], p = plane.at(g);
            const auto pair = colors.at(p);
            if (labels[i] != pair[0] && labels[i] != pair[1])
                throw std::runtime_error("Midpoint ordering disagrees with switch groups");
            slot[i] = 2 * g + (labels[i] == pair[1]);
            ++sizes[g];
        }
        // One group per plane supplies an inexpensive first pass in the shell scan.
        std::array<int, 3> dominant = {-1, -1, -1};
        for (int g = 0; g < 32; ++g) {
            auto& best = dominant.at(plane[g]);
            if (best < 0 || sizes[g] > sizes[best]) best = g;
        }
        largest.resize(group.size());
        for (size_t i = 0; i < group.size(); ++i)
            largest[i] = group[i] == dominant[plane[group[i]]];
    }

    template<class Table>
    void fold(Table& table, int base) const {
        std::array<int, 6> powers = {1};
        for (int i = 1; i < 6; ++i) powers[i] = base * powers[i - 1];
        for (auto pair : colors) {
            const int a = powers[pair[0]], b = powers[pair[1]];
            for (int key = 0; key < static_cast<int>(table.size()); ++key) {
                const int swapped = key + ((key / a % base) - (key / b % base)) * (b - a);
                table[key] |= table[swapped];
            }
        }
    }
};

struct GroupBounds {
    std::array<uint8_t, 64> limits;
    explicit GroupBounds(int maximum) { limits.fill(maximum); }

    template<class Table>
    bool update(const SwitchGroups& groups, const Table& allowed, const std::array<int, 6>& powers,
                int base, int row, int bound, int& key) {
        if (bound < 0) return false;
        int slot = groups.slot[row];
        if (bound >= limits[slot]) return true;
        limits[slot] = bound;
        const int other = limits[slot ^ 1];
        auto pair = groups.colors[groups.plane[slot / 2]];
        std::array<int, 2> sorted = {std::min(bound, other), std::max(bound, other)};
        for (int i = 0; i < 2; ++i) {
            int place = powers[pair[i]], current = key / place % base;
            if (sorted[i] < current) key -= (current - sorted[i]) * place;
        }
        return allowed[key];
    }
};
