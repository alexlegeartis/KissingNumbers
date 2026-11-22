#pragma once

#ifdef NDEBUG
#error "The enumeration uses assertions to certify its inputs and orbit coverage."
#endif
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <immintrin.h>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

using Word = uint64_t;
using Permutation = std::array<uint8_t, 48>;
constexpr int syndromes = 1 << 24;
constexpr int graph_limit = 2048;
using VertexSet = std::bitset<graph_limit>;

Word bit(int i) { return Word{1} << i; }
int highest_bit(Word x) { return 63 - __builtin_clzll(x); }

Word permute(Word x, const Permutation& permutation) {
    Word result = 0;
    while (x) {
        result |= bit(permutation[__builtin_ctzll(x)]);
        x &= x - 1;
    }
    return result;
}

struct Code {
    std::array<Word, 48> pivots{}, columns{};

    Word reduce(Word word) const {
        for (int i = 47; i >= 0; --i)
            if (word & bit(i)) word ^= pivots[i];
        return word;
    }

    Code() {
        const Word mask = bit(47) - 1;
        Word residues = 0;
        for (int i = 1; i < 47; ++i) residues |= bit(i * i % 47);
        for (int i = 0; i < 47; ++i) {
            Word orbit_word = ((residues << i) | (residues >> (47 - i))) & mask;
            Word row = reduce(orbit_word | bit(47));
            if (row) pivots[highest_bit(row)] = row;
        }

        Word free_columns = 0;
        std::vector<Word> basis;
        for (int j = 0; j < 48; ++j) {
            if (pivots[j]) basis.push_back(pivots[j]);
            else free_columns |= bit(j);
        }
        assert(basis.size() == 24);
        Word word = 0;
        int minimum = 48;
        for (int i = 1; i < syndromes; ++i) {
            word ^= basis[__builtin_ctz(i)];
            minimum = std::min(minimum, __builtin_popcountll(word));
        }
        assert(minimum == 12);
        for (int j = 0; j < 48; ++j)
            columns[j] = _pext_u64(reduce(bit(j)), free_columns);
    }
};

std::vector<Word> representatives(const Code& code, int removed) {
    assert(removed == 8 || removed == 9);
    Permutation identity{}, translation{}, inversion{};
    std::iota(identity.begin(), identity.end(), 0);
    for (int j = 0; j < 47; ++j) {
        translation[j] = (j + 1) % 47;
        for (int k = 1; k < 47; ++k)
            if (j * k % 47 == 46) inversion[j] = k;
    }
    translation[47] = 47;
    inversion[0] = 47;
    inversion[47] = 0;
    const std::array<Permutation, 2> generators{translation, inversion};
    for (const auto& generator : generators)
        for (auto row : code.pivots)
            assert(code.reduce(permute(row, generator)) == 0);

    std::vector<Permutation> group{identity};
    std::set<Permutation> seen{identity};
    for (size_t i = 0; i < group.size(); ++i) {
        for (const auto& generator : generators) {
            Permutation product;
            for (int j = 0; j < 48; ++j) product[j] = generator[group[i][j]];
            if (seen.insert(product).second) group.push_back(product);
        }
    }
    assert(group.size() == 51888);

    // Index group elements by the ordered preimages of 0, 1 and infinity.
    std::vector<int> triple(48 * 48 * 48, -1);
    for (size_t i = 0; i < group.size(); ++i) {
        Permutation inverse;
        for (int j = 0; j < 48; ++j) inverse[group[i][j]] = j;
        int key = (inverse[0] * 48 + inverse[1]) * 48 + inverse[47];
        assert(triple[key] == -1);
        triple[key] = i;
    }

    std::unordered_set<Word> marked;
    marked.reserve(removed == 8 ? 1400000 : 8500000);
    std::vector<Word> cuts;
    uint64_t covered = 0;
    for (Word subset = bit(removed - 3) - 1; subset < bit(45); ) {
        Word cut = bit(47) | 3 | (subset << 2);
        if (!marked.count(cut)) {
            cuts.push_back(cut);
            std::vector<int> coordinates;
            for (int j = 0; j < 48; ++j)
                if (cut & bit(j)) coordinates.push_back(j);
            unsigned stabilizer = 0;
            for (int a : coordinates) for (int b : coordinates) for (int c : coordinates) {
                int k = triple[(48 * a + b) * 48 + c];
                if (k < 0) continue;
                Word image = permute(cut, group[k]);
                marked.insert(image);
                stabilizer += image == cut;
            }
            assert(stabilizer && group.size() % stabilizer == 0);
            covered += group.size() / stabilizer;
        }
        // Next fixed-size subset in increasing integer order.
        Word low = subset & -subset, next = subset + low;
        subset = next | (((subset ^ next) >> 2) / low);
    }
    if (removed == 8) {
        assert(marked.size() == 1221759 && covered == 377348994 && cuts.size() == 7497);
    } else {
        assert(marked.size() == 8145060 && covered == 1677106640 && cuts.size() == 32345);
    }
    return cuts;
}

struct SyndromeColumns {
    std::vector<int> extra;
    std::vector<int> deleted;
};

// Choose coordinates in which the first 24 independent remaining columns are units.
SyndromeColumns syndrome_columns(const Code& code, Word cut) {
    int rank = 0;
    std::array<int, 24> pivots{}, labels{};
    SyndromeColumns columns;
    for (int i = 0; i < 48; ++i) {
        if (cut & bit(i)) continue;
        int v = code.columns[i], label = 0;
        while (v && pivots[highest_bit(v)]) {
            int j = highest_bit(v);
            v ^= pivots[j];
            label ^= labels[j];
        }
        if (v) {
            int j = highest_bit(v);
            pivots[j] = v;
            labels[j] = label ^ (1 << rank++);
        } else {
            columns.extra.push_back(label);
        }
    }
    assert(rank == 24);
    for (int i = 0; i < 48; ++i) {
        if (!(cut & bit(i))) continue;
        int v = code.columns[i], label = 0;
        while (v) {
            int j = highest_bit(v);
            v ^= pivots[j];
            label ^= labels[j];
        }
        columns.deleted.push_back(label);
    }
    assert(columns.extra.size() + columns.deleted.size() == 24);
    return columns;
}

// In-place SIMD implementation of d'(x) = min(d(x), 1 + d(x xor column)).
void include_column(uint8_t* distances, int column) {
    alignas(32) uint8_t indices[32];
    for (int j = 0; j < 32; ++j) indices[j] = (j ^ column) & 15;
    auto perm = _mm256_load_si256(reinterpret_cast<const __m256i*>(indices));
    auto shuffle = [&](const __m256i value) {
        auto shuffled = _mm256_shuffle_epi8(value, perm);
        return column & 16 ? _mm256_permute2x128_si256(shuffled, shuffled, 1) : shuffled;
    };
    auto* blocks = reinterpret_cast<__m256i*>(distances);
    const auto one = _mm256_set1_epi8(1);
    for (int i = 0; i < syndromes / 32; ++i) {
        int j = i ^ (column >> 5);
        if (j < i) continue;
        auto a = _mm256_loadu_si256(blocks + i);
        auto b = _mm256_loadu_si256(blocks + j);
        _mm256_storeu_si256(blocks + i, _mm256_min_epu8(a, _mm256_add_epi8(shuffle(b), one)));
        if (i != j)
            _mm256_storeu_si256(blocks + j, _mm256_min_epu8(b, _mm256_add_epi8(shuffle(a), one)));
    }
}

void cross_check_distances(std::vector<int> generators, const std::vector<uint8_t>& expected) {
    for (int i = 0; i < 24; ++i) generators.push_back(1 << i);
    std::vector<uint8_t> distances(syndromes, 255);
    std::vector<int> queue(syndromes);
    distances[0] = 0;
    int end = 1;
    for (int i = 0; i < end; ++i) {
        int x = queue[i];
        for (int generator : generators) {
            int y = x ^ generator;
            if (distances[y] == 255) {
                distances[y] = distances[x] + 1;
                queue[end++] = y;
            }
        }
    }
    assert(end == syndromes && distances == expected);
}
