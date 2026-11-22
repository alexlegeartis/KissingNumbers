// Enumerate translates of every eight-coordinate shortening of the binary QR(48) code.
// Build with GNU libstdc++, AVX2, BMI2 and OpenMP, keeping assertions enabled.
#include "qr_shortening.hpp"
#include <unordered_map>

struct Clique {
    std::vector<VertexSet> neighbours;
    std::vector<int> chosen;

    bool find(VertexSet available) {
        // A tenth coset has been translated to zero.
        if (chosen.size() == 9) return true;
        std::vector<int> order, colours;
        VertexSet left = available;
        int colour = 0;
        while (left.any()) {
            ++colour;
            VertexSet independent = left;
            while (independent.any()) {
                int v = independent._Find_first();
                order.push_back(v);
                colours.push_back(colour);
                left.reset(v);
                independent.reset(v);
                independent &= ~neighbours[v];
            }
        }
        for (int j = int(order.size()) - 1; j >= 0; --j) {
            if (chosen.size() + colours[j] < 9) return false;
            int v = order[j];
            chosen.push_back(v);
            if (find(available & neighbours[v])) return true;
            chosen.pop_back();
            available.reset(v);
        }
        return false;
    }
};


// Choose distinct targets at forbidden Hamming distance by augmenting paths.
bool match(int candidate, const std::vector<int>& blockers,
           std::vector<int>& owner, int& visited) {
    int options = blockers[candidate] & ~visited;
    while (options) {
        int j = __builtin_ctz(options);
        options &= options - 1;
        if (visited & (1 << j)) continue;
        visited |= 1 << j;
        if (owner[j] == -1 || match(owner[j], blockers, owner, visited)) {
            owner[j] = candidate;
            return true;
        }
    }
    return false;
}

size_t check_partial_refills(const std::vector<int>& deleted, const std::vector<int>& deep,
                           const std::vector<uint8_t>& distances) {
    std::vector<int> standard = deleted;
    standard.push_back(0);
    std::unordered_map<int, int> compatible;
    for (int i = 0; i < 9; ++i)
        for (int x : deep) compatible[x ^ standard[i]] |= 1 << i;

    // Translate an arbitrary opposite-parity pair to zero and an odd deep
    // syndrome. Remaining candidate cosets come in opposite-parity pairs.
    std::vector<int> odd;
    for (int x : deep)
        if (distances[x] % 2) odd.push_back(x);
    for (int y : odd) {
        std::vector<int> candidates;
        for (int x : odd)
            if (distances[x ^ y] >= 10) candidates.push_back(x);
        assert(odd.size() <= 10);
        std::vector<int> blockers(candidates.size());
        for (size_t i = 0; i < candidates.size(); ++i)
            for (size_t j = 0; j < candidates.size(); ++j)
                if (distances[candidates[i] ^ candidates[j] ^ y] < 10)
                    blockers[i] |= 1 << j;
        std::vector<int> owner(candidates.size(), -1);
        int bound = 2 + 2 * int(candidates.size());
        for (size_t i = 0; i < candidates.size(); ++i) {
            int visited = 0;
            if (match(int(i), blockers, owner, visited)) --bound;
        }
        if (bound > 9) throw std::runtime_error("opposite-parity matching failed");
    }
    std::vector<int> blockers;
    int units = 0;
    for (const auto& [x, mask] : compatible) {
        if (__builtin_popcount(unsigned(mask)) < 3) continue;
        if (std::find(standard.begin(), standard.end(), x) != standard.end()) continue;
        if (distances[x] == 1) {
            ++units;
        } else {
            blockers.push_back(511 ^ mask);
        }
    }
    // The zero coset and its one-bit translates partition into radius-one
    // balls about its codewords. Each ball contains at most one selected word.
    assert(units == 40);
    for (int retained = 0; retained < 512; ++retained) {
        if (__builtin_popcount(unsigned(retained)) != 3) continue;
        std::vector<int> eligible;
        for (int mask : blockers)
            if (!(mask & retained)) eligible.push_back(mask & 255);
        std::vector<int> owner(8, -1);
        for (size_t i = 0; i < eligible.size(); ++i) {
            int visited = 0;
            // A minimum error vector matches the entire candidate coset
            // bijectively to its target, at forbidden Hamming distance.
            if (!match(int(i), eligible, owner, visited))
                throw std::runtime_error("partial-refill matching failed");
        }
    }
    return odd.size();
}

struct Result {
    size_t deep_count, opposite_count;
    std::vector<int> witness;
};

Result solve(const Code& code, Word cut, const std::vector<uint8_t>& initial, bool cross_check) {
    auto columns = syndrome_columns(code, cut);
    auto distances = initial;
    for (int column : columns.extra) include_column(distances.data(), column);
    if (cross_check) cross_check_distances(columns.extra, distances);

    std::vector<int> deep;
    for (int x = 1; x < syndromes; ++x)
        if (distances[x] >= 10) deep.push_back(x);
    size_t opposite = check_partial_refills(columns.deleted, deep, distances);
    if (deep.size() > graph_limit) throw std::runtime_error("increase graph_limit");
    Clique graph;
    graph.neighbours.resize(deep.size());
    VertexSet available;
    for (size_t i = 0; i < deep.size(); ++i) {
        available.set(i);
        for (size_t j = 0; j < i; ++j) {
            if (distances[deep[i] ^ deep[j]] >= 10) {
                graph.neighbours[i].set(j);
                graph.neighbours[j].set(i);
            }
        }
    }
    Result result{deep.size(), opposite, {}};
    if (graph.find(available))
        for (int v : graph.chosen) result.witness.push_back(deep[v]);
    return result;
}

int main(int argc, char** argv) {
    Code code;
    auto cuts = representatives(code, 8);
    size_t count = argc > 1 ? std::min(cuts.size(), size_t(std::stoul(argv[1]))) : cuts.size();
    std::cout << "orbits " << cuts.size() << " running " << count << std::endl;
    std::vector<uint8_t> initial(syndromes);
    for (int i = 0; i < syndromes; ++i) initial[i] = __builtin_popcount(i);
    size_t done = 0, wins = 0, largest = 0, opposite = 0;
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < count; ++i) {
        Result result = solve(code, cuts[i], initial, i == 0);
        #pragma omp critical
        {
            ++done;
            opposite += result.opposite_count;
            largest = std::max(largest, result.deep_count);
            if (!result.witness.empty()) {
                ++wins;
                std::cout << "WIN " << cuts[i];
                for (int x : result.witness) std::cout << ' ' << x;
                std::cout << std::endl;
            }
            if (done % 250 == 0) std::cout << "checked " << done << std::endl;
        }
    }
    std::cout << "checked " << done << " winning orbits " << wins
              << " maximum deep syndromes " << largest << std::endl;
    std::cout << "three-retained partial-refill cases " << 84 * done << std::endl;
    std::cout << "opposite-parity partial-refill cases " << opposite << std::endl;
    return wins ? 1 : 0;
}
