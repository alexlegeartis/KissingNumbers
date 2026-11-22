// Replace up to three of the eighteen usual signed fibres with rational heights.
#include "qr_shortening.hpp"
#include <map>

struct Fibre {
    int syndrome;
    int height; // in quarters of a unit body coordinate
};

class Heights {
    std::array<std::array<int, 21>, 21> threshold{};
public:
    Heights() {
        for (int h = -10; h <= 10; ++h) {
            for (int k = -10; k <= 10; ++k) {
                long long norm_product = (39 * 16 + h * h) * (39 * 16 + k * k);
                int d = 0;
                long long product = 39 * 16 + h * k;
                while (product > 0 && 4 * product * product > norm_product) {
                    ++d;
                    product -= 32;
                }
                threshold[h + 10][k + 10] = d;
            }
            assert(std::max(distance(h, -4), distance(h, 4)) >= 10);
            assert(distance(h, h) <= 12);
        }
    }
    int distance(int h, int k) const { return threshold[h + 10][k + 10]; }
};

bool extend(const std::vector<Fibre>& options, const std::vector<int>& available,
            const Heights& heights, const std::vector<uint8_t>& distance, int target,
            std::vector<int>& chosen) {
    int needed = target - static_cast<int>(chosen.size());
    if (needed == 0) return true;
    for (int i = 0; i + needed <= static_cast<int>(available.size()); ++i) {
        int v = available[i];
        chosen.push_back(v);
        std::vector<int> compatible;
        if (needed > 1) {
            for (size_t j = i + 1; j < available.size(); ++j) {
                auto a = options[v], b = options[available[j]];
                if (distance[a.syndrome ^ b.syndrome] >= heights.distance(a.height, b.height))
                    compatible.push_back(available[j]);
            }
        }
        if (extend(options, compatible, heights, distance, target, chosen)) return true;
        chosen.pop_back();
    }
    return false;
}

bool test_mark(const std::vector<uint8_t>& distance, const std::vector<int>& deep,
               const std::vector<int>& deleted, int mark, const Heights& heights) {
    std::vector<int> shell{0};
    for (int i = 0; i < 9; ++i)
        if (i != mark) shell.push_back(deleted[i]);
    std::vector<Fibre> base;
    for (int s : shell) base.push_back({s, -4});
    for (int s : shell) base.push_back({s ^ deleted[mark], 4});
    for (size_t i = 0; i < base.size(); ++i)
        for (size_t j = 0; j < i; ++j)
            assert(distance[base[i].syndrome ^ base[j].syndrome] >= heights.distance(base[i].height, base[j].height));

    // Each height needs distance at least ten from one of the two old levels.
    // At least six of the nine old fibres on that level survive any allowed exchange.
    std::vector<int> bodies;
    for (auto a : base)
        for (auto x : deep) bodies.push_back(a.syndrome ^ x);
    std::sort(bodies.begin(), bodies.end());
    bodies.erase(std::unique(bodies.begin(), bodies.end()), bodies.end());

    std::vector<Fibre> options;
    std::map<int, std::vector<int>> by_blockers;
    for (int x : bodies) {
        std::array<std::array<int, 13>, 2> misses{};
        for (int j = 0; j < 18; ++j) {
            int d = distance[x ^ base[j].syndrome];
            for (int need = 8; need <= 11; ++need)
                if (d < need) misses[j / 9][need] |= 1 << j;
        }
        for (int h = -10; h <= 10; ++h) {
            int lower = heights.distance(h, -4), upper = heights.distance(h, 4);
            assert(lower >= 8 && lower <= 11 && upper >= 8 && upper <= 11);
            int mask = misses[0][lower] | misses[1][upper];
            if (__builtin_popcount(mask) > 3) continue;
            auto is_original = [=](Fibre a) { return a.syndrome == x && a.height == h; };
            if (std::any_of(base.begin(), base.end(), is_original)) continue;
            by_blockers[mask].push_back(options.size());
            options.push_back({x, h});
        }
    }
    std::vector<int> removals{0};
    for (int j = 0; j < 18; ++j) {
        size_t end = removals.size();
        for (size_t k = 0; k < end; ++k)
            if (__builtin_popcount(removals[k]) < 3) removals.push_back(removals[k] | (1 << j));
    }
    for (int mask : removals) {
        std::vector<int> available;
        for (int sub = mask; ; sub = (sub - 1) & mask) {
            auto iter = by_blockers.find(sub);
            if (iter != by_blockers.end())
                available.insert(available.end(), iter->second.begin(), iter->second.end());
            if (sub == 0) break;
        }
        std::vector<int> chosen;
        if (extend(options, available, heights, distance, 2 + __builtin_popcount(mask), chosen))
            return true;
    }
    return false;
}

int main(int argc, char** argv) {
    Code code;
    auto cuts = representatives(code, 9);
    size_t count = argc > 1 ? std::min(cuts.size(), size_t(std::stoul(argv[1]))) : cuts.size();
    std::vector<uint8_t> initial(syndromes);
    for (int i = 0; i < syndromes; ++i) initial[i] = __builtin_popcount(i);
    Heights heights;
    size_t wins = 0, done = 0;
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < count; ++i) {
        auto columns = syndrome_columns(code, cuts[i]);
        auto distance = initial;
        for (int c : columns.extra) include_column(distance.data(), c);
        if (i == 0) cross_check_distances(columns.extra, distance);
        std::vector<int> deep;
        for (int x = 1; x < syndromes; ++x)
            if (distance[x] >= 10) deep.push_back(x);
        for (int mark = 0; mark < 9; ++mark) {
            if (test_mark(distance, deep, columns.deleted, mark, heights)) {
                #pragma omp critical
                { ++wins; std::cout << "win at cut " << cuts[i] << " mark " << mark << std::endl; }
            }
        }
        #pragma omp critical
        {
            ++done;
            if (done % 2000 == 0) std::cout << "checked " << done << std::endl;
        }
    }
    std::cout << "cuts " << done << " marked cases " << 9 * done << " wins " << wins << std::endl;
    return wins ? 1 : 0;
}
