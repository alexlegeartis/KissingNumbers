// Finite search only. The NPZ certificate is verified by two separate programs.
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

struct Support { uint64_t lo, hi; };
int pop(uint64_t x) { return __builtin_popcountll(x); }
int intersection(Support x, Support y) {
    return pop(x.lo & y.lo) + pop(x.hi & y.hi);
}
int section_weight(Support x) {
    const int weights[] = {16384, 0, 8192, 0, 6144, 0, 5120, 0, 8960};
    return weights[pop(x.lo & 255)];
}
int row(Support x, int r) {
    return r < 16 ? (x.lo >> (4*r)) & 15 : (x.hi >> (4*(r-16))) & 15;
}
void put(Support &x, int r, int v) {
    if (r < 16)
        x.lo = (x.lo & ~(15ULL << (4*r))) | (uint64_t(v) << (4*r));
    else
        x.hi = (x.hi & ~(15ULL << (4*(r-16)))) | (uint64_t(v) << (4*(r-16)));
}
bool ordinary(Support x) {
    for (int r = 0; r < 24; r++) {
        int v = row(x, r);
        if (v != 0 && v != 3 && v != 12 && v != 15) return false;
    }
    return true;
}
void save(const std::vector<Support> &family, const char *path) {
    std::ofstream stream(path);
    stream << family.size() << '\n';
    for (auto s : family) stream << s.lo << ' ' << s.hi << '\n';
    if (!stream) throw std::runtime_error("Cannot write checkpoint");
}

int main(int argc, char **argv) {
    if (argc != 6) {
        std::cerr << "usage: fill_correct INPUT OUTPUT SEED TRIALS NEUTRAL_PERIOD\n";
        return 2;
    }
    std::ifstream input(argv[1]);
    int n = 0;
    input >> n;
    if (!input || n < 1 || n > 1000000) return 2;
    std::vector<Support> family(n);
    for (auto &s : family) {
        input >> s.lo >> s.hi;
        int occupied = 0;
        for (int r = 0; r < 24; r++) occupied += row(s, r) != 0;
        if (!input || s.hi >= (1ULL << 32) ||
            pop(s.lo) + pop(s.hi) != 24 || occupied > 12) return 2;
    }
    std::mt19937_64 rng(std::stoull(argv[3]));
    uint64_t trials = std::stoull(argv[4]);
    unsigned neutral_period = std::stoul(argv[5]);
    if (!neutral_period) return 2;
    long long score = 0;
    for (auto s : family) score += section_weight(s);
    std::vector<std::vector<Support>> bucket(family.size());
    uint64_t swaps = 0, adds = 0, trades = 0;
    auto start = std::chrono::steady_clock::now();
    for (uint64_t it = 0; it < trials; it++) {
        int parent = rng() % family.size();
        if (ordinary(family[parent]) && rng() % 4) continue;
        Support x = family[parent];
        int mutations = 1 + rng() % 60;
        // Change pairings, or move a two-cell row into an empty row.
        // Both operations preserve weight 24 and at most twelve occupied rows.
        for (int k = 0; k < mutations; k++) {
            int r = rng() % 24, v = row(x, r);
            if (pop(v) != 2) continue;
            if (rng() % 12 == 0) {
                int dest = rng() % 24;
                if (row(x, dest) == 0) { put(x, dest, v); put(x, r, 0); }
            } else {
                const int mixed[] = {5, 6, 9, 10};
                put(x, r, mixed[rng() % 4]);
            }
        }
        if (x.lo == family[parent].lo && x.hi == family[parent].hi) continue;
        int blocker = -1;
        bool reject = false;
        for (int j = 0; j < int(family.size()); j++) {
            if (intersection(x, family[j]) <= 12) continue;
            if (blocker >= 0) { reject = true; break; }
            blocker = j;
        }
        if (reject) continue;
        if (blocker < 0) {
            family.push_back(x);
            bucket.emplace_back();
            score += section_weight(x);
            adds++;
            save(family, argv[2]);
            std::cout << "add " << it << " supports " << family.size()
                      << " score95 " << score << std::endl;
        } else {
            bool traded = false;
            // Candidates sharing one blocker may together replace it.
            // Recheck old bucket entries: earlier exchanges can invalidate them.
            for (Support y : bucket[blocker]) {
                if (intersection(x, y) > 12 || section_weight(x) + section_weight(y)
                    < section_weight(family[blocker])) continue;
                bool ok = true;
                for (int j = 0; j < int(family.size()); j++) {
                    if (j != blocker && intersection(y, family[j]) > 12) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) continue;
                score += section_weight(x) + section_weight(y)
                       - section_weight(family[blocker]);
                family[blocker] = x;
                family.push_back(y);
                trades++;
                traded = true;
                bucket[blocker].clear();
                bucket.emplace_back();
                save(family, argv[2]);
                std::cout << "trade " << it << " supports " << family.size()
                          << " score95 " << score << std::endl;
                break;
            }
            if (!traded) {
                auto &b = bucket[blocker];
                if (b.size() < 64) b.push_back(x);
                else b[rng() % b.size()] = x;
                // Occasional one-for-one moves explore a plateau. Never reduce
                // either the support count or the exact 95D section count.
                if (rng() % neutral_period == 0 &&
                    section_weight(x) >= section_weight(family[blocker])) {
                    score += section_weight(x) - section_weight(family[blocker]);
                    family[blocker] = x;
                    swaps++;
                    b.clear();
                }
            }
        }
        if ((it+1) % 100000 == 0) {
            save(family, argv[2]);
            std::cout << "trial " << it+1 << " swaps " << swaps << " adds " << adds
                      << " supports " << family.size() << " score95 " << score << std::endl;
        }
    }
    save(family, argv[2]);
    std::cout << "done swaps " << swaps << " adds " << adds << " trades " << trades
              << " supports " << family.size() << " score95 " << score << " seconds "
              << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
              << std::endl;
}
