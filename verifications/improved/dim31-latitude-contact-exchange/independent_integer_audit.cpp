// Independent full point-pair enumeration. No antipodal compression, BLAS,
// floating point dot products, search imports, or previous verifier imports.
// Input: contiguous signed-int16 rows padded with zero to 32 coordinates.
// Each real coordinate is checked in [-8,8], so every pairwise product is
// in [-64,64] and every partial 32-term sum is in [-2048,2048]. All SIMD
// multiplication/addition lanes and reductions are exact signed int32.
#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || defined(_M_X64)
#include <immintrin.h>
#endif
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

static inline int dot32(const int16_t* a, const int16_t* b) {
#if defined(__AVX512BW__) && defined(__AVX512F__)
    const __m512i x = _mm512_loadu_si512(a);
    const __m512i y = _mm512_loadu_si512(b);
    const __m512i products = _mm512_madd_epi16(x, y);
    return _mm512_reduce_add_epi32(products);
#else
    int sum = 0;
    for (int k = 0; k < 32; ++k) sum += int(a[k]) * int(b[k]);
    return sum;
#endif
}

int main(int argc, char** argv) {
    if (argc != 4) throw std::runtime_error("usage: audit input.bin row_count retained_count");
    const size_t n = std::stoull(argv[2]), retained = std::stoull(argv[3]);
    if (retained > n) throw std::runtime_error("retained count");
    std::vector<int16_t> x(n * 32);
    std::ifstream file(argv[1], std::ios::binary);
    file.read(reinterpret_cast<char*>(x.data()), x.size() * sizeof(int16_t));
    if (!file || file.peek() != std::char_traits<char>::eof()) throw std::runtime_error("file size/read");
    for (size_t i = 0; i < n; ++i) {
        int norm = 0;
        for (size_t k = 0; k < 32; ++k) {
            const int v = x[i*32+k];
            if (v < -8 || v > 8 || (k == 31 && v != 0)) throw std::runtime_error("coefficient");
            norm += v*v;
        }
        if (norm != 128 || dot32(x.data()+32*i, x.data()+32*i) != 128)
            throw std::runtime_error("norm");
    }
    // Category 0 retained/retained, 1 retained/added, 2 added/added.
    uint64_t count[3] = {}, contacts[3] = {}, violations[3] = {};
    int maximum[3] = {-128,-128,-128}, minimum[3] = {128,128,128};
    const auto start = std::chrono::steady_clock::now();
    for (size_t i = 0; i < n; ++i) {
        const int16_t* a = x.data()+32*i;
        for (size_t j = i+1; j < n; ++j) {
            const int category = i >= retained ? 2 : (j >= retained ? 1 : 0);
            const int dot = dot32(a, x.data()+32*j);
            ++count[category];
            contacts[category] += dot == 64;
            violations[category] += dot > 64;
            maximum[category] = std::max(maximum[category], dot);
            minimum[category] = std::min(minimum[category], dot);
        }
        if ((i+1)%8192 == 0) {
            const auto sec = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            std::cerr << "rows=" << i+1 << "/" << n << " pairs="
                      << count[0]+count[1]+count[2] << " seconds=" << sec << '\n';
        }
    }
    const auto sec = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout << "{\"rows\":" << n << ",\"retained\":" << retained << ",\"seconds\":" << sec << ",\"categories\":[";
    for (int c=0;c<3;++c) {
        if(c) std::cout << ',';
        std::cout << "{\"pairs\":" << count[c] << ",\"contacts\":" << contacts[c]
                  << ",\"violations\":" << violations[c] << ",\"maximum\":" << maximum[c]
                  << ",\"minimum\":" << minimum[c] << '}';
    }
    std::cout << "]}\n";
    return (violations[0]+violations[1]+violations[2]) ? 1 : 0;
}
