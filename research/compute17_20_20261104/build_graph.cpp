// Exact graph construction; 16 integer coordinates plus a sqrt(2) axis.
#include <immintrin.h>
#include <omp.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

static inline int dot16(const int16_t* a, const int16_t* b) {
    __m256i p = _mm256_madd_epi16(_mm256_loadu_si256((const __m256i*)a),
                                _mm256_loadu_si256((const __m256i*)b));
    __m128i t = _mm_add_epi32(_mm256_castsi256_si128(p), _mm256_extracti128_si256(p,1));
    t = _mm_hadd_epi32(t,t); t = _mm_hadd_epi32(t,t);
    return _mm_cvtsi128_si32(t);
}

int main(int argc, char** argv) {
    if (argc != 4) { std::cerr << "usage: build_graph vertices.bin graph.bin threads\n"; return 2; }
    omp_set_num_threads(std::atoi(argv[3]));
    std::ifstream in(argv[1], std::ios::binary);
    uint32_t n,d; in.read((char*)&n,4); in.read((char*)&d,4);
    if (d != 17 || !in) return 3;
    std::vector<int16_t> v(size_t(n)*d);
    in.read((char*)v.data(),v.size()*2);
    if (!in) return 4;
    std::vector<std::vector<uint32_t>> adj(n);
    double start=omp_get_wtime();
    #pragma omp parallel for schedule(dynamic,32)
    for (uint32_t i=0; i<n; ++i) {
        auto a=&v[size_t(i)*17];
        auto& row=adj[i]; row.reserve(3000);
        for (uint32_t j=0;j<n;++j) {
            if(i==j) continue;
            auto b=&v[size_t(j)*17];
            // Opposite norm-six tiers never conflict, by Cauchy-Schwarz.
            if (a[16]*b[16]==-36) continue;
            if(dot16(a,b)+2*a[16]*b[16]>144) row.push_back(j);
        }
    }
    std::vector<uint64_t> off(n+1);
    uint32_t mind=n,maxd=0;
    for(uint32_t i=0;i<n;++i) {
        off[i+1]=off[i]+adj[i].size();
        mind=std::min(mind,uint32_t(adj[i].size())); maxd=std::max(maxd,uint32_t(adj[i].size()));
    }
    std::ofstream out(argv[2],std::ios::binary);
    out.write((char*)&n,4); uint32_t magic=0x4b4e4731;out.write((char*)&magic,4);
    out.write((char*)off.data(),off.size()*8);
    for(auto& row:adj)out.write((char*)row.data(),row.size()*4);
    if (!out) return 5;
    std::cout<<"{\"vertices\":"<<n<<",\"directed_edges\":"<<off[n]
             <<",\"min_degree\":"<<mind<<",\"max_degree\":"<<maxd
             <<",\"seconds\":"<<omp_get_wtime()-start<<"}"<<std::endl;
}
