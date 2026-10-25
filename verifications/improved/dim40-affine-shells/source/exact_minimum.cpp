// Exhaustive, integer-only enumeration in an enclosing triangular ellipsoid.
// verify.py proves that every integer vector of Gram norm <= 4 is included.
// Compile WITHOUT -DNDEBUG. No floating-point decision is used.
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
using I = __int128_t;
using U = __uint128_t;
constexpr int N = 40;
long long G[N][N], M[N][N], D[N], z[N], P;
uint64_t nodes = 0, leaves = 0, short_vectors = 0;

I parse(const std::string& s) {
    I v = 0;
    for (char c : s) {
        assert(c >= '0' && c <= '9' && v < (I(1)<<116));
        v = 10*v + c-'0';
    }
    return v;
}
long long floor_div(long long a, long long b) {
    long long q=a/b, r=a%b;
    return q-(r<0);
}
long long ceil_div(long long a, long long b) { return -floor_div(-a,b); }
uint64_t isqrt(U n) {
    // The caller guarantees n < 2^120. Binary search is entirely integral.
    uint64_t lo=0, hi=(uint64_t(1)<<60);
    while (lo < hi) {
        uint64_t m=lo+(hi-lo+1)/2;
        if (U(m)*m <= n) lo=m; else hi=m-1;
    }
    return lo;
}
void visit(int k, I remaining) {
    ++nodes;
    if (k < 0) {
        ++leaves;
        I norm=0;
        bool nonzero=false;
        for (int i=0;i<N;++i) {
            nonzero |= z[i]!=0;
            for (int j=0;j<N;++j) norm+=I(z[i])*G[i][j]*z[j];
        }
        assert(norm>=0);
        if (nonzero && norm<=4) ++short_vectors;
        return;
    }
    I offset=0;
    for (int j=k+1;j<N;++j) offset+=I(M[j][k])*z[j];
    assert(offset>-(I(1)<<60) && offset<(I(1)<<60));
    assert(remaining>=0 && remaining<(I(1)<<120));
    const uint64_t b=isqrt(U(remaining/D[k]));
    assert(b<(uint64_t(1)<<60));
    const long long off=(long long)offset;
    const long long lo=ceil_div(-(long long)b-off,P);
    const long long hi=floor_div((long long)b-off,P);
    assert(lo>-(1LL<<20) && hi<(1LL<<20));
    for (long long v=lo;v<=hi;++v) {
        z[k]=v;
        const I t=I(P)*v+off;
        assert(t>-(I(1)<<44) && t<(I(1)<<44));
        const I used=I(D[k])*t*t;
        assert(used>=0 && used<=remaining);
        visit(k-1,remaining-used);
    }
}
int main(int argc, char** argv) {
    if (argc!=2) return 2;
    std::ifstream in(argv[1]);
    int n;
    std::string radius_text;
    in>>n>>P>>radius_text;
    assert(in && n==N && P>0 && P<(1LL<<40));
    const I radius=parse(radius_text);
    assert(radius>0 && radius<(I(1)<<120));
    for (auto& x:D) { in>>x; assert(x>0 && x<(1LL<<36)); }
    for (auto& row:M) for (auto& x:row) {
        in>>x; assert(x>-(1LL<<40) && x<(1LL<<40));
    }
    for (auto& row:G) for (auto& x:row) {
        in>>x; assert(x>-(1LL<<20) && x<(1LL<<20));
    }
    assert(in);
    const auto start=std::chrono::steady_clock::now();
    visit(N-1,radius);
    const double seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
    std::cout<<"{\"nonzero_vectors_of_norm_at_most_4\":"<<short_vectors
             <<",\"nodes\":"<<nodes<<",\"leaves\":"<<leaves
             <<",\"seconds\":"<<seconds<<"}\n";
    return short_vectors ? 1 : 0;
}
