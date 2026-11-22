// Classify arbitrary real tails over every chord joining two deleted Leech roots.
#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;
using Vec = std::array<int, 3>;
using Six = std::array<int, 6>;
using Signature = std::array<uint8_t, 21>;
constexpr int roots = 13248;
constexpr std::array<Vec, 6> normals = {{
    {0,1,-1}, {0,1,1}, {1,-1,0}, {1,0,-1}, {1,0,1}, {1,1,0}
}};

template<class T>
std::vector<T> read(const fs::path& path, size_t size) {
    std::vector<T> data(size);
    std::ifstream input(path, std::ios::binary);
    input.read(reinterpret_cast<char*>(data.data()), sizeof(T) * size);
    if (!input) throw std::runtime_error("Incomplete input: " + path.string());
    return data;
}

struct Rational { int num, den; };
bool before(Rational a, Rational b) { return int64_t(a.num)*b.den < int64_t(b.num)*a.den; }
bool restrict_interval(int offset, int slope, Rational& lo, Rational& hi) {
    if (!slope) return offset >= 0;
    if (slope > 0) lo = std::max(lo, Rational{-offset,slope}, before);
    else hi = std::min(hi, Rational{offset,-slope}, before);
    return !before(hi, lo);
}

int dot(Vec a, Vec b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
Vec cross(Vec a, Vec b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
struct Vertex { Vec planes; std::array<Vec,3> columns; int den; };
std::vector<Vertex> make_vertices() {
    std::vector<Vertex> result;
    for (int i=0; i<6; ++i) for (int j=0; j<i; ++j) for (int k=0; k<j; ++k) {
        Vec c0=cross(normals[j],normals[k]), c1=cross(normals[k],normals[i]);
        Vec c2=cross(normals[i],normals[j]);
        int det=dot(normals[i],c0);
        if (!det) continue;
        for (int a : {-1,1}) for (int b : {-1,1}) for (int c : {-1,1}) {
            Vertex vertex{{i,j,k}, {c0,c1,c2}, std::abs(det)};
            Vec signs{a,b,c};
            for (int col=0; col<3; ++col) for (int& x : vertex.columns[col])
                x *= signs[col] * (det>0 ? 1 : -1);
            result.push_back(vertex);
        }
    }
    return result;
}
const auto vertices = make_vertices();
enum Result { empty, axial_midpoint, other };

Result examine(const Six& base, const Six& slope, Rational lo, Rational hi, int factor, bool first) {
    for (int i=0; i<6; ++i)
        if (!restrict_interval(base[i],slope[i],lo,hi)) return empty;
    Result best=empty;
    for (const auto& vertex : vertices) {
        Vec x{}, y{};
        for (int i=0; i<3; ++i) for (int k=0; k<3; ++k) {
            x[k] += vertex.columns[i][k]*base[vertex.planes[i]];
            y[k] += vertex.columns[i][k]*slope[vertex.planes[i]];
        }
        Rational left=lo, right=hi;
        bool valid=true;
        for (int j=0; j<6 && valid; ++j) for (int s : {-1,1}) {
            valid=restrict_interval(vertex.den*base[j]-s*dot(x,normals[j]),
                                    vertex.den*slope[j]-s*dot(y,normals[j]), left,right);
            if (!valid) break;
        }
        if (!valid) continue;
        for (auto f : {left,right}) {
            int64_t norm=0;
            int support=0;
            for (int k=0; k<3; ++k) {
                int64_t v=x[k]*f.den+y[k]*f.num;
                norm+=v*v;
                support+=(v!=0);
            }
            int64_t target=int64_t(factor)*f.num*(f.den-f.num)*vertex.den*vertex.den;
            if (norm < target) continue;
            if (norm > target || factor!=16 || 2*f.num!=f.den || support!=1) return other;
            best=axial_midpoint;
            if (first) return best;
        }
    }
    return best;
}

// Within each interval the largest body product in each of the six groups is affine.
using Line = std::pair<int,int>;  // products at lambda = 1 and 0
using Groups = std::array<std::vector<Line>,6>;
Result intervals(const Groups& groups, Rational lo, Rational hi, int factor, bool first) {
    if (before(hi,lo)) return empty;
    std::vector<Rational> breaks{lo,hi};
    for (const auto& group : groups) for (auto [a,b] : group) for (auto [c,d] : group) {
        if (a<=c || b>=d) continue;
        Rational f{d-b, a-c+d-b};
        if (before(lo,f) && before(f,hi)) breaks.push_back(f);
    }
    std::sort(breaks.begin(),breaks.end(),before);
    Result result=empty;
    for (size_t i=1; i<breaks.size(); ++i) {
        auto l=breaks[i-1], r=breaks[i];
        int n=l.num*r.den+r.num*l.den, d=2*l.den*r.den;
        Six base{}, slope{};
        for (int c=0; c<6; ++c) {
            int highest=-1000000;
            for (auto [x,y] : groups[c]) {
                int score=n*x+(d-n)*y;
                if (score>highest) { highest=score; base[c]=4-y; slope[c]=y-x; }
            }
        }
        result=std::max(result,examine(base,slope,l,r,factor,first));
        if (result==other || (first && result!=empty)) return result;
    }
    return result;
}

bool coarse(const Signature& sig) {
    Rational lo{2,14-sig[18]}, hi{12-sig[19],14-sig[19]}, mid{1,2};
    int factor=16-4*(int(sig[20])-4);
    for (int h=0; h<2; ++h) {
        Groups groups;
        for (int c=0; c<6; ++c) {
            int a=int(sig[3*c+1+h])-8;
            int b=int(sig[3*c])-16-a;
            groups[c]={{a,b},{2,2}};
        }
        auto left=h ? std::max(lo,mid,before) : lo;
        auto right=h ? hi : std::min(hi,mid,before);
        if (intervals(groups,left,right,factor,true)!=empty) return true;
    }
    return false;
}

Result full(const int8_t* a, const int8_t* b, const std::vector<int>& colors, int factor) {
    std::array<std::array<int,13>,6> maximum;
    for (auto& row : maximum) { row.fill(-100); row[8]=2; } // The poles supply the line (2,2).
    for (int j=0; j<roots; ++j) {
        int& high=maximum[colors[j]][a[j]+6];
        high=std::max(high,int(b[j]));
    }
    Groups groups;
    for (int c=0; c<6; ++c) {
        int previous=-100;
        for (int x=6; x>=-6; --x) {
            int y=maximum[c][x+6];
            if (y>previous) { groups[c].emplace_back(x,y); previous=y; }
        }
    }
    return intervals(groups,{0,1},{1,1},factor,false);
}

struct Hash {
    size_t operator()(const Signature& s) const {
        size_t h=1469598103934665603ULL;
        for (auto c : s) h=(h^c)*1099511628211ULL;
        return h;
    }
};
struct Cache {
    std::mutex lock;
    std::unordered_map<Signature,bool,Hash> map;
};

int main(int argc, char** argv) {
    if (argc!=2) return 1;
    fs::path dir=argv[1];
    auto scores=read<int8_t>(dir/"products",size_t(roots)*roots);
    auto gram=read<int8_t>(dir/"root_products",size_t(roots)*roots);
    auto incident=read<int>(dir/"incident",2*roots);
    auto strongest=read<int>(dir/"strongest",288*roots);
    auto colors=read<int>(dir/"colors",roots);
    std::array<Cache,1024> cache;
    int64_t checked=0, examined=0, axial=0, exceptions=0;
    #pragma omp parallel for schedule(dynamic,4) reduction(+:checked,examined,axial,exceptions)
    for (int i=0; i<roots; ++i) {
        const auto* a=scores.data()+size_t(i)*roots;
        for (int j=i+1; j<roots; ++j) {
            ++checked;
            const auto* b=scores.data()+size_t(j)*roots;
            int left=std::max(a[incident[2*j]],a[incident[2*j+1]]);
            int right=std::max(b[incident[2*i]],b[incident[2*i+1]]);
            if (left>=6 || right>=6 || 2*(6-right)>(4-right)*(6-left)) continue;
            Signature sig{};
            for (int c=0; c<6; ++c) {
                int best=-30, low=9, high=-9;
                for (int root : {i,j}) for (int rank=0; rank<12; ++rank) {
                    int m=strongest[root*288+rank*6+c], sum=a[m]+b[m];
                    if (sum>best) { best=sum; low=high=a[m]; }
                    else if (sum==best) { low=std::min(low,int(a[m])); high=std::max(high,int(a[m])); }
                }
                sig[3*c]=best+16; sig[3*c+1]=low+8; sig[3*c+2]=high+8;
            }
            int product=gram[size_t(i)*roots+j];
            sig[18]=left+8; sig[19]=right+8; sig[20]=product+4;
            auto& shard=cache[Hash{}(sig)%cache.size()];
            bool possible;
            {
                std::lock_guard<std::mutex> guard(shard.lock);
                auto it=shard.map.find(sig);
                if (it!=shard.map.end()) possible=it->second;
                else { possible=coarse(sig); shard.map.emplace(sig,possible); }
            }
            if (!possible) continue;
            ++examined;
            Result result=full(a,b,colors,16-4*product);
            axial+=(result==axial_midpoint);
            exceptions+=(result==other);
        }
    }
    std::cout << "{\"chords\":" << checked << ",\"full_envelopes_checked\":" << examined
              << ",\"axial_midpoint_chords\":" << axial << ",\"other_feasible_chords\":" << exceptions << "}\n";
    return exceptions!=0;
}
