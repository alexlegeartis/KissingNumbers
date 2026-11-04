// Exact independent F2 ranks of combined coordinate projections.
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

static unsigned rank(const std::vector<uint32_t>& generator, uint32_t mask) {
    std::array<uint32_t,32> basis{};
    unsigned count=0;
    for (uint32_t row:generator) {
        uint32_t word=row&mask;
        while(word) {
            unsigned pivot=31-__builtin_clz(word);
            if(basis[pivot]) word^=basis[pivot];
            else {basis[pivot]=word; ++count; break;}
        }
    }
    return count;
}

int main() {
    unsigned k,n;
    if(!(std::cin>>k>>n) || k!=17 || n!=3584) return 1;
    std::vector<uint32_t> generator(k),supports(n);
    for(auto& value:generator) if(!(std::cin>>value)) return 2;
    for(auto& value:supports) {
        if(!(std::cin>>value) || __builtin_popcount(value)!=12) return 3;
        if(rank(generator,value)!=11) return 4;
    }
    std::array<uint64_t,18> histogram{};
    std::vector<std::array<unsigned,2>> four_intersections;
    for(unsigned i=0;i<n;++i)
        for(unsigned j=0;j<i;++j) {
            unsigned dimension=k-rank(generator,supports[i]|supports[j]);
            ++histogram[dimension];
            if(dimension==4)four_intersections.push_back({i,j});
        }
    std::cout<<"{";
    bool first=true;
    for(unsigned i=0;i<histogram.size();++i) if(histogram[i]) {
        if(!first)std::cout<<",";
        first=false;
        std::cout<<"\""<<i<<"\":"<<histogram[i];
    }
    std::cout<<"}\n";
    std::cout<<"[";
    for(unsigned i=0;i<four_intersections.size();++i) {
        if(i)std::cout<<",";
        std::cout<<"["<<four_intersections[i][0]<<","<<four_intersections[i][1]<<"]";
    }
    std::cout<<"]\n";
}
