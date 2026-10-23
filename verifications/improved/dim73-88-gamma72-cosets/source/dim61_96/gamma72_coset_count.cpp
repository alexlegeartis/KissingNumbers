// Exact XOR convolutions count norm-11 triples, retaining block costs <=5.
// Omitting longer block vectors reduces the count; the result is a lower bound.
#define main screen_main_unused
#include "gamma72_coset_screen.cpp"
#undef main
#include <cassert>
void fwht(vector<int64_t>&v,bool inverse=false){
    size_t n=v.size();
    for(size_t len=1;len<n;len*=2)for(size_t i=0;i<n;i+=2*len)for(size_t j=0;j<len;++j){
        int64_t x=v[i+j],y=v[i+j+len];
        if(inverse){if((x+y)%2 || (x-y)%2)throw runtime_error("nonintegral inverse");v[i+j]=(x+y)/2;v[i+j+len]=(x-y)/2;}
        else{v[i+j]=x+y;v[i+j+len]=x-y;}
    }
}
int main(int argc,char**argv){
    if(argc!=7&&argc!=8){cerr<<"map list1 list2 list3 candidates output-prefix [total-norm=11]\n";return 2;}
    int total_norm=argc==8?stoi(argv[7]):11;
    uint32_t mask;int base;maps(argv[1],mask,base);
    vector<uint32_t> grouped[2][6][4096],third[6];
    for(int block=0;block<3;++block){auto raw=readbin(argv[2+block]);
        for(uint32_t rec:raw){int cost=rec>>24;uint32_t z=rec&MASK;
            if(cost<2||cost>5)continue;
            if(block<2)grouped[block][cost][apply(nmap,z)].push_back(z);
            else third[cost].push_back(z);
        }
    }
    vector<int64_t> p(1<<24),d3(1<<24),transformed(1<<24,0);
    uint64_t total_triples=0;auto started=chrono::steady_clock::now();
    for(int sum=max(4,total_norm-5);sum<=min(10,total_norm-2);++sum){
        fill(p.begin(),p.end(),0);fill(d3.begin(),d3.end(),0);
        uint64_t pair_count=0;
        for(int c1=2;c1<=5;++c1){int c2=sum-c1;if(c2<2||c2>5)continue;
            for(int r=0;r<4096;++r)for(uint32_t z1:grouped[0][c1][r]){
                uint32_t az1=apply(amap,z1);
                for(uint32_t z2:grouped[1][c2][r]){++p[az1^z2];++pair_count;}
            }
        }
        for(uint32_t z3:third[total_norm-sum])++d3[z3];
        total_triples+=pair_count*third[total_norm-sum].size();
        fwht(p);fwht(d3);
        for(size_t i=0;i<p.size();++i)transformed[i]+=p[i]*d3[i];
        cout<<"pair norm="<<sum<<" pairs="<<pair_count<<" total_triples="<<total_triples
            <<" seconds="<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<endl;
    }
    fwht(transformed,true);
    auto candidates=readbin(argv[5]);
    vector<pair<int64_t,uint32_t>> scored;scored.reserve(candidates.size());
    int64_t minimum=INT64_MAX,maximum=0,total=0;uint64_t above=0;
    for(uint32_t h:candidates){int64_t count=transformed[h];
        if(count<0)throw runtime_error("negative count");
        scored.push_back({count,h});minimum=min(minimum,count);maximum=max(maximum,count);total+=count;
        if(count>=2108)++above;
    }
    sort(scored.rbegin(),scored.rend());
    ofstream out(string(argv[6])+".json");
    out<<"{\n  \"construction\": \"Gamma72 half-coset heads, block costs at most 5\",\n"
        <<"  \"total_norm_unscaled\": "<<total_norm<<",\n"
        <<"  \"candidates\": "<<candidates.size()<<",\n  \"minimum_count\": "<<minimum
        <<",\n  \"maximum_count\": "<<maximum<<",\n  \"total_count\": "<<total
        <<",\n  \"at_least_2108\": "<<above<<",\n  \"all_triple_products\": "<<total_triples
        <<",\n  \"top\": [";
    for(size_t i=0;i<min(size_t(32),scored.size());++i)out<<(i?", ":"")<<"{\"h\": "<<scored[i].second<<", \"count\": "<<scored[i].first<<"}";
    out<<"]\n}\n";
    ofstream counts(string(argv[6])+".counts.bin",ios::binary);counts.write((char*)transformed.data(),8*transformed.size());
    cout<<"candidate range="<<minimum<<".."<<maximum<<" above2108="<<above<<endl;
}
