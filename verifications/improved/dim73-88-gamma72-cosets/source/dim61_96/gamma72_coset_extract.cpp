// Extract individual norm-11 triples for a single screened coset.
#define main screen_main_unused
#include "gamma72_coset_screen.cpp"
#undef main
#include <cassert>
int main(int argc,char**argv){
    if(argc!=7){cerr<<"map list1 list2 list3 h output-prefix\n";return 2;}
    uint32_t mask;int base;maps(argv[1],mask,base);
    vector<uint32_t> raw[3];for(int i=0;i<3;++i)raw[i]=readbin(argv[2+i]);
    vector<uint32_t> grouped[2][6][4096];
    vector<int32_t> third(1<<24,-1);
    for(int block=0;block<3;++block)for(uint32_t index=0;index<raw[block].size();++index){
        uint32_t rec=raw[block][index],z=rec&MASK;int cost=rec>>24;
        if(cost<2||cost>5)throw runtime_error("block cost");
        if(block<2)grouped[block][cost][apply(nmap,z)].push_back(index);
        else{if(third[z]!=-1)throw runtime_error("duplicate key needs multimap");third[z]=index;}
    }
    uint32_t h=stoul(argv[5]);vector<array<uint32_t,3>> triples;
    for(int c1=2;c1<=5;++c1)for(int c2=2;c2<=5;++c2){
        int c3=11-c1-c2;if(c3<2||c3>5)continue;
        for(int r=0;r<4096;++r)for(uint32_t i1:grouped[0][c1][r]){
            uint32_t az1=apply(amap,raw[0][i1]&MASK);
            for(uint32_t i2:grouped[1][c2][r]){
                int32_t i3=third[h^az1^(raw[1][i2]&MASK)];
                if(i3>=0 && (raw[2][i3]>>24)==c3)triples.push_back({i1,i2,uint32_t(i3)});
            }
        }
    }
    ofstream out(string(argv[6])+".triples.bin",ios::binary);
    for(auto&triple:triples)out.write((char*)triple.data(),12);
    for(int block=0;block<3;++block){vector<uint32_t> indices;for(auto&t:triples)indices.push_back(t[block]);
        sort(indices.begin(),indices.end());indices.erase(unique(indices.begin(),indices.end()),indices.end());
        ofstream wanted(string(argv[6])+".wanted"+to_string(block)+".bin",ios::binary);
        wanted.write((char*)indices.data(),4*indices.size());cout<<"block "<<block<<" selected="<<indices.size()<<endl;
    }
    cout<<"triples="<<triples.size()<<endl;
    return 0;
}
