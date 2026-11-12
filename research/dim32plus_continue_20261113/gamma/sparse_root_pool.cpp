// Enumerate a bounded pool: signed sums of at most four distinct basis vectors.
#include <array>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
int gram[72][72];
vector<array<int,4>> pool;
int ids[4];
void enumerate(int size,int depth,int start){
    if(depth<size){
        for(int i=start;i<=72-(size-depth);++i){ids[depth]=i;enumerate(size,depth+1,i+1);}
        return;
    }
    for(int bits=0;bits<(1<<(size-1));++bits){
        int signs[4]={1,1,1,1};
        for(int i=1;i<size;++i)signs[i]=(bits&(1<<(i-1)))?1:-1;
        int norm=0;
        for(int i=0;i<size;++i){
            norm+=gram[ids[i]][ids[i]];
            for(int j=0;j<i;++j)norm+=2*signs[i]*signs[j]*gram[ids[i]][ids[j]];
        }
        if(norm==8){
            array<int,4> a{0,0,0,0};
            for(int i=0;i<size;++i)a[i]=signs[i]*(ids[i]+1);
            pool.push_back(a);
        }
    }
}
int main(int argc,char**argv){
    if(argc!=3)return 2;
    ifstream in(argv[1]);for(auto &r:gram)for(int &x:r)in>>x;
    for(int size=1;size<=4;++size)enumerate(size,0,0);
    ofstream out(argv[2]);
    for(auto a:pool){for(int x:a)out<<x<<' ';out<<'\n';}
    cout<<"Antipodal representatives "<<pool.size()<<'\n';
}
