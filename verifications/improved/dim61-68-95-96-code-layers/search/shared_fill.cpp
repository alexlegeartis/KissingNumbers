#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>
#include <string>
using namespace std;
struct Rec{uint64_t a;array<uint32_t,8> old;};
struct Choice{int idx;uint64_t phase;array<int,4> pm;};
struct Eq{vector<pair<uint64_t,int>>rows;pair<uint64_t,int>reduce(uint64_t m,int rhs)const{for(auto[r,b]:rows)if(m&(r&-r)){m^=r;rhs^=b;}return{m,rhs};}bool add(uint64_t m,int rhs){auto[r,b]=reduce(m,rhs);if(!r)return !b;auto it=rows.begin();while(it!=rows.end()&&__builtin_ctzll(it->first)<__builtin_ctzll(r))it++;rows.insert(it,{r,b});return true;}uint64_t solve()const{uint64_t x=0;for(auto it=rows.rbegin();it!=rows.rend();it++)if((__builtin_popcountll(it->first&x)&1)!=it->second)x^=it->first&-it->first;return x;}};
bool compatible(const vector<uint64_t>&sp,uint64_t h){if(find(sp.begin(),sp.end(),h)!=sp.end())return true;if(sp.size()>=32)return false;for(auto x:sp){int w=__builtin_popcountll(x^h);if(w!=0&&w!=8&&w!=16)return false;}return true;}
void extend(vector<uint64_t>&sp,uint64_t h){if(find(sp.begin(),sp.end(),h)!=sp.end())return;int n=sp.size();for(int i=0;i<n;i++)sp.push_back(sp[i]^h);}
struct Shape{array<int,4>L,R;array<uint64_t,4>P,T;uint64_t O,Q;};
int main(int argc,char**argv){if(argc<5){cerr<<"usage: shared_fill seconds candidates.bin output_prefix old_supports.bin\n";return 2;}double budget=stod(argv[1]);string file=argv[2],prefix=argv[3];vector<uint64_t>old(30828);ifstream(argv[4],ios::binary).read((char*)old.data(),old.size()*8);ifstream f(file,ios::binary|ios::ate);size_t n=f.tellg()/sizeof(Rec);f.seekg(0);vector<Rec>pool(n);f.read((char*)pool.data(),n*sizeof(Rec));
vector<Shape>shapes(n);for(size_t i=0;i<n;i++){auto&r=pool[i];auto&s=shapes[i];array<uint64_t,8>E;for(int j=0;j<8;j++)E[j]=r.a&~old[r.old[j]];vector<int>L={0},R;for(int j=1;j<8;j++)if(E[0]&E[j])L.push_back(j);else R.push_back(j);assert(L.size()==4&&R.size()==4);s.O=0;for(int j:L)s.O|=E[j];s.Q=r.a^s.O;for(int j=0;j<4;j++){s.P[j]=s.O^E[L[j]];s.T[j]=s.Q^E[R[j]];s.L[j]=L[j];s.R[j]=R[j];}for(int j=0;j<4;j++)for(int k=j+1;k<4;k++){if(s.P[j]>s.P[k]){swap(s.P[j],s.P[k]);swap(s.L[j],s.L[k]);}if(s.T[j]>s.T[k]){swap(s.T[j],s.T[k]);swap(s.R[j],s.R[k]);}}}
auto hs=[&](int i,const array<int,4>&pm){auto&s=shapes[i];array<uint64_t,8>H;for(int j=0;j<4;j++){H[s.L[j]]=s.P[j]|(s.Q^s.T[pm[j]]);H[s.R[pm[j]]]=s.T[pm[j]]|(s.O^s.P[j]);}return H;};
vector<array<int,4>>perms;array<int,4>p={0,1,2,3};do{perms.push_back(p);}while(next_permutation(p.begin(),p.end()));
vector<Choice>cur;ifstream seed(prefix+"_best.txt");if(seed){int count;seed>>count;for(int j=0;j<count;j++){uint64_t a;Choice c;seed>>a>>c.phase>>c.idx;for(auto&v:c.pm)seed>>v;assert(pool[c.idx].a==a);cur.push_back(c);}}size_t best=cur.size();mt19937 rng(2026102817);auto start=chrono::steady_clock::now();int iteration=0;
while(chrono::duration<double>(chrono::steady_clock::now()-start).count()<budget){vector<vector<uint64_t>>sp(old.size());for(size_t j=0;j<old.size();j++)sp[j]={0,old[j]};vector<Eq>eq(old.size());vector<Choice>next;vector<uint64_t>selected;shuffle(cur.begin(),cur.end(),rng);size_t keep=cur.size();if(keep){int mode=rng()%100;size_t remove=mode<65?1+rng()%max(size_t(1),keep/30):mode<95?keep/10:keep/2;keep-=min(keep,remove);}if(iteration%1000==999)keep=0;
for(size_t k=0;k<keep;k++){auto c=cur[k];auto&r=pool[c.idx];auto H=hs(c.idx,c.pm);for(int j=0;j<8;j++){assert(compatible(sp[r.old[j]],H[j]));extend(sp[r.old[j]],H[j]);assert(eq[r.old[j]].add(H[j],1^(__builtin_popcountll(c.phase&H[j])&1)));}next.push_back(c);selected.push_back(r.a);}
vector<int>order(n);iota(order.begin(),order.end(),0);shuffle(order.begin(),order.end(),rng);shuffle(perms.begin(),perms.end(),rng);
for(int i:order){auto&r=pool[i];bool okay=true;for(auto a:selected)if(__builtin_popcountll(a&r.a)>8){okay=false;break;}if(!okay)continue;for(auto pm:perms){auto H=hs(i,pm);bool fit=true;for(int j=0;j<8;j++)if(!compatible(sp[r.old[j]],H[j])){fit=false;break;}if(!fit)continue;Eq phase;for(int j=0;j<8;j++){auto[rem,rhs]=eq[r.old[j]].reduce(H[j],0);if(!rem&&!phase.add(H[j],1^rhs)){fit=false;break;}}if(!fit)continue;uint64_t shift=phase.solve();for(int j=0;j<8;j++){extend(sp[r.old[j]],H[j]);assert(eq[r.old[j]].add(H[j],1^(__builtin_popcountll(shift&H[j])&1)));}next.push_back({i,shift,pm});selected.push_back(r.a);break;}}
if(next.size()>best){best=next.size();ofstream out(prefix+"_best.txt");out<<best<<"\n";for(auto c:next){out<<pool[c.idx].a<<" "<<c.phase<<" "<<c.idx;for(int x:c.pm)out<<" "<<x;out<<"\n";}ofstream shifts(prefix+"_old.txt");for(size_t j=0;j<old.size();j++)shifts<<j<<" "<<eq[j].solve()<<"\n";cout<<"iteration "<<iteration<<" best "<<best<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;}
if(next.size()+2>=cur.size()||rng()%100<3)cur=move(next);iteration++;}
cout<<"finished best "<<best<<" iterations "<<iteration<<endl;}
