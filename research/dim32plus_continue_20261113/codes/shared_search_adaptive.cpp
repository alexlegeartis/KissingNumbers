#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
using U=uint64_t;
int wt(U x){return __builtin_popcountll(x);}
struct Rec{U a;array<uint32_t,8> old;};
struct Choice{int idx;U phase;array<int,4> pm;U extra=0;};
struct Eq {
 vector<pair<U,int>>rows;
 pair<U,int>reduce(U m,int rhs)const{for(auto[r,b]:rows)if(m&(r&-r)){m^=r;rhs^=b;}return{m,rhs};}
 bool add(U m,int rhs){auto[r,b]=reduce(m,rhs);if(!r)return !b;auto it=rows.begin();while(it!=rows.end()&&__builtin_ctzll(it->first)<__builtin_ctzll(r))it++;rows.insert(it,{r,b});return true;}
 U solve()const{U x=0;for(auto it=rows.rbegin();it!=rows.rend();it++)if((wt(it->first&x)&1)!=it->second)x^=it->first&-it->first;return x;}
};
bool compatible(const vector<U>&sp,U h){if(find(sp.begin(),sp.end(),h)!=sp.end())return true;if(sp.size()>=32)return false;for(U x:sp){int w=wt(x^h);if(w!=0&&w!=8&&w!=16)return false;}return true;}
void extend(vector<U>&sp,U h){if(find(sp.begin(),sp.end(),h)!=sp.end())return;int n=sp.size();for(int i=0;i<n;i++)sp.push_back(sp[i]^h);}
struct Shape{array<int,4>L,R;array<U,4>P,T;U O,Q;};
struct Selected{U a,phase;vector<U>span;};

int main(int argc,char**argv){
 assert(argc==9);
 double seconds=stod(argv[1]);string prefix=argv[3];U removed=stoull(argv[5]);
 mt19937 rng(stoul(argv[6]));bool odd=stoi(argv[7]),repair=stoi(argv[8]);
 vector<U>old(30828);ifstream of(argv[4],ios::binary);assert(of.read((char*)old.data(),old.size()*8));
 ifstream f(argv[2],ios::binary|ios::ate);assert(f);size_t bytes=f.tellg();assert(bytes%sizeof(Rec)==0);f.seekg(0);
 vector<Rec>pool(bytes/sizeof(Rec));assert(f.read((char*)pool.data(),bytes));size_t n=pool.size();
 for(auto&r:pool)assert(!(r.a&removed));
 vector<Shape>shapes(n);
 for(size_t i=0;i<n;i++){
  auto&r=pool[i];auto&s=shapes[i];array<U,8>E;for(int j=0;j<8;j++){assert(r.old[j]<old.size());E[j]=r.a&~old[r.old[j]];}
  vector<int>L={0},R;for(int j=1;j<8;j++)if(E[0]&E[j])L.push_back(j);else R.push_back(j);
  assert(L.size()==4&&R.size()==4);s.O=0;for(int j:L)s.O|=E[j];s.Q=r.a^s.O;
  for(int j=0;j<4;j++){s.P[j]=s.O^E[L[j]];s.T[j]=s.Q^E[R[j]];s.L[j]=L[j];s.R[j]=R[j];}
  for(int j=0;j<4;j++)for(int k=j+1;k<4;k++){
   if(s.P[j]>s.P[k]){swap(s.P[j],s.P[k]);swap(s.L[j],s.L[k]);}
   if(s.T[j]>s.T[k]){swap(s.T[j],s.T[k]);swap(s.R[j],s.R[k]);}
  }
 }
 auto checks=[&](int i,const array<int,4>&pm){auto&s=shapes[i];array<U,8>H;for(int j=0;j<4;j++){H[s.L[j]]=s.P[j]|(s.Q^s.T[pm[j]]);H[s.R[pm[j]]]=s.T[pm[j]]|(s.O^s.P[j]);}return H;};
 auto basic=[&](U a,const array<U,8>&H){vector<U>sp={0,a};for(U h:H){assert(compatible(sp,h));extend(sp,h);}assert(sp.size()==16);return sp;};
 auto complete=[&](U a,vector<U>&sp){
  assert(sp.size()==16);vector<pair<unsigned,int>>columns;
  for(int p=0;p<64;p++)if(a>>p&1){unsigned code=0;for(size_t j=0;j<sp.size();j++)code|=((sp[j]>>p)&1)<<j;columns.push_back({code,p});}
  sort(columns.begin(),columns.end());assert(columns.size()==16);U h=0;
  for(int j=0;j<16;j+=2){assert(columns[j].first==columns[j+1].first);if(j)assert(columns[j].first!=columns[j-1].first);h|=U(1)<<columns[j+(rng()&1)].second;}
  assert(compatible(sp,h));extend(sp,h);assert(sp.size()==32);return h;
 };
 vector<array<int,4>>perms;array<int,4>p={0,1,2,3};do{perms.push_back(p);}while(next_permutation(p.begin(),p.end()));
 vector<Choice>cur;ifstream seed(prefix+"_best.txt");
 if(seed){int count;seed>>count;string line;getline(seed,line);for(int j=0;j<count;j++){getline(seed,line);istringstream row(line);U a;Choice c;row>>a>>c.phase>>c.idx;for(auto&x:c.pm)row>>x;assert(row);row>>c.extra;assert(c.idx>=0&&size_t(c.idx)<n&&pool[c.idx].a==a);assert(!odd||(wt(a&c.phase)&1));auto sp=basic(a,checks(c.idx,c.pm));if(c.extra){assert(compatible(sp,c.extra));extend(sp,c.extra);assert(sp.size()==32);}else c.extra=complete(a,sp);cur.push_back(c);}}
 size_t best=cur.size();int iteration=0;auto start=chrono::steady_clock::now();
 cout<<"configuration candidates "<<n<<" removed "<<removed<<" odd "<<odd<<" new_pair_repairs "<<repair<<" seed_size "<<best<<endl;
 while(chrono::duration<double>(chrono::steady_clock::now()-start).count()<seconds){
  vector<vector<U>>sp(old.size());for(size_t j=0;j<old.size();j++)sp[j]={0,old[j]};vector<Eq>eq(old.size());
  if(odd)for(size_t j=0;j<old.size();j++)if(!(old[j]&removed))assert(eq[j].add(old[j],1));
  vector<Choice>next;vector<Selected>selected;shuffle(cur.begin(),cur.end(),rng);size_t keep=cur.size();
  if(keep){int mode=rng()%100;size_t erase=mode<65?1+rng()%max(size_t(1),keep/30):mode<95?keep/10:keep/2;keep-=min(keep,erase);}if(iteration%1000==999)keep=0;
  for(size_t k=0;k<keep;k++){auto c=cur[k];auto&r=pool[c.idx];auto H=checks(c.idx,c.pm);
   for(int j=0;j<8;j++){if(old[r.old[j]]&removed)continue;assert(compatible(sp[r.old[j]],H[j]));extend(sp[r.old[j]],H[j]);assert(eq[r.old[j]].add(H[j],1^(wt(c.phase&H[j])&1)));}
   auto own=basic(r.a,H);extend(own,c.extra);assert(own.size()==32);sort(own.begin(),own.end());next.push_back(c);selected.push_back({r.a,c.phase,move(own)});
  }
  vector<int>order(n);iota(order.begin(),order.end(),0);shuffle(order.begin(),order.end(),rng);shuffle(perms.begin(),perms.end(),rng);
  for(int i:order){auto&r=pool[i];bool fit=true;vector<int>ten;
   for(size_t j=0;j<selected.size();j++){int w=wt(r.a&selected[j].a);if(w>8){if(w!=10||!repair){fit=false;break;}ten.push_back(j);}}if(!fit)continue;
   for(auto pm:perms){auto H=checks(i,pm);fit=true;
    for(int j=0;j<8;j++)if(!(old[r.old[j]]&removed)&&!compatible(sp[r.old[j]],H[j])){fit=false;break;}if(!fit)continue;
    Eq phase;if(odd)assert(phase.add(r.a,1));
    for(int j=0;j<8;j++){if(old[r.old[j]]&removed)continue;auto[rem,rhs]=eq[r.old[j]].reduce(H[j],0);if(!rem&&!phase.add(H[j],1^rhs)){fit=false;break;}}if(!fit)continue;
    auto own=basic(r.a,H);U extra=0;
    for(int j:ten){auto&other=selected[j];U overlap=r.a&other.a,common=0;
     for(U h:other.span)if(wt(h)==8&&(h&overlap)==h){assert(!common);common=h;}
     if(!common||!compatible(own,common)||!phase.add(common,1^(wt(other.phase&common)&1))){fit=false;break;}
     if(find(own.begin(),own.end(),common)==own.end()){extra=common;extend(own,common);}
    }if(!fit)continue;
    if(own.size()==16)extra=complete(r.a,own);assert(extra&&own.size()==32);sort(own.begin(),own.end());U shift=phase.solve();
    for(int j=0;j<8;j++){if(old[r.old[j]]&removed)continue;extend(sp[r.old[j]],H[j]);assert(eq[r.old[j]].add(H[j],1^(wt(shift&H[j])&1)));}
    next.push_back({i,shift,pm,extra});selected.push_back({r.a,shift,move(own)});break;
   }
  }
  if(next.size()>best){best=next.size();ofstream out(prefix+"_best.txt");out<<best<<"\n";for(auto c:next){out<<pool[c.idx].a<<" "<<c.phase<<" "<<c.idx;for(int x:c.pm)out<<" "<<x;out<<" "<<c.extra<<"\n";}cout<<"iteration "<<iteration<<" best "<<best<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;}
  if(next.size()+2>=cur.size()||rng()%100<3)cur=move(next);iteration++;
 }
 cout<<"finished best "<<best<<" iterations "<<iteration<<endl;
}
