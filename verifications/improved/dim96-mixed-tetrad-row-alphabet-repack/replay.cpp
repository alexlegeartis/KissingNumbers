#include <bits/stdc++.h>

using namespace std;
using Bits=bitset<4096>;
unsigned take(unsigned x,unsigned m){
  unsigned t=0,k=1;
  for(;m;m&=m-1,k<<=1)if(x&(m&-m))t|=k;
  return t;
}
int pc(unsigned x){return __builtin_popcount(x);}
int ct(unsigned x){return __builtin_ctz(x);}
vector<unsigned>D{0};
struct T{unsigned a,b,m,s; array<int,24> ord,tag;};
T prep(unsigned a,unsigned b,unsigned s){
 T t{a,b,a|b,s,{},{}};
 array<unsigned,3> A{a&~b,b&~a,a&b};
 int j=0;
 for(int q=0;q<3;q++)for(unsigned m=A[q];m;m&=m-1){
   int i=ct(m);t.ord[i]=j++;t.tag[i]=q;
 }
 assert(j==12);
 return t;
}
unsigned spread(unsigned x){unsigned y=0;for(;x;x&=x-1)y|=3u<<(2*ct(x));return y;}
struct Bucket { vector<uint16_t> start,ids; };
map<unsigned,Bucket> cache;
const Bucket& getp(unsigned m){
 auto o=cache.find(m);if(o!=cache.end())return o->second;
 int n=1<<(2*pc(m)); unsigned q=spread(m);
 Bucket b;b.start.assign(n+1,0);b.ids.resize(4096);
 for(unsigned x:D)b.start[take(x,q)+1]++;
 for(int j=1;j<=n;j++)b.start[j]+=b.start[j-1];
 auto next=b.start;
 for(int i=0;i<4096;i++)b.ids[next[take(D[i],q)]++]=i;
 return cache.emplace(m,move(b)).first->second;
}
int tab[3][3][4][4];
void table(){
 int H[3][4][4]{};int p[3][4]={{0,1,2,3},{0,2,1,3},{0,3,1,2}};
 for(int g=0;g<3;g++)for(int x=0;x<4;x++){H[g][x][p[g][x&2]]=1;H[g][x][p[g][(x&2)+1]]=(x&1)?1:-1;}
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int x=0;x<4;x++)for(int y=0;y<4;y++){
  for(int k=0;k<4;k++)tab[i][j][x][y]+=H[i][x][k]*H[j][y][k];
  assert(i==j?(tab[i][j][x][y]==2*(x==y)):abs(tab[i][j][x][y])==1);
 }
}
struct Signs{
 array<unsigned,12> col{};
 array<Bits,16>M;
};
Signs assign(const vector<pair<int,int>>&loc){
 unsigned used=0;for(auto [i,j]:loc)used|=1u<<j;
 array<unsigned,4>B{};
 auto red=[&](unsigned x){for(int k=3;k>=0;k--)if(x>>k&1)x^=B[k];return x;};
 for(int j=0;j<12;j++)if(!(used>>j&1)){
  unsigned z=red(j+1);if(z)B[31-__builtin_clz(z)]=z;
 }
 Signs out;
 for(auto[i,j]:loc)out.col[i]=red(j+1);
 unsigned S[4096]{};
 for(int x=0;x<4096;x++){
  if(x)S[x]=S[x&(x-1)]^out.col[ct(x)];
  out.M[S[x]].set(x);
 }
 return out;
}
struct Interface{
 unsigned mn=0,mo=0,xn=0,xo=0,shift;
 vector<pair<int,int>> shared;
 vector<array<int,4>> off;
 vector<Signs> mats;
 const Bucket *bucket;
 Interface(T const& X,T const& Y):shift(Y.s){
   unsigned q=X.m&Y.m,r=q&~((X.a^Y.a)|(X.b^Y.b));
   assert(pc(q)+pc(r)==16);
   for(unsigned m=r;m;m&=m-1){
     int i=ct(m);shared.emplace_back(X.ord[i],Y.ord[i]);
     mn|=1u<<X.ord[i];mo|=1u<<Y.ord[i];
   }
   sort(shared.begin(),shared.end());
   for(size_t j=1;j<shared.size();j++)assert(shared[j-1].second<shared[j].second);
   for(unsigned m=q&~r;m;m&=m-1){int i=ct(m);off.push_back({X.ord[i],Y.ord[i],X.tag[i],Y.tag[i]});}
   xn=spread(mn);xo=spread(mo);
   bucket=&getp(mo);
   for(int miss=-1;miss<(int)shared.size();miss++){
    vector<pair<int,int>> rest;
    for(int j=0;j<(int)shared.size();j++)if(j!=miss)rest.push_back(shared[j]);
    for(auto t:off)rest.emplace_back(t[0],t[1]);
    mats.push_back(assign(rest));
   }
 }
 void filter(unsigned x,Bits&S)const{
   unsigned v=take(x,xn)^take(shift,xo);
   for(int m=-1;m<(int)shared.size();m++){
    const auto& C=mats[m+1];
    unsigned Q[9]{},bl=0,n=0;
    if(m<0){
      for(auto t:off)Q[++n]=C.col[t[0]];
    }
    for(int h=1;h<=(m<0?1:3);h++){
     unsigned key=v^(m<0?0u:unsigned(h)<<(2*m));
     for(unsigned at=bucket->start[key];at<bucket->start[key+1];at++){
      unsigned y=D[bucket->ids[at]]^shift;
      unsigned val=0;
      for(auto t:off)if(tab[t[2]][t[3]][(x>>(2*t[0]))&3][(y>>(2*t[1]))&3]<0)val^=C.col[t[0]];
      for(unsigned p=0;p<=n;p++)bl|=1u<<(val^Q[p]);
     }
    }
    for(unsigned p=bl;p;p&=p-1)S &= ~C.M[ct(p)];
    if(S.none())return;
   }
 }
};


array<Bits,12> bits;
Bits flip(const Bits&S,int i){
 int f=1<<i;
 return ((S&bits[i])>>f) | ((S&~bits[i])<<f);
}
int colored(const Bits &S){
 int E[13]{},C[13];
 for(int i=0;i<13;i++)C[i]=-1;
 auto edge=[&](int i,int j){E[i]|=1<<j;E[j]|=1<<i;};
 for(int i=0;i<12;i++){
   Bits A=flip(S,i);
   if((A&S).any())edge(i,12);
   for(int j=0;j<i;j++)if((flip(A,j)&S).any())edge(i,j);
 }
 C[12]=0;
 for(int sh=0;sh<12;sh++){
  int i=-1,best=-1,L=0;
  for(int j=0;j<12;j++)if(C[j]<0){
    int B=0;for(unsigned p=E[j];p;p&=p-1)if(C[ct(p)]>=0)B|=1<<C[ct(p)];
    int v=pc(B)*100+pc(E[j]);if(v>best){best=v;i=j;L=B;}
  }
  C[i]=ct(unsigned(~L));
 }
 int cnt[16]{};
 for(size_t i=S._Find_first();i<4096;i=S._Find_next(i)){
   int y=0;for(unsigned j=i;j;j&=j-1)y^=C[ct(j)];cnt[y]++;
 }
 return *max_element(cnt,cnt+16);
}
int main(){
 string seed="210100100101";
 for(int i=0;i<12;i++){
   unsigned x=0;
   for(int j=0;j<12;j++)x|=unsigned(seed[j]-'0')<<(2*((i+j)%12));
   int n=D.size();for(int j=0;j<n;j++)D.push_back(x^D[j]);
 }
 assert(set<unsigned>(D.begin(),D.end()).size()==4096);
 int md=12;
 for(unsigned x:D)if(x)md=min(md,pc((x|(x>>1))&0x555555));
 assert(md==6);
 table();
 for(int j=0;j<12;j++)for(int x=0;x<4096;x++)if(x>>j&1)bits[j].set(x);
 array<Bits,16> SG;
 int sd=12;
 for(int x=0;x<4096;x++){
   int t=0;
   for(int j=0;j<12;j++)if(x>>j&1)t^=j+1;
   if(x&&!t)sd=min(sd,pc(x));
   SG[t].set(x);
 }
 assert(sd==3);for(auto &s:SG)assert(s.count()==256);
 int N;cin>>N;assert(N==2970);
 vector<T>P;
 for(int j=0;j<N;j++){unsigned a,b,s;cin>>a>>b>>s;P.push_back(prep(a,b,s));}
 int K;cin>>K;assert(K>0);
 long long all=0,axiscount=0,edges=0;
 for(int j=0;j<K;j++){
   unsigned a,b,s;long long expected;string p, m;
   cin>>a>>b>>p>>s>>m>>expected;
   assert(p.size()==12 && s<0x1000000);
   array<int,12>A;int used=0;
   for(int h=0;h<12;h++){
    A[h]=(p[h]<'a'?p[h]-'0':p[h]-'a'+10);
    assert(0<=A[h]&&A[h]<12&&!(used>>A[h]&1));used|=1<<A[h];
   }
   assert(m.size()==36);
   unsigned B[12][4]{};
   for(int h=0;h<12;h++){
    unsigned mask=1;
    for(int q=1;q<4;q++){
      unsigned t=m[3*h+q-1]-'0';
      assert(t>0 && t<4 && !(mask>>t&1));
      mask|=1u<<t;B[h][q]=t;
    }
   }
   T X=prep(a,b,0);
   vector<Interface>I;
   for(const auto &Y:P){
     unsigned q=X.m&Y.m,r=q&~((X.a^Y.a)|(X.b^Y.b));
     if(pc(q)+pc(r)<=12)continue;
     assert(pc(q)+pc(r)==16);
     I.emplace_back(X,Y);
   }
   assert(!I.empty());edges+=I.size();
   long long val=0;
   for(unsigned u:D){
    unsigned x=s;
    for(int h=0;h<12;h++)x^=B[h][(u>>(2*A[h]))&3]<<(2*h);
    Bits S;S.set();
    for(const auto& T:I){T.filter(x,S);if(S.none())break;}
    int best=0;
    for(const auto& C:SG)best=max(best,int((S&C).count()));
    int t=S.any()?colored(S):0;
    val+=max(best,t);
    axiscount+=(best!=0);
   }
   assert(val==expected);
   all+=val;
 }
 cout<<K<<' '<<all<<' '<<axiscount<<' '<<edges<<'\n';
}
