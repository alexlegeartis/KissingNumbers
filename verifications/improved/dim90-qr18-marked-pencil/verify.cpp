#include <bits/stdc++.h>
using namespace std;
using Q=uint64_t;
struct T{Q a=0,b=0;Q mask()const{return a|b;}};
constexpr Q E24=(1u<<24)-1,M48=(1ull<<48)-1;
inline T add(T x,T y){return {(x.a&~(y.a|y.b))|(~(x.a|x.b)&y.a)|(x.b&y.b),(x.b&~(y.a|y.b))|(~(x.a|x.b)&y.b)|(x.a&y.a)};}
inline T neg(T x){return {x.b,x.a};}
inline int pop(Q a){return __builtin_popcountll(a);}
inline int digit(T x,int i){return (x.a>>i&1)+2*(x.b>>i&1);}
inline void put(T& x,int i,int v){if(v==1)x.a|=1ull<<i;if(v==2)x.b|=1ull<<i;}
inline int mod(int a){return (a+30)%3;}
inline Q neq(T x,T y){return (x.a^y.a)|(x.b^y.b);}
inline Q sameNonzero(T x,T y){return (x.a&y.a)|(x.b&y.b);}
inline Q opp(T x,T y){return (x.a&y.b)|(x.b&y.a);}
inline int tdot(T x,T y){return mod(pop(sameNonzero(x,y))-pop(opp(x,y)));}
inline T cat(T x,T y){return{x.a|(y.a<<24),x.b|(y.b<<24)};}
const T R[24]={{9573376,295655},{10758144,591310},{13127680,1182620},{13511,2261272},{27022,4522544},{9453944,790659},{10519280,1581318},{12649952,3162636},{265511,7344856},{10094760,6598929},{2230066,14282885},{4233218,2188681},{467169,13651220},{1623840,2230285},{3247680,4460570},{15796260,964179},{7374345,616512},{7805557,8404994},{15213065,1310946},{4506161,2756866},{1638529,12759138},{2364868,1303073},{4729736,2606146},{147827,13175296}}; const T IR[24]={{2141729,10424584},{6209635,8391960},{13511,2261272},{27022,4522544},{54044,9045088},{3407889,9091040},{5603587,10637544},{265511,7344856},{531022,14689712},{10518933,5203008},{4233218,2188681},{8466436,4377362},{1623840,2230285},{3247680,4460570},{6495360,8921140},{7374345,616512},{14748690,1233024},{7905292,8525313},{4506161,2756866},{9012322,5513732},{2364868,1303073},{4729736,2606146},{9459472,5212292},{9573376,295655}};

static void must(bool ok, const char* s){if(!ok){cerr<<"ERROR: "<<s<<"\n";exit(2);}}
T body(T x){T r;for(int i=0;i<24;i++){int z=digit(x,i);if(z)r=add(r,z==1?R[i]:neg(R[i]));}return r;}
// Explicit parentheses to avoid bitwise/operator-precedence ambiguity.
bool member(T x){T z=body(x);return !(x.a&x.b)&&((x.a>>24)==z.a)&&((x.b>>24)==z.b)&&(x.mask()<=M48);}
vector<T> V;
int fullMin=0;
void enum_min(T x,T y,int next,int w,const T A[24],bool rev){
 if(w+pop(y.mask())==15){fullMin++;T z=rev?cat(y,x):cat(x,y);if(!(z.mask()&7))V.push_back(z);}
 if(w==7)return;
 for(int i=next;i<24;i++){
  const Q p=1ull<<i;
  enum_min({x.a|p,x.b},add(y,A[i]),i+1,w+1,A,rev);
  enum_min({x.a,x.b|p},add(y,neg(A[i])),i+1,w+1,A,rev);
 }
}
struct Z{T u,v,p;};
constexpr int BLOCK=(7906+63)/64;
Q bitset0[15][3][BLOCK];
int trid[15][15][15];
vector<Q> tri;
inline Q samedigit(T p,T q){return~neq(p,q);}
inline int twiceDotDiv3(T a,T p,T x,T q){
 Q A=sameNonzero(a,x), B=opp(a,x), D=samedigit(p,q);
 return 3*(pop(A&D)-pop(B&D))+pop(B)-pop(A);
}

void canonical_QR(){
 auto poly=[](initializer_list<pair<int,int>> a){vector<int> z(24);for(auto [i,j]:a)z[i]=mod(j);return z;};
 const auto g=poly({{23,1},{20,1},{18,-1},{17,1},{15,-1},{12,1},{10,1},{9,-1},{7,-1},{6,-1},{5,-1},{2,-1},{1,-1},{0,-1}});
 const auto h=poly({{23,1},{22,1},{21,1},{18,1},{17,1},{16,1},{14,1},{13,-1},{11,-1},{8,1},{6,-1},{5,1},{3,-1},{0,-1}});
 auto mul=[](const vector<int>& a,const vector<int>& b){vector<int> c(a.size()+b.size()-1);for(unsigned i=0;i<a.size();i++)for(unsigned j=0;j<b.size();j++)c[i+j]=mod(c[i+j]+a[i]*b[j]);return c;};
 auto f=mul({2,1},g), z=mul(f,h);
 must(z.size()==48&&z[0]==2&&z[47]==1,"QR polynomial product degree");
 for(int i=1;i<47;i++)must(z[i]==0,"QR factor (x-1)gh = x^47-1");
 must(mod(accumulate(g.begin(),g.end(),0))&&mod(accumulate(h.begin(),h.end(),0)),"nontrivial factors");
 int order=1,r=3;while(r!=1){++order;r=3*r%47;}
 must(order==23,"order of 3 mod 47 gives the two residue coset factors");
 int C[24][48]={};
 for(int i=0;i<23;i++)for(int j=0;j<25;j++)C[i][i+j]=f[j];
 for(int j=0;j<48;j++)C[23][j]=1;
 for(int i=0;i<24;i++){
  int j=i;while(j<24&&!C[j][i])j++;
  must(j<24,"canonical QR systematic pivot");
  if(j!=i)for(int k=0;k<48;k++)swap(C[i][k],C[j][k]);
  if(C[i][i]==2)for(int k=0;k<48;k++)C[i][k]=mod(C[i][k]*2);
  for(int j=0;j<24;j++)if(j!=i&&C[j][i]){
   int s=C[j][i];for(int k=0;k<48;k++)C[j][k]=mod(C[j][k]-s*C[i][k]);
  }
 }
 for(int i=0;i<24;i++)for(int j=0;j<24;j++){
  must(C[i][j]==(i==j)&&C[i][24+j]==digit(R[i],j),"literal R equals canonical QR with unchanged coordinate order");
  must(digit(IR[i],j)==mod(-digit(R[j],i)),"IR = -R^T");
  must(tdot(R[i],R[j])==mod(-(i==j)),"R R^T = -I");
  int x=0;for(int k=0;k<24;k++)x=mod(x+digit(R[i],k)*digit(IR[k],j));
  must(x==(i==j),"IR two-half inverse");
 }
}

int main(int argc,char**argv){
 canonical_QR();
 must(argc==2,"usage: audit witness.bin < phase.txt");
 ifstream f(argv[1],ios::binary);must((bool)f,"witness readable");
 vector<T>K;string w;
 while(cin>>w){
  must(w.size()==15,"phase length");T x;for(int i=0;i<15;i++){must(w[i]>='0'&&w[i]<='2',"phase digits"); put(x,i,w[i]-'0');}K.push_back(x);
 }
 must(K.size()==7906,"phase count");
 for(int i=0;i<7906;i++)for(int j=0;j<i;j++)must(pop(neq(K[i],K[j]))>=5,"phase distance");

 char magic[8];f.read(magic,8);must((bool)f&&string(magic,8)=="QR1890v1","witness magic");
 auto rd=[&](int n){Q x=0;for(int i=0;i<n;i++){int k=f.get();must(k!=EOF,"truncated witness");x|=Q(k)<<(8*i);}return x;};
 int expN=(int)rd(4);must(expN>=1&&expN<=100000,"witness count");Q expTotal=21956420122ULL+expN;
 vector<Z>N;N.reserve(expN);
 for(int i=0;i<expN;i++){Z z;z.u={rd(6),rd(6)};z.v={rd(6),rd(6)};z.p={rd(6),rd(6)};N.push_back(z);}
 vector<T> shifts;int nonzero=0;
 for(int i=0;i<130944;i++){
  int n=(int)rd(3);T p;for(int j=14;j>=0;j--){put(p,j,n%3);n/=3;}
  must(n==0,"shift digits");nonzero+=p.mask()!=0;shifts.push_back(p);
 }
 must(f.get()==EOF,"no trailing data");
 for(int i=0;i<(int)N.size();i++){
  auto z=N[i];must(member(z.u)&&member(z.v),"new QR membership");
  must(!(z.u.mask()&7)&&pop(z.u.mask())==18&&(z.v.mask()&7)==7&&pop(z.v.mask())==15,"new QR weight");
  must(!(z.v.mask()&~(z.u.mask()|7)),"pencil containment");
  must(!(z.p.a&z.p.b)&&!(z.p.mask()&~z.u.mask()),"phase support");
  for(Q s:{z.u.a,z.u.b}){
   T p=add(z.p,{0,s});must(tdot(z.u,p)==0&&tdot(z.v,p)!=0,"both dense syndromes");
  }
  for(int j=0;j<i;j++)must(twiceDotDiv3(z.u,z.p,N[j].u,N[j].p)<=18,"new/new dot");
 }
 enum_min({},{},0,0,R,false);enum_min({},{},0,0,IR,true);
 must(fullMin==415104&&V.size()==130944,"exact parent and shortened population");
 sort(V.begin(),V.end(),[](T a,T b){return pair{a.a,a.b}<pair{b.a,b.b};});
 for(int i=1;i<(int)V.size();i++)must(neq(V[i-1],V[i])!=0,"old uniqueness");
 for(int i=0;i<7906;i++)for(int j=0;j<15;j++)bitset0[j][digit(K[i],j)][i>>6]|=1ull<<(i&63);
 int nt=0;
 for(int a=0;a<15;a++)for(int b=a+1;b<15;b++)for(int c=b+1;c<15;c++)trid[a][b][c]=nt++;
 tri.resize(nt*27*BLOCK);
 for(int a=0;a<15;a++)for(int b=a+1;b<15;b++)for(int c=b+1;c<15;c++)
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++){
   auto* p=tri.data()+((trid[a][b][c]*27+i*9+j*3+k)*BLOCK);
   for(int l=0;l<BLOCK;l++)p[l]=bitset0[a][i][l]&bitset0[b][j][l]&bitset0[c][k][l];
  }
 long long nearCount=0;
 for(int ix=0;ix<(int)V.size();ix++){
  T x=V[ix],t=shifts[ix];
  vector<int>P;for(Q h=x.mask();h;h&=h-1)P.push_back(__builtin_ctzll(h));
  for(auto z:N){
   int d=pop(neq(x,z.u));if(d>=18)continue;
   must(d==15,"near nonmultiple of 3");nearCount++;
   Q eq=sameNonzero(x,z.u), no=opp(x,z.u);
   int A[15],B[15],a=0,b=0;
   int D[15];
   for(int i=0;i<15;i++){
    D[i]=mod(digit(z.p,P[i])-digit(t,i));
    if(eq>>P[i]&1)A[a++]=i;
    if(no>>P[i]&1)B[b++]=i;
   }
   must(a>=6&&a<=9&&b<=6,"rank-two intersection profile");
   auto* F=tri.data()+((trid[A[0]][A[1]][A[2]]*27+D[A[0]]*9+D[A[1]]*3+D[A[2]])*BLOCK);
   auto* G=tri.data()+((trid[A[3]][A[4]][A[5]]*27+D[A[3]]*9+D[A[4]]*3+D[A[5]])*BLOCK);
   for(int j=0;j<BLOCK;j++){
    Q p=F[j]&G[j];if(!p)continue;
    for(int k=6;k<a&&p;k++)p&=bitset0[A[k]][D[A[k]]][j];
    for(int k=0;k<b&&p;k++)p&=~bitset0[B[k]][D[B[k]]][j];
    must(!p,"old/new near phases");
   }
  }
 }
 cout<<"AUDIT PASS count="<<N.size()<<" old_shifts="<<nonzero<<" near="<<nearCount<<" exact_total="<<expTotal<<"\n";
}

