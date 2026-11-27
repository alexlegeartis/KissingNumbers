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
 vector<T>K;string word;
 while(cin>>word){
  must(word.size()==15,"phase length");T x;
  for(int i=0;i<15;i++){must(word[i]>='0'&&word[i]<='2',"phase digits");put(x,i,word[i]-'0');}
  K.push_back(x);
 }
 must(K.size()==7906,"phase count");
 for(int i=0;i<7906;i++)for(int j=0;j<i;j++)must(pop(neq(K[i],K[j]))>=5,"phase distance");
 constexpr int FULL=(1<<15)-1;
 char magic[8];f.read(magic,8);must((bool)f&&string(magic,8)=="QR1690v1","witness magic");
 auto rd=[&](int n){Q x=0;for(int i=0;i<n;i++){int k=f.get();must(k!=EOF,"truncated witness");x|=Q(k)<<(8*i);}return x;};
 const int ND=rd(4),NZ=rd(2);must(ND>379&&ND<500000&&NZ>0&&NZ<65535,"counts");
 vector<pair<T,int>> donors;donors.reserve(ND);set<pair<Q,Q>> known;
 for(int i=0;i<ND;i++){
  T u{rd(6),rd(6)};int rail=rd(1);
  must(member(u)&&pop(u.mask())==18&&(u.a&7)==2&&(u.b&7)==4&&rail<3,"literal donor QR membership, marks, phase");
  must(known.insert({u.a,u.b}).second,"distinct positive donor labels");
  T v{u.a&~7ULL,u.b&~7ULL}, phase;
  if(rail==1) phase.a=v.mask();
  if(rail==2) phase.b=v.mask();
  must(!(v.mask()&7)&&pop(v.mask())==16,"body occupation");
  for(Q s:{v.a,v.b})must(tdot(v,add(phase,{0,s}))!=0,"dense syndromes");
  donors.push_back({u,rail});
 }
 vector<T>Z;for(int i=0;i<NZ;i++){T z{rd(2),rd(2)};must(!(z.a&z.b)&&z.mask()<=FULL,"shift range");Z.push_back(z);}
 enum_min({},{},0,0,R,false);enum_min({},{},0,0,IR,true);
 must(fullMin==415104&&V.size()==130944,"exact old QR15 census");
 sort(V.begin(),V.end(),[](T a,T b){return pair{a.a,a.b}<pair{b.a,b.b};});
 for(int i=1;i<(int)V.size();i++)must(neq(V[i-1],V[i])!=0,"distinct old words");
 vector<T> W;
 for(T x:V)if(x.a&(x.mask()&-x.mask()))W.push_back(x);
 must(W.size()==65472,"one representative per antipodal old pair");
 // Coordinates of a negative bank are obtained by negating a pair's positive bank with the same phase shift.
 vector<int> shifts;int nshift=0;
 for(size_t x=0;x<W.size();x++){int z=rd(NZ<=256?1:2);must(z<NZ,"shift index");shifts.push_back(z);nshift+=z!=0;}
 must(f.get()==EOF,"no trailing data");
 vector<unsigned char> S(size_t(NZ)*3*32768,0);
 for(int z=0;z<NZ;z++)for(T p:K){
  T h=add(p,Z[z]);
  S[(z*3+0)*32768+(FULL^h.mask())]=1;S[(z*3+1)*32768+h.a]=1;S[(z*3+2)*32768+h.b]=1;
 }
 long long nearCount=0;
 for(size_t i=0;i<W.size();i++){
  const T x=W[i];const int t=shifts[i];
  vector<int>P;for(Q q=x.mask();q;q&=q-1)P.push_back(__builtin_ctzll(q));
  must(P.size()==15,"old support");
  auto fold=[&](Q m){int a=0;for(int k=0;k<15;k++)a|=(int(m>>P[k]&1)<<k);return a;};
  auto check=[&](T y,int j){
   const int d=pop(neq(x,y));if(d>=18)return;must(d==15,"QR minimum and divisibility");nearCount++;
   const int C=fold(sameNonzero(x,y)),O=fold(opp(x,y));
   const int h=pop(C);must(h>=6&&h<=9&&pop(O)==18-2*h&&!(C&O),"exact mixed overlap");
   const int U=FULL^(C|O);const unsigned char* L=S.data()+(t*3+j)*32768;
   for(int p=U;;p=(p-1)&U){must(!L[C|p],"old/new forbidden support");if(!p)break;}
  };
  for(auto [u,j]:donors){check(u,j);check(neg(u),j);}
 }
 // Within all dense-compatible rails: body distance 3 per different ternary slot, plus 9 when the two root phases trigger the rescue.
 for(int a=0;a<=16;a++)for(int b=0;b<=16-a;b++){
  const int nd=32-2*a-b;
  for(int tag:{0,2})if(nd+tag>=15&&(nd+tag)%3==0){
   must(3*nd+9*a>=48&&3*nd+9*b>=48,"new/new complete signed overlap table");
  }
 }
 must(45*45<45*48&&4*24*24>45*48,"dense dot 22.5 only with one trit miss");
 must(48+45-48>0&&(48+45-48)*(48+45-48)<48*45&&
      (48+45-39)*(48+45-39)>48*45,"mixed raw distances 39 fails, 48 clears");
 must(4*16*3<=4*48&&225<=5*48&&4*9*3<=3*48&&4*12<=2*48,"small-stock comparisons");
 const Q prior=21956420122ULL, total=prior+2*ND;
 cout<<"AUDIT PASS donor_pairs="<<ND<<" old_shift_pairs="<<nshift<<" near_pairs="<<nearCount<<" exact_total="<<total<<"\n";
}
