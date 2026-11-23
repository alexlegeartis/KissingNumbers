#include<bits/stdc++.h>
#include<immintrin.h>
using namespace std;using U=uint64_t;using V=uint16_t;
struct Code{V row[32]={};unsigned n=2,s=0;Code(){row[1]=65535;}
 int ins(V a,int p){
  for(unsigned i=0;i<n;i++)if(row[i]==a)return ((s>>i)&1)==unsigned(p)?0:-1;
  if(n==32)return -1;
  for(unsigned i=0;i<n;i++){V x=row[i]^a;if(popcount(x)!=8)return -1;row[i+n]=x;}
  s|=(s^(p?((1u<<n)-1):0))<<n;n*=2;return 1;
 }
 void complete(){
   for(V a=1;n<32;){a=0;unsigned r=countr_zero(n);array<unsigned,32> cnt={};for(int i=0;i<16;i++){unsigned label=0;for(unsigned j=0;j<r;j++)label|=((row[1u<<j]>>i)&1u)<<j;if(cnt[label]++&1)a|=V(1u<<i);}assert(ins(a,0)==1);}
   for(int i=0;i<16;i++)for(int j=i+1;j<16;j++){int cnt=0;for(unsigned k=0;k<32;k++)cnt+=popcount(unsigned(row[k]&((1u<<i)|(1u<<j))))&1;assert(cnt==16);}
 }
};
struct Support{U b;unsigned t;Code h;};
int main(int argc,char**argv){
 vector<U> Q(1<<22),cb,cub;for(int i=0;i<64;i++){
  U z=0;for(int j=0;j<64;j++)if((j&i)==i)z|=1ull<<j;
  if(popcount(unsigned(i))<=2)cb.push_back(z);if(popcount(unsigned(i))==3)cub.push_back(z);
 }
 for(unsigned i=1;i<Q.size();i++)Q[i]=Q[i-1]^cb[countr_zero(i)];
 // Monomials on subsets are linearly independent (Boolean Mobius transform).
 auto eval=[&](unsigned z){U a=0;for(;z;z&=z-1)a^=cub[countr_zero(z)];return a;};
 vector<unsigned>C;ifstream fc("cubic80.txt");for(unsigned c;fc>>c;)C.push_back(c);assert(fc.eof()&&C.size()==80&&set<unsigned>(C.begin(),C.end()).size()==80);
 set<unsigned>D;for(unsigned a:C)for(unsigned b:C)D.insert(a^b);
 for(unsigned d:D){U z=eval(d);for(unsigned i=(d==0);i<Q.size();i++)assert(popcount(z^Q[i])>=16);}
 array<unsigned,7>S={0,822310,1025326,213640,27520,779738,1041562};
 array<unsigned,7>T={0,3,12,5,10,6,9};
 array<unsigned,7>ct={37996,3904,3776,3904,3520,3904,3520};
 vector<Support>B;
 for(int i=0;i<7;i++){
  vector<U> cs;for(unsigned c:C){U z=eval(c^S[i]);for(U q:Q)if(popcount(z^q)==16-popcount(T[i]))cs.push_back(z^q);}
  assert(cs.size()==ct[i]);if(!i)sort(cs.begin(),cs.end());for(U x:cs)B.push_back({x,T[i],Code()});
 }
 int n0=B.size();assert(n0==60524);
 map<pair<int,int>,pair<V,V>> guards;
 ifstream f("signed_witness.txt");
 for(int i,j,a,b;f>>i>>j>>a>>b;){
  if(j==-1){assert(i==int(B.size()));U v=U(uint32_t(a))|(U(uint32_t(b))<<32);assert(popcount(v)==16);B.push_back({v,0,Code()});}
  else if(j==-2){assert(i==int(B.size())-1&&(b==0||b==1)&&B[i].h.n==2&&B[i].h.s==0);B[i].h.s=2*b;}
  else{
   assert(i==int(B.size())-1&&0<=j&&j<i&&0<=a&&a<=65535&&0<=b&&b<=131071);
   int p=b>>16;b&=65535;
   assert(popcount(unsigned(a))==8&&popcount(unsigned(b))==8);
   assert(_pdep_u64(a,B[i].b)==_pdep_u64(b,B[j].b));assert(B[i].t==0&&B[j].t==0);
   assert(B[i].h.ins(a,1-p)>=0&&B[j].h.ins(b,p)>=0);
   assert(guards.insert({{i,j},{a,b}}).second);
  }
 }
 assert(f.eof());
 assert(int(B.size())==n0+2354);
 for(auto& z:B)z.h.complete();
 size_t checked=0;for(int i=0;i<(int)B.size();i++)for(int j=0;j<i;j++){
  auto&a=B[i];auto&b=B[j];int k=popcount(a.b&b.b)+popcount(a.t&b.t);
  if(k<=8)continue;
  assert(k<=10&&guards.count({i,j}));++checked;
 }
 assert(checked==guards.size());cout<<"Checked signed overlap edges "<<checked<<"\n";
 // All sparse vectors have norm16. Each B-support uses 2048 signs at min Hamming distance >=4.
 // Dense signs norm64; cross dot <=16. Signed affine 4-flats norm4; cross dot <=4.
 // Distinct affine 2-flats meet in <=2 points. All coordinate poles are compatible.
 const U dense=(U(C.size())<<22), aff=10416*16;
 set<U> F;for(int a=0;a<64;a++)for(int b=a+1;b<64;b++)for(int c=b+1;c<64;c++)F.insert((1ull<<a)|(1ull<<b)|(1ull<<c)|(1ull<<(a^b^c)));
 assert(F.size()==10416);
 for(int d=64;d<=68;d++){
  U sparse=0;for(auto&v:B)sparse+=v.t<(1u<<(d-64));
  cout<<"DIM "<<d<<" sparse "<<sparse<<" N "<<dense+aff+2*d+2048*sparse<<"\n";
 }
 // Dimension 63: coordinate equality retains half dense signs (independent generator columns),
 // half signed points on sparse supports containing both positions (their dual has no weight2),
 // all on supports containing neither position.
 U mx=0;int p0=0,p1=0;
 vector<U> Q0;for(auto&z:B)if(z.t==0)Q0.push_back(z.b);
 int si[64]={},fi[64]={},sc[64][64]={},ac[64][64]={};
 auto inc=[&](U u,int* a,int(*b)[64]){for(U x=u;x;x&=x-1){int i=countr_zero(x);a[i]++;for(U y=x&(x-1);y;y&=y-1)b[i][countr_zero(y)]++;}};
 for(U x:Q0)inc(x,si,sc);for(U x:F)inc(x,fi,ac);
 for(int p=0;p<64;p++)for(int q=p+1;q<64;q++){
   U n=dense/2+2048*(Q0.size()-si[p]-si[q]+sc[p][q])+1024*sc[p][q]+16*(F.size()-fi[p]-fi[q]+ac[p][q])+8*ac[p][q]+124;
   if(n>mx){mx=n;p0=p;p1=q;}
 }
 cout<<"DIM 63 N "<<mx<<" at "<<p0<<","<<p1<<"\n";

 // Independently replay the inherited three small BCH cells. Their only
 // nonzero new coordinate is +/-sqrt(17); all sparse old vectors have L1<=32
 // when rescaled to squared norm64 and have zero last coordinate in 65D.
 auto degree=[](U v){
  for(int a=0;a<6;a++){U m=0;for(int j=0;j<64;j++)if(!(j>>a&1))m|=1ull<<j;v^=(v&m)<<(1<<a);}
  int r=0;for(int j=0;j<64;j++)if(v>>j&1)r=max(r,popcount(unsigned(j)));return r;
 };
 ifstream fp("polar.txt");int layers;fp>>layers;assert(fp&&layers==3);
 vector<U>H={0};for(int j=0;j<16;j++){U t;fp>>t;assert(degree(t)<=2);int m=H.size();for(int k=0;k<m;k++)H.push_back(H[k]^t);}
 assert(H.size()==65536&&set<U>(H.begin(),H.end()).size()==65536);
 for(U u:H)if(u)assert(popcount(u)>=22);
 vector<U>R(3);for(U &u:R){fp>>u;assert(degree(u)<=3);}
 U s;fp>>s;assert(fp&&degree(s)<=2);
 for(int i=0;i<layers;i++){
  for(int a:C){U z=R[i]^eval(a);for(U q:Q)assert(popcount(z^q)>=14);}
  for(int j=0;j<layers;j++)for(U h:H){
   if(i!=j)assert(popcount(R[i]^R[j]^h)>=22);
   assert(popcount(R[i]^R[j]^h^s)>=8);
  }
 }
 assert(64-2*22+17<=40&&64-2*8-17<=40&&64-2*14<=36&&32<=36&&64*17<=36*36);
 cout<<"POLAR 65 ADD "<<2*layers*H.size()<<endl;
}
