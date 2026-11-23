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
#include "atlas46.inc"
// Reconstruction of all published sparse masks and their shared overlap guards.
int main(){
 ifstream A("assignment.txt");
 int nh;A>>nh;assert(nh==46);
 U used=0;
 for(int i=0;i<32;i++){
   int j,x,y;A>>j>>x>>y;
   assert(0<=j&&j<80&&selected_sites[j]==0&&0<=x&&x<64&&0<=y&&y<64&&x!=y);
   U w=(U(1)<<x)|(U(1)<<y);
   assert(!(used&w));used|=w;selected_sites[j]=w;
 }
 assert(used==~U(0)); for(auto &m:selected_sites)m=0;
 for(int j;A>>j;){assert(0<=j&&j<126);assert(selected_heads.insert(j).second);}
 assert(A.eof()&&int(selected_heads.size())==nh);
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
 int n1=B.size();
 ifstream z("tagged_witness.txt");
 auto pdep=[&](int a,const Support& b)->pair<U,U>{
  return {_pdep_u64(a,b.b),_pdep_u64(unsigned(a)>>popcount(b.b),b.t)};
 };
 for(int i,j,tag,a,b;z>>i>>j>>tag>>a>>b;){
  if(j==-1){
   assert(i==int(B.size())&&0<=tag&&tag<16);U v=U(uint32_t(a))|(U(uint32_t(b))<<32);
   assert(popcount(v)+popcount(unsigned(tag))==16);B.push_back({v,unsigned(tag),Code()});
  }else if(j==-2){
   assert(i==int(B.size())-1&&(b==0||b==1)&&B[i].h.n==2&&B[i].h.s==0);B[i].h.s=2*b;
  }else{
   assert(i==int(B.size())-1&&0<=j&&j<i&&0<=a&&a<=65535&&0<=b&&b<=131071);
   int p=b>>16;b&=65535;
   assert(popcount(unsigned(a))==8&&popcount(unsigned(b))==8);
   assert(pdep(a,B[i])==pdep(b,B[j]));
   assert(B[i].h.ins(a,1-p)>=0&&B[j].h.ins(b,p)>=0);
   assert(guards.insert({{i,j},{a,b}}).second);
  }
 }
 assert(z.eof());assert(int(B.size())-n1==4981);
 cout<<"New tagged signed supports "<<int(B.size())-n1<<"\n";


 int base=B.size();
 ifstream m("joint_witness.txt");
 for(int i,j,tag,a,b;m>>i>>j>>tag>>a>>b;){
  if(j==-1){
   assert(i==int(B.size())&&0<=tag&&tag<16);U v=U(uint32_t(a))|(U(uint32_t(b))<<32);
   assert(popcount(v)+popcount(unsigned(tag))==16);B.push_back({v,unsigned(tag),Code()});
  }else if(j==-2){
   assert(i==int(B.size())-1&&(b==0||b==1)&&B[i].h.n==2&&B[i].h.s==0);B[i].h.s=2*b;
  }else if(j==-3){
   assert(i==int(B.size())-1&&0<=a&&a<=65535&&(b==0||b==1)&&popcount(unsigned(a))==8);
   assert(B[i].h.ins(a,b)>=0);
  }else{
   assert(i==int(B.size())-1&&0<=j&&j<i&&0<=a&&a<=65535&&0<=b&&b<=131071);
   int p=b>>16;b&=65535;
   assert(popcount(unsigned(a))==8&&popcount(unsigned(b))==8);
   assert(pdep(a,B[i])==pdep(b,B[j]));
   assert(B[i].h.ins(a,1-p)>=0&&B[j].h.ins(b,p)>=0);
   assert(guards.insert({{i,j},{a,b}}).second);
  }
 }
 assert(m.eof()); assert(int(B.size())-base==64);

 ifstream rs("restored68.txt");int restore=0;
 for(U u;rs>>u;){assert(popcount(u)==12);B.push_back({u,15,Code()});++restore;}
 assert(rs.eof()&&restore==352);
 std::set<unsigned> recovered;
 std::ifstream rr("recover.txt");
 for(unsigned j;rr>>j;){assert(j<B.size()&&B[j].t==3&&(B[j].h.s&2u));
                        assert(recovered.insert(j).second);}
 assert(rr.eof()&&recovered.size()==51);
 vector<Code> next(B.size());
 for(size_t i=0;i<B.size();i++)
   if(B[i].t&&!recovered.count(i))next[i].s=B[i].h.s&2u;
 auto value=[](const Code&h,V g)->int{
   for(unsigned i=0;i<h.n;i++)if(h.row[i]==g)return (h.s>>i)&1;
   assert(false);return 0;
 };
 ifstream phase("phase.txt");string p;phase>>p;
 assert(p.size()==guards.size()&&p.size()==31512);
 unsigned e=0;
 for(auto [ij,ab]:guards){
   auto [i,j]=ij; auto [a,b]=ab;
   int y=value(B[i].h,a),z=value(B[j].h,b);
   assert(y!=z&&(p[e]=='0'||p[e]=='1'));
   int d=p[e++]-'0';
   assert(next[i].ins(a,y^d)>=0&&next[j].ins(b,z^d)>=0);
 }
 string tail;assert(phase&&!(phase>>tail));
 for(size_t i=0;i<B.size();i++)B[i].h=next[i];
 auto extra=add_heads(B,C,cub);
 for(size_t i=0;i<B.size();i++)for(auto [v,b]:HG[i])assert(next[i].ins(v,b)>=0);
 for(size_t i=0;i<B.size();i++)B[i].h=next[i];
 // All zero-tail supports have the new even signs required by the spike donors.
 vector<U> dense_rows=cb;
 for(auto c:C)dense_rows.push_back(eval(c));
 set<pair<U,unsigned>> seen;
 for(auto s:B){
   assert(seen.insert({s.b,s.t}).second);
   assert(popcount(s.b)+popcount(s.t)==16);
   if(s.t)assert(popcount(s.t)==2||s.t==15);
   else{
     assert((s.h.s&2u)==0);
     for(U d:dense_rows)assert(!(popcount(d&s.b)&1));
   }
 }
 cout<<"extra signed overlap layers "<<B.size()-base<<endl;

 { std::ofstream zz("tag3.txt");
 for(size_t i=0;i<B.size();i++)if(B[i].t==3){
 zz<<i<<' '<<B[i].b<<' '<<B[i].h.n<<' '<<B[i].h.s;
 for(unsigned j=0;j<B[i].h.n;j++)zz<<' '<<B[i].h.row[j];
 zz<<'\n';}}
 for(auto& z:B)z.h.complete();
 size_t checked=0;for(int i=0;i<(int)B.size();i++)for(int j=0;j<i;j++){
  auto&a=B[i];auto&b=B[j];int k=popcount(a.b&b.b)+popcount(a.t&b.t);
  if(k<=8)continue;
  assert(k<=10&&guards.count({i,j}));++checked;
 }
 assert(checked==guards.size());cout<<"Checked signed overlap edges "<<checked<<"\n";
 // Sparse pair gates, balanced completion, and head-to-head support gates
 // have now all been replayed. Print reconstructed sparse population.
 assert(B.size()==68275&&guards.size()==31512);
 for(int n=65;n<=68;n++){
   int k=0;for(auto b:B)k+=(b.t<(1u<<(n-64)));
   cout<<"ATLAS "<<n<<" "<<k<<" "<<extra[n-64]<<"\n";
 }
}
