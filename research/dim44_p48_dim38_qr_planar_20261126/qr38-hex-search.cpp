
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <string>
using namespace std;
constexpr unsigned Q=1u<<24;
using U=uint64_t;
struct P{unsigned q;int phase;};
struct C{unsigned q;int h,i=-1,j=-1;};
vector<uint8_t> D(Q,255);
vector<array<int,3>> H;
const array<array<int,3>,6> B={{{-16,8,8},{8,-16,8},{-8,-8,16},{8,8,-16},{-8,16,-8},{16,-8,-8}}};
int dot(array<int,3> a,array<int,3> b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
int need(array<int,3>a,array<int,3>b){long long s=5184+dot(a,a),t=5184+dot(b,b);for(int d=0;d<=12;d++){long long g=5184-288*d+dot(a,b);if(g<=0||4*g*g<=s*t)return d;}assert(false);return 1000;}
bool relaxed(array<int,6> m){
 // Scaled-head coordinates z=3*head, ||z||²<=96. A bright ±12 coordinate
 // functional gives exact outer intervals when the radius cap is substituted.
 int L[3]={50-6*m[0],50-6*m[1],50-6*m[3]};
 int U[3]={6*m[5]-50,6*m[4]-50,6*m[2]-50};
 for(int i=0;i<3;i++){L[i]=max(L[i],-8);U[i]=min(U[i],8);if(L[i]>U[i])return false;}
 if(L[0]+L[1]+L[2]>0||U[0]+U[1]+U[2]<0)return false;
 for(int c=0;c<27;c++){
  int act[3]={c%3,c/3%3,c/9},F=0,S=0,n=0;bool ok=true;
  for(int i=0;i<3;i++){if(act[i]==1){n++;continue;}int y=act[i]==0?L[i]:U[i];S+=y;F+=y*y;}
  if(n){for(int i=0;i<3;i++)if(act[i]==1&&(-S<L[i]*n||-S>U[i]*n))ok=false;
   if(ok&&F*n+S*S<=96*n)return true;
  }else if(S==0&&F<=96)return true;
 }
 return false;
}
int main(){
 ios::sync_with_stdio(false);cin.tie(nullptr);
 int dim,a,p;cin>>dim>>a>>p;assert(dim==36&&a==62&&p<=64&&p>0);
 vector<unsigned> col(dim);for(auto&x:col)cin>>x;
 vector<P> base(a);for(auto&x:base){cin>>x.q>>x.phase;assert(x.phase>=0&&x.phase<6);}
 H.resize(p);for(auto&x:H){cin>>x[0]>>x[1]>>x[2];assert(x[0]+x[1]+x[2]==0&&dot(x,x)<=384);}
 assert(cin); vector<unsigned> queue(Q);unsigned tail=1;queue[0]=0;D[0]=0;
 for(unsigned i=0;i<tail;i++){unsigned x=queue[i];uint8_t n=D[x]+1;for(unsigned g:col){unsigned t=x^g;if(D[t]==255){D[t]=n;queue[tail++]=t;}}}
 assert(tail==Q);int cov=*max_element(D.begin(),D.end());vector<unsigned>().swap(queue);
 vector<long long> hist(cov+1);for(int q:D)hist[q]++;
 cerr<<"BFS complete "<<tail<<" cover "<<cov<<" distances";for(auto w:hist)cerr<<" "<<w;cerr<<'\n';
 for(int i=0;i<a;i++)for(int j=0;j<i;j++)assert(D[base[i].q^base[j].q]>=need(B[base[i].phase],B[base[j].phase]));
 // Interleave the six hues so a failure in any tier is seen early.
 vector<P> order;for(int d=0;d<62;d++)for(int t=0;t<6;t++){int k=0;for(P x:base)if(x.phase==t&&k++==d){order.push_back(x);break;}}
 assert(order.size()==62);
 U mask[6][40]={};vector<vector<int>> nd(p,vector<int>(6));
 for(int j=0;j<6;j++)for(int k=0;k<p;k++){nd[k][j]=need(H[k],B[j]);for(int d=0;d<nd[k][j];d++)mask[j][d]|=U(1)<<k;}
 vector<vector<int>> DD(p,vector<int>(p));for(int i=0;i<p;i++)for(int j=0;j<p;j++)DD[i][j]=need(H[i],H[j]);
 const U all=p==64?~U(0):(U(1)<<p)-1;
 const size_t LIM=200000;vector<C> cands;unsigned co[3]={};unsigned relaxedN=0;vector<pair<unsigned,array<int,6>>> rs;
 for(unsigned q=0;q<Q;q++){
  U u0=all,u1=0,u2=0;bool okrel=true,skip=false;array<int,6> m={10,10,10,10,10,10};
  for(P x:order){
   int d=D[q^x.q];m[x.phase]=min(m[x.phase],d);if(d<7)okrel=false;
   U k=mask[x.phase][d];u2=(u2&~k)|(u1&k);u1=(u1&~k)|(u0&k);u0&=~k;
   if(!(u0|u1|u2)&&!okrel){skip=true;break;}
  }
  if(skip)continue;
  if(okrel&&relaxed(m)){relaxedN++;if(rs.size()<500)rs.push_back({q,m});}
  U u=u0|u1|u2;
  while(u){int k=__builtin_ctzll(u);u&=u-1;int len=(u0>>k&1)?0:(u1>>k&1)?1:2;
   C c{q,k};int n=0;bool identical=false;
   for(int i=0;i<a;i++){
    P x=base[i];if(q==x.q&&H[k]==B[x.phase]){identical=true;break;}
    if(D[q^x.q]<nd[k][x.phase]){if(n++==0)c.i=i;else c.j=i;}
   }
   if(identical)continue;
   assert(n==len);co[len]++;
   if(cands.size()<LIM)cands.push_back(c);
  }
 }
 cerr<<"scan "<<Q<<" plane outer "<<relaxedN<<" menu candidates "<<co[0]<<" "<<co[1]<<" "<<co[2]<<" saved "<<cands.size()<<'\n';
 auto compat=[&](const C&a,const C&b){return D[a.q^b.q]>=DD[a.h][b.h];};
 vector<C> best;vector<int> del;
 auto test=[&](vector<int> out,vector<C> in){if((int)in.size()-(int)out.size()>(int)best.size()-(int)del.size()){best=in;del=out;}};
 vector<C> zero;vector<vector<C>> one(a),two(a*a);
 for(C x:cands){if(x.i<0)zero.push_back(x);else if(x.j<0)one[x.i].push_back(x);else two[x.i*a+x.j].push_back(x);}
 if(!zero.empty()){
  for(auto x:zero){vector<C> v{x};for(auto y:zero)if(all_of(v.begin(),v.end(),[&](const C&z){return compat(y,z);}))v.push_back(y);test({},v);if(v.size()>=10)break;}
 }
 for(int i=0;i<a;i++)for(size_t j=0;j<one[i].size();j++)for(size_t k=0;k<j;k++)
  if(compat(one[i][j],one[i][k])){vector<C>d={one[i][j],one[i][k]};for(C x:one[i])if(all_of(d.begin(),d.end(),[&](C y){return compat(x,y);}))d.push_back(x);test({i},d);}
 // Restrict the seeded 2→3 exchange to the first 250 distinct witnesses per outgoing pair; terminal output scopes this as a search when truncated.
 long long groupN=0,trunc=0,bin0=0;
 for(int i=0;i<a;i++)for(int j=i+1;j<a;j++){
  vector<C> q=two[i*a+j];if(q.empty()&&one[i].empty()&&one[j].empty()&&zero.empty())continue;
  q.insert(q.end(),one[i].begin(),one[i].end());q.insert(q.end(),one[j].begin(),one[j].end());q.insert(q.end(),zero.begin(),zero.end());
  if(q.size()<3)continue;groupN++;if(q.size()>250){q.resize(250);trunc++;}
  for(size_t k=0;k<q.size();k++)for(size_t t=0;t<k;t++)if(compat(q[k],q[t]))for(size_t s=0;s<t;s++)if(compat(q[k],q[s])&&compat(q[t],q[s])){
   vector<C>d={q[k],q[t],q[s]};for(C x:q)if(all_of(d.begin(),d.end(),[&](C y){return compat(x,y);}))d.push_back(x);test({i,j},d);
  }
 }
 for(C x:best){for(int i=0;i<a;i++)if(find(del.begin(),del.end(),i)==del.end())assert(D[x.q^base[i].q]>=nd[x.h][base[i].phase]);
  for(C y:best)if(x.q!=y.q||x.h!=y.h)assert(compat(x,y));}
 cout<<"{\"syndromes\":"<<Q<<",\"cover\":"<<cov<<",\"relaxed\":"<<relaxedN<<",\"menu\":["<<co[0]<<','<<co[1]<<','<<co[2]<<"],\"saved\":"<<cands.size()<<",\"outgoing_pairs_tested\":"<<groupN<<",\"truncated_pairs\":"<<trunc<<",\"deleted\":[";
 for(size_t i=0;i<del.size();i++){if(i)cout<<',';cout<<del[i];}
 cout<<"],\"added\":[";for(size_t i=0;i<best.size();i++){if(i)cout<<',';C x=best[i];const auto z=H[x.h];cout<<'['<<x.q<<','<<z[0]<<','<<z[1]<<','<<z[2]<<']';}
 cout<<"],\"relaxed_examples\":[";for(size_t i=0;i<rs.size();i++){if(i)cout<<',';cout<<'['<<rs[i].first;for(int d:rs[i].second)cout<<','<<d;cout<<']';}cout<<"]}\n";
}

