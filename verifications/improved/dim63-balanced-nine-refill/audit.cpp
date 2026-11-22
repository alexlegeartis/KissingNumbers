/* Exact scalar cross-check of all supports and signed phase comparisons. */
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <map>
using U=uint64_t;using V=__uint128_t;
int pc(V x){return __builtin_popcountll(U(x))+__builtin_popcountll(U(x>>64));}
U read(std::istream& s){U a=0;s.read(reinterpret_cast<char*>(&a),8);assert(s);return a;}
V readmask(std::istream&s){U a=read(s),b=read(s);assert(b==0);return (V(b)<<64)|a;}
int main(int argc,char**argv){
 assert(argc==2);
 std::ifstream f(argv[1], std::ios::binary);
 U nw=read(f),nt=read(f),ns=read(f); V H=readmask(f); assert(pc(H)==8);
 assert(nw==32801&&nt==10416);
 std::vector<V> W(nw),T(nt),S(ns);
 std::vector<std::vector<V>> P(ns);
 for(auto&x:W){x=readmask(f);assert(pc(x)==16);}
 for(auto&x:T){x=readmask(f);assert(pc(x)==4);}
 std::set<V> seen;std::map<V,std::set<V>> removed;
 U added=0,full_comparisons=0,relevant=0;
 for(U i=0;i<ns;i++){
  V s=S[i]=readmask(f);assert(pc(s)==9&&seen.insert(s).second);
  U n=read(f);assert(n>0&&n<=64);added+=n;
  for(auto&v:W)assert(pc(s&v)<=6);
  std::vector<V> blocks;
  for(auto v:T)if((s&v)==v)blocks.push_back(v);
  auto&u=P[i]; for(U j=0;j<n;j++){
   auto x=readmask(f);assert((s&x)==x && 2*pc(x&H)==pc(s&H));
   for(V y:u)assert(pc(x^y)>=3);
   u.push_back(x);
   for(V v:blocks)if(2*pc(v&x&H)==pc(v&H))removed[v].insert(x&v);
  }
  for(U j=0;j<i;j++){
   V q=s&S[j];int k=pc(q);if(k<=4)continue;
   relevant++;
   for(V a:u)for(V b:P[j]){assert(k-2*pc((a^b)&q)<=4);full_comparisons++;}
  }
 }
 assert(f.peek()==std::char_traits<char>::eof());
 U loss=0;for(auto&x:removed){assert(x.second.size()<=16);loss+=x.second.size();}
 std::cout<<"{\"dimension\":63,\"added\":"<<added-2084<<",\"removed\":"<<loss<<",\"total\":"<<184417784+added-loss-2084<<",\"relevant_support_pairs\":"<<relevant<<",\"all_signed_pairs_on_relevant_supports\":"<<full_comparisons<<"}\n";
}

