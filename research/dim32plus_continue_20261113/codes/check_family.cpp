// Independent scalar verification of every support interaction and parent word.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;using U=uint64_t;
void require(bool x){if(!x)throw runtime_error("certificate rejected");}
int wt(U x){return __builtin_popcountll(x);}
U expand(U support,unsigned word){U out=0;unsigned j=0;for(int i=0;i<64;i++)if(support>>i&1){if(word>>j&1)out|=U(1)<<i;j++;}return out;}
struct Code{U support,phase;vector<U>span;};
Code read_code(istream&in,bool odd){
 U s;unsigned shift;array<unsigned,5>checks;in>>s;for(auto&h:checks)in>>h;in>>shift;
 require(bool(in)&&wt(s)==16&&shift<65536&&checks[0]==65535);
 for(unsigned h:checks)require(h<65536);
 set<unsigned>columns;
 for(int j=0;j<16;j++){unsigned v=0;for(int k=0;k<5;k++)v|=((checks[k]>>j)&1)<<k;require(v&1);columns.insert(v);}
 require(columns.size()==16);vector<U>span;set<unsigned>seen;
 for(int mask=0;mask<32;mask++){unsigned v=0;for(int k=0;k<5;k++)if(mask>>k&1)v^=checks[k];require(wt(v)==0||wt(v)==8||wt(v)==16);seen.insert(v);span.push_back(expand(s,v));}
 require(seen.size()==32&&(!odd||(wt(shift)&1)));sort(span.begin(),span.end());
 return{s,expand(s,shift),span};
}
void repair(const Code&a,const Code&b){
 U overlap=a.support&b.support;require(wt(overlap)==10);int valid=0;
 for(U h:a.span)if(wt(h)==8&&(h&overlap)==h&&binary_search(b.span.begin(),b.span.end(),h)){
  require((wt(h&a.phase)&1)!=(wt(h&b.phase)&1));valid++;
 }
 require(valid==1);
}
int main(int argc,char**argv){
 try{
  require(argc==2);ifstream in(argv[1]);size_t no,nn,nm,ng;bool odd;U removed;
  in>>no>>nn>>nm>>ng>>odd>>removed;require(bool(in)&&no==30828&&ng==28&&nn>0&&nn<112897&&nm<=no);
  vector<U>gen(ng),old(no);for(U&x:gen)in>>x;for(U&x:old)in>>x;require(bool(in));
  vector<U>generated;array<U,65>hist{};hist[0]=1;U word=0;
  for(uint32_t i=1;i<(uint32_t(1)<<28);i++){
   word^=gen[__builtin_ctz(i)];int w=wt(word);require(w>=16);hist[w]++;if(w==16)generated.push_back(word);
  }
  sort(generated.begin(),generated.end());auto sorted_old=old;sort(sorted_old.begin(),sorted_old.end());require(generated==sorted_old);
  require(set<U>(old.begin(),old.end()).size()==no);
  map<U,Code>modified;
  for(size_t i=0;i<nm;i++){Code c=read_code(in,odd);require(!(c.support&removed)&&binary_search(sorted_old.begin(),sorted_old.end(),c.support));require(modified.emplace(c.support,c).second);}
  vector<Code>added;set<U>unique;
  for(size_t i=0;i<nn;i++){Code c=read_code(in,odd);require(!(c.support&removed)&&!binary_search(sorted_old.begin(),sorted_old.end(),c.support)&&unique.insert(c.support).second);for(U row:gen)require(!(wt(row&c.support)&1));added.push_back(c);}
  U old_pairs=0,new_pairs=0,old_repairs=0,new_repairs=0;set<U>touched;
  for(size_t i=0;i<added.size();i++){
   const auto&a=added[i];
   for(U b:old)if(!(b&removed)){old_pairs++;int overlap=wt(a.support&b);require(overlap<=10);if(overlap>8){auto it=modified.find(b);require(it!=modified.end());repair(a,it->second);touched.insert(b);old_repairs++;}}
   for(size_t j=0;j<i;j++){new_pairs++;int overlap=wt(a.support&added[j].support);require(overlap<=10);if(overlap>8){repair(a,added[j]);new_repairs++;}}
  }
  require(touched.size()==modified.size());
  cout<<"{\"status\":\"PASS\",\"parent_words\":268435456,\"parent_minimum\":16,\"parent_weight16_words\":"<<generated.size()<<",\"old_new_pairs\":"<<old_pairs<<",\"new_new_pairs\":"<<new_pairs<<",\"old_new_repairs\":"<<old_repairs<<",\"new_new_repairs\":"<<new_repairs<<"}"<<endl;
 }catch(const exception&e){cerr<<e.what()<<endl;return 1;}
}
