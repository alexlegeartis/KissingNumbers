#include<bits/stdc++.h>
using namespace std;
using V=uint64_t;
vector<V> make(int vars,int level) {
 vector<V> a{0};
 for(int m=0;m<(1<<vars);m++)if(__builtin_popcount(m)<=level){
  V u=0;for(int x=0;x<(1<<vars);x++)if((x&m)==m)u|=V(1)<<x;
  int n=a.size();for(int i=0;i<n;i++)a.push_back(a[i]^u);
 }
 return a;
}
int main(int argc,char**argv){
 vector<int> cubic;
 for(int x=0;x<64;x++)if(__builtin_popcount(x)==3)cubic.push_back(x);
 auto l=make(5,1),q=make(5,2);
 assert(l.size()==64&&q.size()==65536);
 ifstream f(argc>1?argv[1]:"cubic80.txt");vector<V>C;
 for(unsigned u;f>>u;){
  V s=0;for(int i=0;i<64;i++){int b=0;for(int j=0;j<20;j++)if((u>>j&1)&&(i&cubic[j])==cubic[j])b^=1; s|=V(b)<<i;}
  C.push_back(s);
 }
 assert(C.size()>=80&&set<V>(C.begin(),C.end()).size()==C.size());
 V t=argc>=3?stoull(argv[2]):0;
 vector<V> words;
 for(V c:C)for(V a:l)for(V u:q){
  V s=((u^a)<<32)|u;
  V z=s^c^t;
  if(__builtin_popcountll(z)==16)words.push_back(z);
 }
 sort(words.begin(),words.end());assert(adjacent_find(words.begin(),words.end())==words.end());
 ofstream fw("supports64.txt"); for(V s:words)fw<<s<<"\n";
 const int M=(words.size()+63)/64;
 vector<vector<V>> Z(64,vector<V>(M)),O=Z;
 for(int i=0;i<64;i++)for(int j=0;j<(int)words.size();j++)
  ((words[j]>>i&1)?O:Z)[i][j>>6]|=V(1)<<(j&63);
 int best2=-1,best4=-1;array<int,2>a2{};array<int,4>a4{};int sz2a=0,sz2b=0;
 vector<V>b2(M),b3(M);
 for(int a=0;a<64;a++)for(int b=a+1;b<64;b++){
  int n=0,m=0;
  for(int i=0;i<M;i++){b2[i]=Z[a][i]&Z[b][i];n+=__builtin_popcountll(b2[i]);m+=__builtin_popcountll(O[a][i]&O[b][i]);}
  if(2*n+m>best2)best2=2*n+m,a2={a,b},sz2a=n,sz2b=m;
  for(int c=b+1;c<64;c++){
   for(int i=0;i<M;i++)b3[i]=b2[i]&Z[c][i];
   for(int d=c+1;d<64;d++){
    int nn=0;for(int i=0;i<M;i++)nn+=__builtin_popcountll(b3[i]&Z[d][i]);
    if(nn>best4)best4=nn,a4={a,b,c,d};
   }
  }
 }
 cout<<"64 supports "<<words.size()<<" 63 best "<<a2[0]<<" "<<a2[1]<<" none "<<sz2a<<" both "<<sz2b<<" 62 best ";
 for(int x:a4)cout<<x<<" ";cout<<"none "<<best4<<endl;
}
