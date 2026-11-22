#include <bits/stdc++.h>
using namespace std; int pc(int x){return __builtin_popcount(x);} struct V{int s,x,w;};
int main(){mt19937 g(123); static unsigned char DD[4096][4096]; int best=384;
for(int run=0;run<2000;run++){vector<int> bas,G{0}; for(int i=0;i<6;i++){int r; do{r=g()%4096;}while(find(G.begin(),G.end(),r)!=G.end());bas.push_back(r);for(int j=0;j<(1<<i);j++)G.push_back(G[j]^r);}
for(int s=0;s<4096;s++)if(pc(s)>=5){auto d=DD[s]; for(int x=s;;x=(x-1)&s){d[x]=100;if(!x)break;} for(int h:G)d[h&s]=0;for(int y=s;y;y&=y-1){int b=y&-y;for(int x=s;;x=(x-1)&s){if(x&b){auto a=d[x],c=d[x^b];d[x]=min((int)a,c+1);d[x^b]=min((int)c,a+1);}if(!x)break;}}}
vector<V> pool;for(int s=0;s<4096;s++)if(pc(s)==8||pc(s)==9){int c=0,md=12;for(int h:G){if((h&s)==0)c++;else md=min(md,pc(h&s));}if(pc(s)-2*md>4)continue;for(int x=s;;x=(x-1)&s){bool ok=true;for(int h:G)if(((h&s)^x)<x){ok=false;break;}if(ok)pool.push_back({s,x,64/c});if(!x)break;}}
for(int iter=0;iter<100;iter++){shuffle(pool.begin(),pool.end(),g);vector<V> chosen;int total=0;for(auto a:pool){bool ok=true;for(auto b:chosen){int s=a.s&b.s,p=pc(s);if(p>4&&p-2*DD[s][(a.x^b.x)&s]>4){ok=false;break;}}if(ok){chosen.push_back(a);total+=a.w;}}if(total>best){best=total;cerr<<run<<" "<<iter<<" "<<total<<endl;ofstream f("win.txt");for(int b:bas)f<<b<<" ";f<<endl;for(auto a:chosen)f<<a.s<<" "<<a.x<<" "<<a.w<<endl;}}
if(run%20==0)cerr<<"done "<<run<<" pool "<<pool.size()<<endl;
}}
