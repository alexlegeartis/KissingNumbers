#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
using namespace std;
int main(int argc,char**argv){if(argc<3)return 2;int runs=atoi(argv[2]),target[24];for(int j=0;j<24;j++){char c=argv[1][j];target[j]=(c>='0'&&c<='9')?c-'0':10+c-'a';if(j<4&&__builtin_popcount((unsigned)target[j])!=2)return 3;if(j>=4&&__builtin_popcount((unsigned)target[j])%2!=1)return 3;}mt19937_64 R(argc>3?strtoull(argv[3],nullptr,10):920020);int g[]={7,1,4,4,2,7,3,3,4,7,2,1},G[12][24]={},gray[8];for(int i=0;i<12;i++){for(int j=0;j<12;j++)G[i][i+j]=g[j];G[i][23]=3;}for(int z=0;z<8;z++){gray[z]=0;for(int x=0;x<2;x++)for(int y=0;y<2;y++)gray[z]|=(((z>>2&1)^((z&1)*x)^(((z>>1&1)*y)))<<(2*x+y));}
 int best=96;for(int run=0;run<runs;run++){int c[12]={},a[24]={};for(int i=0;i<12;i++){if(i<4){int valid[6],n=0;for(int z=0;z<8;z++)if(__builtin_popcount((unsigned)gray[z])==2)valid[n++]=z;c[i]=(7*(valid[R()%n]-a[i]))&7;}else c[i]=R()&7;for(int j=0;j<24;j++)a[j]=(a[j]+c[i]*G[i][j])&7;}auto cost=[&](){int d=0;for(int j=0;j<24;j++)d+=__builtin_popcount((unsigned)(gray[a[j]]^target[j]));return d;};int cur=cost();for(int it=0;it<200000;it++){int i=4+R()%8,v=R()&7,dv=v-c[i],nd=0;if(!dv)continue;for(int j=0;j<24;j++)nd+=__builtin_popcount((unsigned)(gray[(a[j]+dv*G[i][j])&7]^target[j]));double T=.3+2.0*(1.0-double(it%20000)/20000);if(nd<=cur||double(R()>>11)*(1.0/9007199254740992.0)<exp((cur-nd)/T)){for(int j=0;j<24;j++)a[j]=(a[j]+dv*G[i][j])&7;c[i]=v;cur=nd;}best=min(best,cur);if(cur<24){printf("WITNESS B ");for(int x:c)printf("%d",x);printf(" distance %d\n",cur);return 0;}}}printf("UNRESOLVED best %d\n",best);
}
