#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
using namespace std;
int main(int argc,char **argv){if(argc<3)return 2;unsigned char map[24][8];ifstream f(argv[1],ios::binary);f.read((char*)map,sizeof(map));if(!f)return 2;int runs=atoi(argv[2]),gray[8],g[]={7,1,4,4,2,7,3,3,4,7,2,1},G[12][24]={};mt19937_64 rng(9322191);for(int i=0;i<12;i++){for(int j=0;j<12;j++)G[i][i+j]=g[j];G[i][23]=3;}for(int z=0;z<8;z++){gray[z]=0;for(int x=0;x<2;x++)for(int y=0;y<2;y++)gray[z]|=(((z>>2&1)^((z&1)*x)^((z>>1&1)*y))<<(2*x+y));}
 int best=96;for(int run=0;run<runs;run++){int ca[12]={},cb[12]={},a[24]={},b[24]={};for(int i=0;i<10;i++){int va[6],vb[6],na=0,nb=0;for(int z=0;z<8;z++){if(__builtin_popcount((unsigned)map[i][z])==2)va[na++]=z;if(__builtin_popcount((unsigned)gray[z])==2)vb[nb++]=z;}if(na!=6||nb!=6)return 3;ca[i]=(7*(va[rng()%6]-a[i]))&7;cb[i]=(7*(vb[rng()%6]-b[i]))&7;for(int j=0;j<24;j++){a[j]=(a[j]+ca[i]*G[i][j])&7;b[j]=(b[j]+cb[i]*G[i][j])&7;}}
 unsigned char xa[64][24],yb[64][24];for(int k=0;k<64;k++)for(int j=0;j<24;j++){xa[k][j]=map[j][(a[j]+(k&7)*G[10][j]+(k>>3)*G[11][j])&7];yb[k][j]=gray[(b[j]+(k&7)*G[10][j]+(k>>3)*G[11][j])&7];}
 for(int u=0;u<64;u++)for(int v=0;v<64;v++){int d=0;for(int j=0;j<24&&d<24;j++)d+=__builtin_popcount((unsigned)(xa[u][j]^yb[v][j]));best=min(best,d);if(d<24){ca[10]=u&7;ca[11]=u>>3;cb[10]=v&7;cb[11]=v>>3;printf("WITNESS A ");for(int x:ca)printf("%d",x);printf(" B ");for(int x:cb)printf("%d",x);printf(" distance %d prefixes %d\n",d,run+1);return 0;}}}
 printf("UNRESOLVED best %d prefixes %d\n",best,runs);
}
