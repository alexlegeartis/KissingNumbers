#include <bits/stdc++.h>
using namespace std;using Z=__int128_t;
Z R,c[24],D;
int L[24][24],z[24];
int bits[24];ofstream out;uint64_t countv[17]{};
Z read(){string s;cin>>s;Z t=0;for(char a:s)t=10*t+a-'0';return t;}
Z flo(Z a,Z b){if(a<0)return -((b-1-a)/b);return a/b;}
void go(int k,Z left){
 if(k<0){Z n=R-left;assert(n%D==0);int d=n/D;assert(0<=d&&d<=16);countv[d]++;if(d>0){
    int32_t buf[25];buf[0]=d;for(int i=0;i<24;i++)buf[i+1]=z[i];out.write((char*)buf,sizeof(buf));
 }return;}
 auto q=left/c[k];uint64_t r=sqrtl((long double)q);
 while(Z(r+1)*(r+1)<=q)r++;while(Z(r)*r>q)r--;
 int sh=0;for(int j=k+1;j<24;j++)sh+=L[k][j]*z[j];
 int step=L[k][k]; int lo=-flo(Z(r)+sh,step),hi=flo(Z(r)-sh,step);
 if(bits[k]>=0&&(lo&1)!=bits[k])lo++;
 for(int v=lo;v<=hi;v+=(bits[k]<0?1:2)){z[k]=v; Z a=sh+step*v; go(k-1,left-c[k]*a*a);}
}
int main(int n,char**a){
 assert(n==4);out.open(a[2],ios::binary);
 int k=stoi(a[1]);for(int i=0;i<24;i++)bits[i]=k<0?-1:((k>>i)&1);
 D=read();R=stoi(a[3])*D;
 for(int i=0;i<24;i++){c[i]=read();for(int j=i;j<24;j++)cin>>L[i][j];}
 go(23,R);
 for(int i=0;i<17;i++)if(countv[i])cerr<<i<<':'<<countv[i]<<' ';cerr<<'\n';
}
