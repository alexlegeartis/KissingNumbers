// Independent direct integer audit. No transform, floating point, or producer imports.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using I = std::int64_t;
using V = std::array<std::int32_t,39>;
static void need(bool c,const char *s){if(!c)throw std::runtime_error(s);}
static I dot(const V &a,const V &b){I s=0;for(int k=0;k<39;++k)s+=I(a[k])*b[k];return s;}
int main(int argc,char **argv){try{
 need(argc==2,"Usage: audit_cpp input.txt");std::ifstream in(argv[1]);need(bool(in),"open input");
 int n,q,m,nr,nh;I L,U;in>>n>>q>>m>>nr>>nh>>L>>U;
 need(n==10269&&q==32&&m==3424&&nr==2908&&nh==64,"unexpected sizes");
 need(L==99999996&&U==100000004,"norm envelope");
 const I limit=L*L/2, crosslimit=8*L*L; // sparse norms are 32, so d^2 <= 8 L^2.
 std::vector<std::uint64_t> G(q);for(auto &g:G)in>>g;
 std::vector<V>A(n);for(auto &a:A)for(auto &x:a)in>>x;
 std::vector<std::array<int,8>>S(m);for(auto &s:S)for(auto &k:s)in>>k;
 std::vector<V>R(nr),H(nh);for(auto &r:R)for(auto &x:r)in>>x;for(auto &h:H)for(auto &x:h)in>>x;
 need(bool(in),"truncated input");std::string extra;need(!(in>>extra),"trailing input");
 for(const auto &a:A)for(auto x:a)need(x!=0 && std::abs(I(x))<=U,"coordinate safety bound");
 I normmin=std::numeric_limits<I>::max(),normmax=0;
 for(auto &a:A){I t=dot(a,a);need(L*L<=t&&t<=U*U,"norm enclosure");normmin=std::min(normmin,t);normmax=std::max(normmax,t);}
 std::vector<std::array<int,2>>rootcoords(nr);
 for(int r=0;r<nr;++r){int z=0;for(int k=0;k<39;++k)if(R[r][k]){need(z<2,"root shape");rootcoords[r][z++]=k;}need(z==2&&dot(R[r],R[r])==32,"root norm");}
 std::vector<std::array<int,5>>heavycoords(nh);
 for(int h=0;h<nh;++h){int z=0;for(int k=0;k<39;++k)if(H[h][k]){need(z<5,"heavy shape");heavycoords[h][z++]=k;}need(z==5&&dot(H[h],H[h])==32,"heavy norm");}
 I maxrr=std::numeric_limits<I>::min(),maxrh=maxrr,maxhh=maxrr,maxhm=maxrr;
 for(int r=0;r<nr;++r)for(int t=r+1;t<nr;++t){I d=0;for(int k:rootcoords[r])d+=I(R[r][k])*R[t][k];maxrr=std::max(maxrr,d);need(d<=16,"root-root");}
 for(int h=0;h<nh;++h){for(int r=0;r<nr;++r){I d=0;for(int k:rootcoords[r])d+=I(H[h][k])*R[r][k];maxrh=std::max(maxrh,d);need(d<=16,"heavy-root");}
 for(int t=h+1;t<nh;++t){I d=dot(H[h],H[t]);maxhh=std::max(maxhh,d);need(d<=16,"heavy-heavy");}
 for(auto &s:S){for(unsigned mask=0;mask<256;++mask){unsigned p=mask;p^=p>>4;p^=p>>2;p^=p>>1;if(p&1)continue;I d=0;for(int k=0;k<8;++k)d+=I(H[h][s[k]])*((mask>>k)&1?-2:2);maxhm=std::max(maxhm,d);need(d<=16,"heavy-middle");}}}
 std::cerr<<"Sparse direct checks passed\n";
 // A group mask restricted to a middle support changes only sign parity.
 // Enumerate every allowed 8-sign vector directly by Gray-code additions.
 std::vector<bool>both(m,false);int supports_both=0;
 for(int s=0;s<m;++s){for(auto g:G){unsigned p=0;for(int k:S[s])p^=unsigned((g>>k)&1);both[s]=both[s]||p;}supports_both+=both[s];}
 I maxmiddle=0;int middlei=-1,middles=-1;std::uint64_t middleterms=0;
 for(int i=0;i<n;++i)for(int s=0;s<m;++s){I cur=0;for(int k:S[s])cur+=A[i][k];I best=cur;unsigned prev=0;
 for(unsigned t=1;t<256;++t){unsigned gray=t^(t>>1),change=gray^prev;int b=0;while(((change>>b)&1)==0)++b;cur+=((gray>>b)&1?-2:2)*I(A[i][S[s][b]]);prev=gray;
 // Gray vector's parity equals the low bit of its sequence index.
 if(both[s]||!(t&1))best=std::max(best,cur);}
 middleterms+=both[s]?256:128;
 need(best>=0&&best<=200000000,"middle bound overflow guard");need(best*best<=2*L*L,"dense-middle");
 if(best>maxmiddle){maxmiddle=best;middlei=i;middles=s;}}
 std::cerr<<"Dense-middle direct sign enumerations passed\n";
 // Test every dense translate against every actual heavy point directly.
 // Roots include all four signs on each retained pair, so their exact maximum
 // is 4 (|a_p|+|a_q|), unchanged by every sign-group element.
 I maxroot=0,maxheavy=0;std::uint64_t rootchecks=0,heavychecks=0;
 for(int i=0;i<n;++i){for(int r=0;r<nr;r+=4){auto xy=rootcoords[r];I d=4*(std::abs(I(A[i][xy[0]]))+std::abs(I(A[i][xy[1]])));need(d<=1000000000,"root overflow guard");need(d*d<=crosslimit,"dense-root");maxroot=std::max(maxroot,d);++rootchecks;}
 for(auto g:G){V a=A[i];for(int k=0;k<39;++k)if((g>>k)&1)a[k]=-a[k];for(int h=0;h<nh;++h){I d=0;for(int k:heavycoords[h])d+=I(a[k])*H[h][k];need(std::abs(d)<=1000000000,"heavy overflow guard");if(d>0)need(d*d<=crosslimit,"dense-heavy");maxheavy=std::max(maxheavy,d);++heavychecks;}}}
 std::cerr<<"Dense-root and dense-heavy direct checks passed\n";
 // Store all transformed right rows once. Dot accumulation is signed int64.
 // Norm enclosure implies |coordinate| <= U = 100,000,004.
 // This is checked before dot calls; 39*U^2 = 390000031200000624 < 2^63. No squared dense dot products are evaluated.
 std::vector<V>B(std::size_t(q)*n);
 for(int g=0;g<q;++g)for(int j=0;j<n;++j){B[std::size_t(g)*n+j]=A[j];for(int k=0;k<39;++k)if((G[g]>>k)&1)B[std::size_t(g)*n+j][k]=-A[j][k];}
 I maxdense=std::numeric_limits<I>::min();std::uint64_t products=0;int wi=-1,wj=-1,wg=-1;
 for(int i=0;i<n;++i){for(int g=0;g<q;++g){int first=i+(G[g]==0);for(int j=first;j<n;++j){I d=dot(A[i],B[std::size_t(g)*n+j]);need(d<=limit,"dense-dense");if(d>maxdense){maxdense=d;wi=i;wj=j;wg=g;}++products;}}
 if(i%1024==0)std::cerr<<"Dense rows "<<i<<" / "<<n<<"\n";}
 need(products==std::uint64_t(q)*n*(n-1)/2+std::uint64_t(q-1)*n,"pair count");
 std::cout<<"{\n\"passed\": true,\n\"arithmetic\": \"direct signed int64, no transforms or floating point\",\n"
 <<"\"norm_squared_min\": "<<normmin<<",\n\"norm_squared_max\": "<<normmax<<",\n"
 <<"\"dense_relative_products\": "<<products<<",\n\"max_dense_dot\": "<<maxdense<<",\n\"dense_dot_limit\": "<<limit<<",\n"
 <<"\"dense_witness\": ["<<wi<<","<<wj<<","<<G[wg]<<"],\n"
 <<"\"max_root_root_dot\": "<<maxrr<<",\n\"max_heavy_root_dot\": "<<maxrh<<",\n\"max_heavy_heavy_dot\": "<<maxhh<<",\n\"max_heavy_middle_dot\": "<<maxhm<<",\n"
 <<"\"middle_supports_allowing_both_parities\": "<<supports_both<<",\n\"middle_sign_sums_evaluated\": "<<middleterms<<",\n"
 <<"\"max_dense_middle_half_dot\": "<<maxmiddle<<",\n\"middle_witness\": ["<<middlei<<","<<middles<<"],\n\"middle_square_slack\": "<<2*L*L-maxmiddle*maxmiddle<<",\n"
 <<"\"max_dense_root_dot\": "<<maxroot<<",\n\"dense_root_coordinate_pair_checks\": "<<rootchecks<<",\n"
 <<"\"max_dense_heavy_dot\": "<<maxheavy<<",\n\"dense_heavy_direct_checks\": "<<heavychecks<<"\n}\n";
 return 0;
}catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}
