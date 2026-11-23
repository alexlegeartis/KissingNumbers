// Atlas geometry and head separation based on the published complete-atlas verifier.
static vector<pair<V,int>> HG[70000];
static U selected_sites[80]={};
static set<int> selected_heads;
struct Fibre {
 U p;
 uint32_t v;
 unsigned axes;
};
vector<U> add_heads(vector<Support>& B,const vector<unsigned>& C,const vector<U>& cub){
 using W=uint32_t;
 auto mon=[](int i){W w=0;for(int k=0;k<32;k++)if((k&i)==i)w|=1u<<k;return w;};
 vector<W> qb,cb;
 for(int i=0;i<32;i++)(popcount(unsigned(i))<=2?qb:cb).push_back(mon(i));
 assert(qb.size()==16&&cb.size()==16);
 auto span=[](const vector<W>&bs){vector<W>A={0};for(auto b:bs){auto n=A.size();for(size_t j=0;j<n;j++)A.push_back(A[j]^b);}return A;};
 const auto Q=span(qb),X=span(cb);
 assert(set<W>(Q.begin(),Q.end()).size()==65536&&set<W>(X.begin(),X.end()).size()==65536);
 for(W x:Q)if(x)assert(popcount(x)>=8);
 auto gf=[](unsigned a,unsigned b){
  unsigned c=0;for(;b;b>>=1){if(b&1)c^=a;a<<=1;if(a&32)a^=37;}
  return c;
 };
 vector<W> K={~W(0)};for(int i=0;i<5;i++)K.push_back(mon(1<<i));
 for(int i=0;i<5;i++){
  W w=0;
  for(unsigned x=0;x<32;x++){
   unsigned z=gf(gf(gf(x,x),x),1u<<i),t=0;
   for(int j=0;j<5;j++){t^=z;z=gf(z,z);}
   assert(t<=1);w|=t<<x;
  }
  K.push_back(w);
 }
 auto H=span(K);
 set<W> qset(Q.begin(),Q.end()), coset(H.begin(),H.end());
 assert(H.size()==2048&&coset.size()==2048);
 for(W h:H)assert(qset.count(h)&&(!h||popcount(h)>=12));
 vector<W> phase;
 for(W q:Q){
  if(phase.size()==24)break;
  bool fresh=true;for(W a:phase)if(coset.count(a^q)){fresh=false;break;}
  if(fresh)phase.push_back(q);
 }
 assert(phase.size()==24);
 for(size_t i=0;i<phase.size();i++)for(size_t j=0;j<i;j++)for(W h:H)
   assert(popcount(phase[i]^phase[j]^h)>=8);
 // Heads are sqrt(2) times the following integer coordinates.
 vector<array<int,4>> tails;
 for(int i=0;i<4;i++)for(int z:{-2,2}){array<int,4> h={};h[i]=z;tails.push_back(h);}
 for(int i=0;i<16;i++){array<int,4>h;for(int j=0;j<4;j++)h[j]=2*((i>>j)&1)-1;tails.push_back(h);}
 auto dot=[](auto v,auto w){int d=0;for(int j=0;j<4;j++)d+=v[j]*w[j];return d;};
 for(int i=0;i<24;i++){
  assert(dot(tails[i],tails[i])==4);
  for(int j=0;j<i;j++)assert(dot(tails[i],tails[j])<=2);
 }
 set<W>F0,F;
 for(int i=1;i<32;i++)for(int j=i+1;j<32;j++)for(int k=j+1;k<32;k++)if(k!=(i^j)){
  W f=0;for(int a:{0,i,j,k,i^j,i^k,j^k,i^j^k})f|=1u<<a;
  assert(popcount(f)==8);F0.insert(f);
 }
 for(W f:F0)for(int t=0;t<32;t++){
  W q=0;for(W u=f;u;u&=u-1)q|=1u<<(countr_zero(u)^t);F.insert(q);
 }
 assert(F.size()==620);
 for(W f:F)for(W q:qb)assert(!(popcount(f&q)&1));
 vector<U> old;
 for(unsigned x:C){
  U y=0;for(;x;x&=x-1)y^=cub[countr_zero(x)];old.push_back(y);
 }
 ifstream in("hyperplane_heads.txt");vector<Fibre>L;
 int which=-1;
 for(int m,k,g,a;in>>m>>k>>g>>a;){
  if(!selected_heads.count(++which))continue;
  assert(1<=m&&m<64&&(k==0||k==1)&&0<=g&&g<65536&&a==15);
  U p=0;for(int x=0;x<64;x++)if((popcount(unsigned(x&m))&1)==k)p|=U(1)<<x;
  assert(popcount(p)==32);
  // Compressing in sorted-coordinate order deletes the lowest nonzero
  // normal coordinate and is an affine parametrization of the hyperplane.
  int b=countr_zero(unsigned(m)),idx=0;
  for(U p0=p;p0;p0&=p0-1){
   unsigned x=countr_zero(p0),y=(x&((1u<<b)-1))|((x>>(b+1))<<b);
   assert(y==unsigned(idx++));
  }
  assert(none_of(L.begin(),L.end(),[&](Fibre f){return f.p==p;}));
  W v=X[g];
  for(size_t j=0;j<old.size();j++){
   W y=_pext_u64(old[j],p)^v, marked=_pext_u64(selected_sites[j],p);
   for(W q:Q){
    W d=y^q;
    int w=popcount(d);
    assert(w>=4);
    if(w==4)assert(!(d&marked));
   }
  }
  L.push_back({p,v,unsigned(a)});
 }
 assert(in.eof()&&which==125&&L.size()==selected_heads.size());
 size_t checks=0;unsigned backtracked=0;
 // A sparse support's extra parity rows are independent of the
 // rows chosen for other supports. Solve its small conjunction exactly.
 for(size_t i=0;i<B.size();i++){
  auto &s=B[i];
  vector<vector<pair<V,int>>> A;
  for(auto f:L){
   U I=f.p&s.b;W J=_pext_u64(I,f.p);int d=popcount(I);
   int maxhead=0;
   for(auto h:tails){int c=0;for(int j=0;j<4;j++)if(s.t>>j&1)c+=abs(h[j]);maxhead=max(maxhead,c);}
   int deficit=d+maxhead-12;
   if(deficit<=0)continue;
   assert(deficit<=4);
   if(deficit>2){
    assert(s.t==0&&d==16);
    for(W q:qb)assert(!(popcount(J&q)&1));
    assert((popcount(J&f.v)&1)==(s.h.s>>1&1));
    HG[i].push_back({65535,popcount(J&f.v)&1});
   }
   vector<pair<V,int>> options;
   for(W flat:F)if((flat&J)==flat){
    V row=_pext_u64(_pdep_u64(flat,f.p),s.b);
    assert(popcount(row)==8);
    Code h=s.h;int p1=1-(popcount(flat&f.v)&1);
    if(h.ins(row,p1)>=0)options.emplace_back(row,p1);
   }
   assert(!options.empty());
   A.push_back(move(options));++checks;
  }
  if(A.empty())continue;
  sort(A.begin(),A.end(),[](auto &q,auto&r){return q.size()<r.size();});
  vector<pair<V,int>> keep(A.size());
  auto pick=[&](auto&& self,size_t a,Code h)->bool{
    if(a==A.size()){s.h=h;return true;}
    for(auto [v,p]:A[a]){Code h1=h;if(h1.ins(v,p)>=0&&self(self,a+1,h1)){keep[a]={v,p};return true;}}
    ++backtracked;return false;
  };
  assert(pick(pick,0,s.h));
  HG[i].insert(HG[i].end(),keep.begin(),keep.end());
 }
 cout<<"HEAD assignment backtracks "<<backtracked<<endl;
 // On every pair of nonparallel hyperplanes the common affine 4-flat
 // has 16 sites. A parity that differs on it or any coordinate half
 // forces a disagreement, even after every quadratic phase and sign.
 U ec[6]={};
 for(int j=0;j<6;j++)for(int x=0;x<64;x++)if(x>>j&1)ec[j]|=U(1)<<x;
 for(size_t i=0;i<L.size();i++)for(size_t j=0;j<i;j++){
  auto f=L[i];auto g=L[j];U I=f.p&g.p;
  if(!I)continue;
  assert(popcount(I)==16);
  bool pass=false;
  auto test=[&](U u){
   W a=_pext_u64(u,f.p),b=_pext_u64(u,g.p);
   if((popcount(a&f.v)^popcount(b&g.v))&1){
    for(W q:qb)assert(!(popcount(a&q)&1)&&!(popcount(b&q)&1));
    pass=true;
   }
  };
  test(I);
  for(U c:ec)test(I&c);
  assert(pass);
 }
 // All dense/affine-2-flat/pole gates, written divided by sqrt(2)
 // where applicable. New-new upper bounds use squared norm 72.
 assert(32-2*4<=24&&4<=6&&2<=3);
 assert(64-4*12+8<=36&&64-4*8+4<=36&&32+4<=36&&32-4+8<=36);
 cout<<"HEADS checked "<<checks<<" sparse guards; "<<L.size()<<" hyperplanes; 24 quadratic phases"<<endl;
 return {0,4096*L.size(),4096*2*L.size(),4096*3*L.size(),2048*24*L.size()};
}

