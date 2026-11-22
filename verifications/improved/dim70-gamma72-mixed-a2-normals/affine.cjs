
const SM=15281788354560n,S8=6218175600n;
function comb(n,k){let x=1n;for(let t=1;t<=k;t++)x=x*BigInt(n+1-t)/BigInt(t);return x}
const E4=[1n,240n,2160n,6720n,17520n,30240n,60480n];
const E6=[1n,-504n,-16632n,-122976n,-532728n,-1575504n,-4058208n];
const dl=[0n,1n,-24n,252n,-1472n,4830n,-6048n];
function mul(a,b,l=7){return Array.from({length:l},(_,i)=>a.slice(0,i+1).reduce((s,x,j)=>s+x*(b[i-j]||0n),0n))}
function pow(a,k,l=7){let z=[1n];for(let i=0;i<k;i++)z=mul(z,a,l);return z}
const theta=pow(E4,9).map((x,i)=>x-2160n*mul(pow(E4,6),dl)[i]+965520n*mul(pow(E4,3),pow(dl,2))[i]-27302400n*pow(dl,3)[i]);
if(theta[4]!==S8||theta[5]!==SM)throw Error("theta");
const vs=[-2,-1,0,1,2];
const gcd=(a,b)=>b?gcd(b,a%b):(a<0n?-a:a);
function pol(r,t,j,s) {
 const k=r+t;let nu=1n,de=1n;let prev=new Map([[r*64+t,1n]]),C=[[r,t,1n,1n]];
 for(let l=1;l<=k/2;l++){
  nu*=BigInt(-s);de*=BigInt(2*l*(2*k-2*l+70));let g=gcd(nu,de);nu/=g;de/=g;
  let n=new Map;
  const add=(a,b,h)=>{const q=a*64+b;n.set(q,(n.get(q)||0n)+h)};
  for(const [a_b,h] of prev){
   const a=Math.floor(a_b/64),b=a_b%64;
   if(a>=2)add(a-2,b,h*BigInt(8*a*(a-1)));
   if(a&&b)add(a-1,b-1,h*BigInt(2*j*a*b));
   if(b>=2)add(a,b-2,h*BigInt(8*b*(b-1)));
  }
  for(const [i,h] of n)C.push([Math.floor(i/64),i%64,nu*h,de]);
  prev=n;
 }
 const den=C.reduce((a,[,,,b])=>a/gcd(a,b)*b,1n);
 return [den,C.map(([a,b,v,w])=>[a,b,v*(den/w)])];
}
const powers=new Map(Array.from({length:21},(_,i)=>[i-10,Array.from({length:45},(_,k)=>BigInt(i-10)**BigInt(k))]));
function evalPoly(t,a,b){const p=powers.get(a),q=powers.get(b);return t.reduce((s,[r,h,v])=>s+v*p[r]*q[h],0n)}
function moment(a,b,j,s){
 const k=(a+b)/2;let w=0n,odd=n=>{let v=1n;for(let t=n;t>0;t-=2)v*=BigInt(t);return v};
 for(let l=0;l<=Math.min(a,b);l++)if((a-l)%2===0&&(b-l)%2===0){
  let f=1n;for(let t=2;t<=l;t++)f*=BigInt(t);
  w+=comb(a,l)*comb(b,l)*f*odd(a-l-1)*odd(b-l-1)*8n**BigInt(k-l)*BigInt(j)**BigInt(l);
 }
 let den=1n;for(let h=0;h<k;h++)den*=BigInt(72+2*h);
 return [theta[s/2]*BigInt(s)**BigInt(k)*w,den];
}
function relations(k,ns){
 const basis=[];
 for(let j=0;j<=2;j++)for(let b=0;b<=1;b++){
  let ap=k-12-6*b-12*j;
  if(ap>=0&&ap%4===0)basis.push(mul(mul(mul([1n,-96n,4464n],pow(E4,ap/4,3),3),pow(E6,b,3),3),pow(dl,j,3),3).slice(0,ns));
 }
 if(!basis.length)return Array.from({length:ns},(_,i)=>Array.from({length:ns},(_,j)=>BigInt(i===j)));
 if(basis.length===1)return Array.from({length:ns-1},(_,i)=>Array.from({length:ns},(_,j)=>j===0?-basis[0][i+1]:BigInt(j===i+1)));
 if(basis.length===ns)return [];
 const [a,b]=basis;
 return [[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]];
}
function system(j,ns){
 const S=ns===3?[8,10,12]:[8,10],P={},E={},all=[];let K=[],idx=new Map;
 const q=(s,a,b,x,y,z)=>8*(x*x+y*y)+s*z*z+2*(j*x*y+a*x*z+b*y*z);
 const key=(s,a,b)=>[s,a<0||a===0&&b<0?-a:a,a<0||a===0&&b<0?-b:b].join(",");
 for(const s of S){
  E[s]=[];
  for(const x of vs)for(const y of vs)if(q(0,0,0,x,y,0)===s)E[s].push([8*x+j*y,j*x+8*y]);
  P[s]=[];
  for(let a=-(s)/2;a<=(s)/2;a++)for(let b=-(s)/2;b<=(s)/2;b++){
   if(E[s].some(([u,v])=>a===u&&b===v))continue;
   if(64*s+2*j*a*b-s*j*j-8*(a*a+b*b)<0)continue;
   if(vs.some(x=>vs.some(y=>vs.some(z=>(x||y)&&z&&q(s,a,b,x,y,z)>0&&q(s,a,b,x,y,z)<8))))continue;
   if(vs.some(x=>vs.some(y=>(x||y)&&q(s,a,b,x,y,1)===0)))continue;
   P[s].push([a,b]);const t=key(s,a,b);
   if(!idx.has(t)){idx.set(t,K.length);K.push(t);}
  }
 }
 function row(){return Array(K.length+1).fill(0n)}
 function add(R,s,a,b,w){const k=idx.get(key(s,a,b));if(k!==undefined)R[k]+=w;else if(E[s].some(([u,v])=>a===u&&b===v))R[K.length]-=w;}
 function push(R){if(R.some(x=>x))all.push(R)}
 for(const s of S){
  for(let k=0;k<=10;k+=2)for(let a=0;a<=k;a++){
   const b=k-a,[v,d]=moment(a,b,j,s),R=row();R[K.length]=v;
   for(const [x,y] of P[s].concat(E[s]))add(R,s,x,y,d*powers.get(x)[a]*powers.get(y)[b]);push(R);
  }
 }
 for(let k=12;k<=(ns===3?38:22);k+=2)for(const rel of relations(k,ns)){
  for(let a=0;a<=k;a++){
   const b=k-a,R=row();
   const polys=S.map(s=>pol(a,b,j,s));
   const de=polys.reduce((a,[b])=>a/gcd(a,b)*b,1n);
   for(const [i,s] of S.entries())if(rel[i]){
    for(const [x,y] of P[s].concat(E[s])) add(R,s,x,y,rel[i]*(de/polys[i][0])*evalPoly(polys[i][1],x,y));
   }
   push(R);
  }
 }
 for(const s of S)for(const [a,b] of P[s].concat(E[s]))for(const x of vs)for(const y of vs)if(x||y){
  const t=q(s,a,b,x,y,1);if(!S.includes(t))continue;
  const c=a+8*x+j*y,d=b+j*x+8*y;
  const R=row();add(R,s,a,b,1n);add(R,t,c,d,-1n);push(R);
 }
 return {K,P,E,all};
}
function solveMod(all,N,p){
 const B=BigInt(p),mod=x=>((x%p)+p)%p;
 const inv=x=>{let e=p-2,r=1;while(e){if(e%2)r=r*x%p;x=x*x%p;e=Math.floor(e/2)}return r};
 const A=all.map(r=>r.map(x=>mod(Number(x%B))));let rank=0;const cols=[];
 for(let c=0;c<N;c++){
  const t=A.findIndex((r,i)=>i>=rank&&r[c]);if(t<0)continue;
  [A[t],A[rank]]=[A[rank],A[t]];
  const v=inv(A[rank][c]);
  for(let i=c;i<=N;i++)A[rank][i]=mod(A[rank][i]*v);
  for(let i=0;i<A.length;i++)if(i!==rank&&A[i][c]){
   const f=A[i][c];
   for(let k=c;k<=N;k++)A[i][k]=mod(A[i][k]-f*A[rank][k]);
  }
  cols.push(c);rank++;
 }
 if(A.slice(rank).some(r=>r[N]))throw Error("inconsistent mod "+p);
 const single=cols.filter((c,i)=>A[i].slice(0,N).filter(x=>x).length===1);
 const free=Array.from({length:N},(_,i)=>i).filter(i=>!cols.includes(i));
 const floating=free.concat(cols.filter((c,i)=>free.some(k=>A[i][k])));
 return {rank, cols,single,floating, val:cols.map((c,i)=>A[i][N]), rows:A.slice(0,rank)};
}

const axis=[-2,-1,0,1,2],C=[];
const gram=(x,y)=>2*x[0]*y[0]-x[0]*y[1]-x[1]*y[0]+2*x[1]*y[1]; // multiply A inverse by 12
function possible(s,x,y){
 const g=gram([x,y],[x,y]);
 if(g>12*s)return false;
 for(const a of axis)for(const b of axis)for(const c of axis)if(c&&(a||b)){
  const q=8*(a*a+b*b+a*b)+s*c*c+2*c*(a*x+b*y);
  if(q>0&&q<8)return false;
 }
 return true;
}
for(let a=-4;a<=4;a++)for(let b=-4;b<=4;b++){
 if(gram([a,b],[a,b])<96&&possible(8,a,b))C.push([a,b]);
}
const naive={},hit={};
for(const a of C)for(const b of C){
 const k=a+';'+b, A=96-gram(a,a),B=96-gram(b,b),ab=gram(a,b);
 const yes=(dot)=>12*dot-ab<=0 || 4*(12*dot-ab)**2<=A*B;
 naive[k]=yes(4);
 let dot=4;while(dot>-8&&!possible(16-2*dot,a[0]-b[0],a[1]-b[1]))dot--;
 hit[k]=yes(dot);
}
function verify(data, xInput){
 const {K,all}=system(4,3), N=K.length, maps=new Map(data.K.map((key,i)=>[key,i]));
 const affine=data.Z.map(v=>v.map(BigInt));
 if(data.den!=="1"||K.length!==141||solveMod(all,N,49999).rank!==140)throw Error("rank");
 const Z=K.map(key=>affine[maps.get(key)]);
 const xcol=Z[K.indexOf("8,4,4")];
 if(xcol[0]!==1n||xcol[1]!==0n)throw Error("corner");
 for(const row of all){
  let u=0n,v=-row[N];
  for(let i=0;i<N;i++)if(row[i]){u+=row[i]*Z[i][0];v+=row[i]*Z[i][1];}
  if(u||v)throw Error("affine relation");
 }
 const ceil=(x,y)=>{const q=x/y,r=x%y;return r>0n?q+1n:q};
 const floor=(x,y)=>-ceil(-x,y);
 let lo=0n,hi=6218175600n;
 for(const [c,b] of Z){
   if(c>0n){const l=ceil(-b,c);if(l>lo)lo=l;}
   if(c<0n){const h=floor(b,-c);if(h<hi)hi=h;}
   if(!c&&b<0n)throw Error("infeasible");
 }
 if(lo!==0n||hi!==8091n)throw Error("interval");
 if(Object.keys(hit).filter(k=>hit[k]&&!naive[k]).length!==1074)throw Error("edges");
 const norm8Keys=new Map(K.map((key,i)=>[key,i]));
 const coeff=(a,b)=>{
  if(a<0||a===0&&b<0){a=-a;b=-b;}
  return Z[norm8Keys.get(["8",a,b].join(","))];
 };
 const chosen=[[-1,-1],[-1,0],[-2,-1]];
 function optimal(edge,V,expected){
   const w=q=>{const [a,b]=coeff(...q);return a*V+b;};
   const L=C.map(q=>[q,w(q)]).sort((a,b)=>a[1]>b[1]?-1:a[1]<b[1]?1:0);
   let best=-1n;
   function go(s,L){
    if(L.reduce((a,b)=>a+b[1],s)<=best)return;
    if(s>best)best=s;
    for(let i=0;i<L.length;i++)if(L[i][1]>0n)go(s+L[i][1],L.slice(i+1).filter(y=>edge[L[i][0]+";"+y[0]]));
   }
   go(0n,L);
   if(best!==expected||chosen.reduce((s,q)=>s+w(q),0n)!==expected||chosen.some(a=>chosen.some(b=>!edge[a+";"+b])))throw Error("weighted clique");
 }
 optimal(hit,lo,1361180296n);
 optimal(hit,hi,1360986112n);
 for(const [x,values] of [[1795n,[610915680n,610915680n,139305856n]],[1579n,[610923456n,610923456n,139295488n]]]){
   if(chosen.some((p,i)=>{const [a,b]=coeff(...p);return a*x+b!==values[i]}))throw Error("publication counts");
 }
 const x=BigInt(xInput);
 if(x<lo||x>hi)throw Error("corner outside proven range");
 const counts=chosen.map(p=>{const [a,b]=coeff(...p);return a*x+b;});
 const lower=counts.reduce((a,b)=>a+b,0n);
 return JSON.stringify({counts:counts.map(String),lower:String(lower)});
}
if(typeof require!=="undefined"&&require.main===module){
 const data=JSON.parse(require("fs").readFileSync(process.argv[2],"utf8"));
 console.log(verify(data,process.argv[3]));
}
