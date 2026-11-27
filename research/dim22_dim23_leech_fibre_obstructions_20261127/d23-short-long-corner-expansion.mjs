// Lower61: fixture-free integer replay of the 23D fixed-stock partial norm-6 corner obstruction.
import assert from "node:assert/strict";

const N=24, G=[0]; const base=((1<<11)|(1<<10)|(1<<6)|(1<<5)|(1<<4)|(1<<2)|1);
function pop(x){let n=0;for(let z=x;z;z&=z-1)n++;return n;}
for(let i=0;i<12;i++){let v=base<<i;v^=(pop(v)&1)<<23;G.push(...G.slice().map(x=>x^v));}
const W={};for(let u of G){let n=pop(u);(W[n]??=[]).push(u);}
if(G.length!==4096||new Set(G).size!==4096||W[8].length!==759||W[12].length!==2576||W[16].length!==759)throw Error("Golay input");


const D=W[12][0], F=(1<<24)-1-D;
const A=(F&-F)|((F^(F&-F))&-(F^(F&-F))), B=F^A;
const part=Array.from({length:24},(_,k)=>A>>k&1?0:D>>k&1?1:2);
const choose=(n,k)=>{if(k<0||k>n)return 0;let z=1;for(let i=1;i<=k;i++){z=z*(n+1-i)/i;if(!Number.isSafeInteger(z))throw Error("binomial");}return z;};
const K=new Map(), tab=(s,a,b)=>K.get([s,a,b].join())||0, odd=x=>Math.abs(x)%2;
function add(s,a,b,n){if(!n)return;if(a%2||b%4)throw Error("coordinate");let q=[s,a/2,b/4].join();K.set(q,(K.get(q)||0)+n);}
for(let i=0;i<N;i++)for(let j=i+1;j<N;j++)for(let x of [-4,4])for(let y of [-4,4])add(4,(part[i]===0)*x+(part[j]===0)*y,(part[i]===1)*x+(part[j]===1)*y,1);
for(let w of [8,12])for(let u of W[w]){
 const n=pop(A&u),m=pop(D&u),l=w-n-m;
 for(let i=0;i<=n;i++)for(let j=0;j<=m;j++)for(let k=0;k<=l;k++){
  const a=2*(n-2*i),b=2*(m-2*j),t=choose(n,i)*choose(m,j)*choose(l,k);
  if(!odd(i+j+k)){add(w===8?4:6,a,b,t);continue;}
  if(w===12)continue;
  for(let z of [-4,4]){add(6,a+z,b,t*(2-n));add(6,a,b+z,t*(12-m));add(6,a,b,t*(10-l));}
 }
}
for(let u of G){
 const n=pop(A&u),m=pop(D&u),l=pop(B&u);
 const t=[[2-n,n],[12-m,m],[10-l,l]],a=2-2*n,b=12-2*m;
 for(let k=0;k<3;k++)for(let j=0;j<2;j++){
  const v=j?-4:4,x=k===0?v:0,y=k===1?v:0;
  add(4,a-x,b-y,t[k][j]); add(6,a+x,b+y,t[k][j]);
 }
 const f=t.flat();
 for(let i=0;i<=3;i++)for(let j=0;j<=3-i;j++)for(let k=0;k<=3-i-j;k++)for(let r=0;r<=3-i-j-k;r++)for(let s=0;s<=3-i-j-k-r;s++){
  const q=3-i-j-k-r-s;
  add(6,a-4*(i-j),b-4*(k-r),[i,j,k,r,s,q].reduce((z,x,g)=>z*choose(f[g],x),1));
 }
}
const T=[4,6].map(s=>[...K].filter(([k,v])=>k.startsWith(s+",")).reduce((z,[k,v])=>z+v,0));
if(T.join()!=="196560,16773120")throw Error("shell count");
function mod(a,p){return ((a%p+p+Math.floor(p/2))%p)-Math.floor(p/2);}
function q(a,b){return 3*a*a+2*b*b;}
const M=new Map();for(let [key,n]of K){
 const [s,a,b]=key.split(",").map(Number),r=12*s-q(a,b);
 if(r<=0||r>48)continue;
 const A=mod(a,4),B=mod(b,6),S=(r+q(A,B))/12;
 if(S!==4&&S!==6)throw Error("even saturation");
 const k=[r,A,B].join(),o=M.get(k);if(o&&o.n!==n)throw Error("duplicate");
 M.set(k,{s:S,a:A,b:B,r,n});
}
const L=[...M.values()].sort((x,y)=>x.s-y.s||x.a-y.a||x.b-y.b);
function pp(i,j){const a=L[i],b=L[j],x=mod(a.a-b.a,4),y=mod(a.b-b.b,6),Q=q(x,y),s=tab(4,x,y)?4:tab(6,x,y)?6:8;return {u:a.r+b.r-12*s+Q,rs:a.r*b.r,s,Q};}

const profiles=new Map(L.map((f,i)=>[[f.s,f.a,f.b].join(),i]));
const data=L.map(x=>new Int8Array(x.n*24));
const cnt=new Int32Array(L.length);
let x=new Int8Array(24), Sa=[],Sd=[];for(let i=0;i<24;i++){if(part[i]===0)Sa.push(i);else if(part[i]===1)Sd.push(i);}
function put(s,aa,bb){
 let a=aa,b=bb;
 if(a===undefined){a=0;b=0;for(let i of Sa)a+=x[i];for(let i of Sd)b+=x[i];}
 const f=profiles.get([s,a/2,b/4].join());if(f===undefined)return;
 let ix=cnt[f]++*24;
 for(let i=0;i<24;i++)data[f][ix+i]=3*x[i]-(part[i]===0?3*a/2:part[i]===1?b/4:0);
}
for(let i=0;i<24;i++)for(let j=i+1;j<24;j++)for(let a of [-4,4])for(let b of [-4,4]){
 x[i]=a;x[j]=b;put(4);x[i]=x[j]=0;
}
for(let w of [8,12])for(let u of W[w]){
 if(w===12&&(pop(A&u)<2||pop(D&u)<6))continue;
 let S=[],O=[],Q=[];for(let i=0;i<24;i++)((u>>i&1)?S:O).push(i);for(let i of S)Q.push(part[i]);
 for(let z=0;z<1<<w;z++){
  let a=0,b=0;for(let i=0;i<w;i++){let k=S[i],v=z>>i&1?-2:2;x[k]=v;if(Q[i]===0)a+=v;else if(Q[i]===1)b+=v;}
  if(!(pop(z)&1)){put(w===8?4:6,a,b);continue;}
  if(w===12||Math.abs(a+4)>4||Math.abs(b+12)>4)continue;
  for(let i of O)for(let q of [-4,4]){
   let A=a+(part[i]===0?q:0),B=b+(part[i]===1?q:0);
   if(A===-4&&B===-12){x[i]=q;put(6,A,B);x[i]=0;}
  }
 }
 for(let i of S)x[i]=0;
}
for(let u of G){
 let a=0,b=0;for(let i=0;i<24;i++){x[i]=u>>i&1?-1:1;if(part[i]===0)a+=x[i];if(part[i]===1)b+=x[i];}
 for(let i=0;i<24;i++){
  let v=x[i],da=part[i]===0?4*v:0,db=part[i]===1?4*v:0;
  x[i]=-3*v;put(4,a-da,b-db);x[i]=v;
  if(a+da===-4&&b+db===-12){x[i]=5*v;put(6,-4,-12);x[i]=v;}
 }
 for(let i=0;i<24;i++)for(let j=i+1;j<24;j++)for(let k=j+1;k<24;k++){
  let A=a-4*((part[i]===0)*x[i]+(part[j]===0)*x[j]+(part[k]===0)*x[k]),
      B=b-4*((part[i]===1)*x[i]+(part[j]===1)*x[j]+(part[k]===1)*x[k]);
  if(A!==-4||B!==-12)continue;
  let a0=x[i],b0=x[j],c0=x[k];x[i]*=-3;x[j]*=-3;x[k]*=-3;put(6,A,B);x[i]=a0;x[j]=b0;x[k]=c0;
 }
}
if(cnt.join()!==L.map(f=>f.n).join())throw Error("population mismatch "+JSON.stringify(L.map((x,i)=>[i,x.n,cnt[i]]).filter(x=>x[1]!==x[2])));
function dp(i,j,a,b){const x=data[i],y=data[j];let d=0;for(let k=0;k<24;k++)d+=x[24*a+k]*y[24*b+k];return d;}

const raw=(i,j)=>Array.from({length:24},(_,k)=>(data[i][24*j+k]+(part[k]===0?3*L[i].a:part[k]===1?L[i].b:0))/3);
const enc=x=>x.map(v=>(v+8).toString(16)).join("");
const find=(s,a,b)=>L.findIndex(f=>[s,a,b].join()===[f.s,f.a,f.b].join());
const rows=(s,a,b)=>{let i=find(s,a,b);return Array.from({length:L[i].n},(_,j)=>raw(i,j));};
const A0=rows(6,-2,-3), lowParents=rows(4,-2,-2),highParents=rows(4,-2,0);
const oldLow=new Set(rows(4,0,-1).map(enc)), oldHigh=new Set(rows(4,0,-3).map(enc));
assert.equal(A0.length,9072);
assert.equal(lowParents.length,162);assert.equal(highParents.length,2004);
assert.equal(oldLow.size,23004);assert.equal(oldHigh.size,324);
const Gset=new Set(G);
function member(x,s,a,b){
 const rem=x[0]&1;let r=0,p=0,q=0,mask=0,sm=0;
 for(let k=0;k<24;k++){
  if((x[k]&1)!==rem)return false;
  sm+=x[k];r+=x[k]*x[k];
  if((x[k]-rem)&3)mask|=1<<k;
  if(part[k]===0)p+=x[k];if(part[k]===1)q+=x[k];
 }
 return ((sm-4*rem)&7)===0&&r===8*s&&p===2*a&&q===4*b&&Gset.has(mask);
}
for(const x of A0)assert(member(x,6,-2,-3));
assert.equal(new Set(A0.map(enc)).size,A0.length);
assert(lowParents.every(z=>member(z,4,-2,-2)));assert(highParents.every(z=>member(z,4,-2,0)));
assert.equal(new Set(lowParents.map(enc)).size,lowParents.length);
assert.equal(new Set(highParents.map(enc)).size,highParents.length);
const degree=new Map();
function rawdot(x,z){let t=0;for(let k=0;k<24;k++)t+=x[k]*z[k];return t;}
for(const x of A0){
 let l=0,h=0;
 for(const z of lowParents)if(rawdot(x,z)===24){
  const y=enc(x.map((w,k)=>w-z[k]));assert(oldLow.has(y));l++;degree.set(y,(degree.get(y)||0)+1);
 }
 for(const z of highParents)if(rawdot(x,z)===24){const y=enc(x.map((w,k)=>w-z[k]));assert(oldHigh.has(y));h++;}
 assert.equal(l,30);assert.equal(h,12);
}
assert.equal(degree.size,22680);assert([...degree.values()].every(v=>v===12));
const eqOwners=rows(4,-2,-3);
assert.equal(eqOwners.length,2);
for(let k=0;k<24;k++)assert.equal(eqOwners[0][k]+eqOwners[1][k],part[k]===0?-4:part[k]===1?-2:0);
assert(25*63>36*36); // at T=sqrt(63), (5/2+T/6)^2-(7/2+T²/6)>0
assert(63>15); // at T>=sqrt(63), T²/12+3T/2-5/4>0
console.log(JSON.stringify({candidateCount:A0.length,lowCandidateDegree:30,lowOldDegree:12,
 oldLowWithPositiveDegree:degree.size, highCandidateDegree:12, equatorOwnerStocks:2},null,2));
