'use strict';
// Replay of the 33D complete four-mark real-head cap. Node v22+,
// optionally pass exact certificate33_370408_explicit.json.
// Cheng–Sloane generator is from the pinned repository certificate blob
// a54ead19f3c683bf2423b6112673583f954b410e.
const assert=require('node:assert/strict');
const G=[1324400769,614502530,1214288004,508579976,1821397136,4197499040,2864709824,3529539840,2120274432,3093304320,3623933952,3428868096,518062080,2998009856,1953759232,3638034432,4278190080],H=[1,2,4,8,16,32,64,127,128,256,512,1024,1822,935,1741,1524,2048,2537,4096,8192,14839,4558,14413,8605,16384,29411,27487,24115,18321,22939,29985,27812];
if(process.argv[2]) {
  const j=JSON.parse(require('node:fs').readFileSync(process.argv[2],'utf8'));
  assert.deepEqual(j.dense_generator,G);
  assert.deepEqual(j.dense_coordinates,Array.from({length:32},(_,i)=>i));
}
function rank(A){const B=new Int32Array(32);for(let x of A)while(x){const j=31-Math.clz32(x);if(B[j])x^=B[j];else{B[j]=x;break;}}return B.filter(v=>v!==0).length}
function syndrome(x){return H.reduce((s,h,i)=>s^((x>>>i&1)?h:0),0)}
assert.equal(rank(G),17);assert.equal(rank(H),15);assert(G.every(g=>syndrome(g)===0));
const w=x=>{let n=0;for(;x;x&=x-1)n++;return n};
assert(G.every(g=>w(g)%2===0)); 
let C=[0];for (const g of G) C=C.concat(C.map(x=>x^g));
assert.equal(C.length,131072);assert.equal(C.slice(1).reduce((m,x)=>Math.min(m,w(x)),32),8);
function auto33(G){
const pc=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return Math.imul((x+(x>>>4))&0xf0f0f0f,0x01010101)>>>24};
let C=[0];for(let g of G)C=C.concat(C.map(x=>x^g));const W=C.filter(x=>pc(x)===8),hash=new Set(C);
const M=new Uint16Array(32768);for(const w0 of W){let L=[];for(let w=w0;w;w&=w-1)L.push(31-Math.clz32(w&-w));for(const x of L)for(const y of L)for(const z of L)M[x*1024+y*32+z]++;}
const M0=(x,y,z)=>M[x*1024+y*32+z],P=new Int8Array(32).fill(-1),A=Array.from({length:32},(_,i)=>i),out=[];let calls=0,limit=false;
function f(uu,V){
++calls;
if(V.length===32){if(G.every(g=>hash.has(A.reduce((s,i)=>s^((g>>>i&1)<<P[i]),0))))out.push(Array.from(P));return;}
const Cand=i=>A.filter(j=>!(uu>>>j&1)&&V.every(v=>M0(i,i,v)===M0(j,j,P[v])&&V.every(w=>M0(i,v,w)===M0(j,P[v],P[w]))));
let best=null,chosen=0;for(const i of A)if(P[i]<0){const L=Cand(i);if(!L.length)return;if(best===null||best.length>L.length){best=L;chosen=i;if(L.length===1)break;}}
for(const j of best){P[chosen]=j;V.push(chosen);f(uu|(1<<j),V);V.pop();P[chosen]=-1;if(limit)return;}
}
f(0,[]);return {out,calls,limit};}
function allmark33(P,H,k,visit){
const pack=A=>A.reduce((s,i)=>s|(1<<i),0),map=(m,p)=>p.reduce((s,v,i)=>s^((m>>>i&1)<<v),0);
const seen=new Set();function go(i,L){if(L.length<k){for(;i<=32-k+L.length;i++)go(i+1,[...L,i]);return;}
const m=pack(L);if(seen.has(m))return;const O=new Set(P.map(p=>map(m,p)));for(const z of O)seen.add(z);
const R=H.filter((v,i)=>!(m>>>i&1)),D=new Uint8Array(32768).fill(255),Q=new Int32Array(32768);let a=0,b=1;D[0]=0;
while(a<b){const v=Q[a++],n=D[v]+1;for(const h of R){const u=v^h;if(D[u]===255){Q[b++]=u;D[u]=n;}}}
if(b!==32768)throw Error('quotient not full rank');
const E=[];for(let s=0;s<32768;s++)if(D[s]===7)E.push(s);else if(D[s]>7)throw Error('radius >7');
visit(L,E,D);
}go(0,[]);return seen.size;}
function graph33(D,E){
const K=new Uint8Array(32768);for(const x of E)for(const y of E)K[x^y]++;
const V=[];for(let x=1;x<32768;x++)if(K[x]&&D[x]>=4)V.push(x);
return {K,V};
}
function nec33(D,K,V,S,n0,n1){
const V0=V.filter(s=>D[s]===6),V1=V.filter(s=>D[s]===4); if(V0.length<n0-2)return {nO:0,p:0,calls:0,samples:[]};
let p=0,nO=0,calls=0,best=0,samples=[],limit=false;const O=[0],paired=[false],target=n1===10?6:5;
function echeck(){
 nO++;const C=new Uint8Array(32768),V=[];
 for(const x of O)for(const s of S){const z=x^s;if(C[z]++===0)V.push(z);}
 const U=V.filter(z=>C[z]>=n0-5);const deg=O.map(x=>U.filter(z=>D[x^z]===7).length);if(Math.min(...deg)<target){samples.push({O:[...O],minDegree:Math.min(...deg),needed:target});return;}
 const b=U.map(z=>C[z]).sort((a,b)=>b-a).slice(0,n1).reduce((a,b)=>a+b,0);if(b<target*n0){samples.push({O:[...O],edgeBound:b,needed:target*n0});return;}
 p++; if(samples.length<20)samples.push({O:[...O],U:U.length,deg,coverage:b});
}
function solve(V){
 ++calls;if(O.length===n0){echeck();return};
 if(O.length+V.length<n0||limit)return;
 for(let j=0;j<=V.length-n0+O.length;j++){
  const u=V[j];let mate=-1,bad=false;
  for(let k=0;k<O.length;k++){let t=D[u^O[k]];if(t===4){if(mate!==-1||paired[k]){bad=true;break;}mate=k;}else if(t!==6){bad=true;break;}}
  if(bad)continue;
  const U=[];for(let k=j+1;k<V.length;k++){const v=V[k],z=v^u;
   if(!K[z]||D[z]<4||D[z]===4&&mate>=0)continue;U.push(v);}
  O.push(u);paired.push(mate!==-1);if(mate>=0)paired[mate]=true;solve(U);if(mate>=0)paired[mate]=false;O.pop();paired.pop();if(limit)return;
 }
}
solve(V0);for(const u of V1){if(limit)break;O.push(u);paired[0]=true;paired.push(true);solve(V0.filter(v=>K[u^v]&&D[u^v]===6));O.pop();paired.pop();paired[0]=false;}
return {nO,p,calls,limit,samples};
}
const aut=auto33(G.map(x=>x|0));assert.equal(aut.out.length,384);
const histogram={},cases=[];let reps=0;
const total=allmark33(aut.out,H,4,(M,S,D)=>{
  ++reps;histogram[S.length]=(histogram[S.length]||0)+1;
  assert(S.every(x=>D[x]===7)); // D is the exact syndrome walk metric.
  for(const s of M)assert.equal(D[H[s]],7); // useful independent check
  if(S.length<5)return; // already excludes the larger group's minimum incidence
  const {K,V}=graph33(D,S);
  const a=nec33(D,K,V,S,8,9),b=nec33(D,K,V,S,7,10);
  assert.equal(a.p,0);assert.equal(b.p,0);
  if(a.nO||b.nO)cases.push({M,S:S.length,a:a.samples,b:b.samples});
});
assert.equal(reps,140);assert.equal(total,35960);
assert.equal(cases.length,1);assert.deepEqual(cases[0].M,[0,1,2,5]);
assert.equal(cases[0].a.length,1);assert.equal(cases[0].b.length,7);
console.log(JSON.stringify({markSets:total,markOrbits:reps,automorphisms:aut.out.length,histogram,cases}));
