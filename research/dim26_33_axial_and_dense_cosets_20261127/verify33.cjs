const pc=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return Math.imul((x+(x>>>4))&0xf0f0f0f,0x01010101)>>>24};
const axes=Array.from({length:32},(_,i)=>i);
function quot(G,marks){
 let A=G.map(x=>x|0),zero=marks.reduce((s,x)=>s|(1<<x),0);
 for(const i of marks){
  const p=A.findIndex(x=>x>>>i&1);
  if(p===-1)throw Error('marks not independent');
  const v=A[p];
  A=A.filter((x,j)=>j!==p).map(x=>(x>>>i&1)?x^v:x);
 }
 const B=Array(32).fill(0);
 for(let x of A)while(x){const j=31-Math.clz32(x);if(B[j])x^=B[j];else {B[j]=x;break;}}
 if(B.filter(x=>x).length!==17-marks.length)throw Error('dim');
 const rest=axes.filter(i=>!marks.includes(i)),free=rest.filter(i=>!B[i]);
 if(free.length!==15)throw Error('codim');
 const map=Array(32).fill(0); for(let j=0;j<32;j++)if(j in B&&B[j])for(let i=j+1;i<32;i++)if(B[i]>>>j&1)B[i]^=B[j];
 for(const i of rest){let z=B[i]||1<<i,w=0;
  for(let j=0;j<15;j++)if(z>>>free[j]&1)w|=1<<j;
  map[i]=w;
 }
 const D=new Uint8Array(32768).fill(255),Q=new Int32Array(32768),L=rest.map(x=>map[x]);let a=0,b=1;
 D[0]=0;Q[0]=0;
 while(a<b){const v=Q[a++],n=D[v]+1;for(let j=0;j<L.length;j++){const u=v^L[j];if(D[u]===255){D[u]=n;Q[b++]=u;}}}
 if(b!==32768)throw Error('incomplete');
 function codeSyndrome(w){let s=0;for(const i of rest)if(w>>>i&1)s^=map[i];return s;}
 const old=new Set([0]);for(const v of G){const s=codeSyndrome(v);for(const a of [...old])old.add(a^s);}
 if(old.size!==1<<marks.length)throw Error('punctured quotient');
 return {D,map,old:[...old].sort((a,b)=>a-b),rest,free};
}

function leaderLists(map,cap=7){
 const L=Array.from({length:32768},()=>[]);
 function go(j,left,s,m){
  if(left===0){L[s].push(m|0);return;}
  for(;j<=32-left;j++)go(j+1,left-1,s^map[j],m|(1<<j));
 }
 for(let k=0;k<=cap;k++)go(0,k,0,0);
 return L;
}
function hittable(L,cap,mask=0){
 let o=L.find(z=>(z&mask)===0);
 if(o===undefined)return mask;
 if(!cap)return null;
 for(let x=o;x;x^=(x&-x)){const t=hittable(L,cap-1,mask|(x&-x));if(t!==null)return t;}
 return null;
}

function puncture(map,marks){
 const B=new Int32Array(15);
 for(const i of marks){let x=map[i];while(x){const j=31-Math.clz32(x);if(B[j])x^=B[j];else{B[j]=x;break;}}}
 const free=[];for(let j=0;j<15;j++)if(!B[j])free.push(j);
 if(free.length!==15-marks.length)throw Error('marked syndrome rank');
 function red(x){for(let j=14;j>=0;j--)if(B[j]&&(x>>j&1))x^=B[j];let s=0;for(let j=0;j<free.length;j++)if(x>>free[j]&1)s|=1<<j;return s;}
 const A=[...new Set(map.map(red))].filter(x=>x!==0),D=new Uint8Array(1<<free.length).fill(255),Q=new Int32Array(D.length);
 let a=0,b=1,max=0;Q[0]=0;D[0]=0;
 while(a<b){const v=Q[a++],n=D[v]+1;for(let j=0;j<A.length;j++){const u=v^A[j];if(D[u]===255){D[u]=n;Q[b++]=u;if(max<n)max=n;}}}
 if(b!==Q.length)throw Error('span');
 return {D,max,map:A,free,B};
}
function sublists(n,k,f){const A=[];function go(i){if(A.length===k){f(A);return;}for(;i<=n-k+A.length;i++){A.push(i);go(i+1);A.pop();}}go(0);}

function r2(H){
const L=leaderLists(H,6),ret={n:0,M:{},singles:0};
for(let s=0;s<32768;s++){const V=L[s],L0=V[0];if(!V.length)throw Error('cover>6 initially');
 function add(a,b){if(a>b)[a,b]=[b,a];const o=ret.M[[a,b]]||(ret.M[[a,b]]=[]);if(o.at(-1)!==s)o.push(s);}
 for(let x=L0;x;x^=x&-x){const a=x&-x,ai=31-Math.clz32(a),v=V.find(m=>!(m&a));
  if(v===undefined){ret.singles++;for(let i=0;i<32;i++)if(i!==ai)add(i,ai);}
  else for(let y=v;y;y^=y&-y){const b=y&-y;if(V.every(m=>m&(a|b)))add(ai,31-Math.clz32(b));}
 }
}
let same=true;
for(let a=0;a<32;a++)for(let b=a+1;b<32;b++)same&&=ret.M[[a,b]]?.length===2&&ret.M[[a,b]].includes(H[a])&&ret.M[[a,b]].includes(H[b]);
return {...ret,same};
}
const assert=require('node:assert/strict');
const fs=require('node:fs');const C=JSON.parse(fs.readFileSync(process.argv[2],'utf8'));const H=quot(C.dense_generator,[]).map;
const PIN=[1,2,4,8,16,32,64,127,128,256,512,1024,1822,935,1741,1524,2048,2537,4096,8192,14839,4558,14413,8605,16384,29411,27487,24115,18321,22939,29985,27812];
assert.deepEqual(H,PIN);assert(C.dense_generator.every(g=>pc(g)%2===0));
const L=leaderLists(H,7);
assert.equal(L.reduce((s,y)=>s+y.length,0),4514873);assert(L.every(y=>hittable(y,4)===null));
const out=r2(H);assert(out.same);assert.equal(out.singles,32);assert.equal(Object.keys(out.M).length,496);
const P={};sublists(32,4,S=>{const p=puncture(H,S);P[p.max]=(P[p.max]||0)+1;});
assert.deepEqual(P,{4:35960});
console.log('33D: all 32768 syndromes survive <=4 deletions at distance 7; all 496 double-deletions have exactly their 2 prescribed radius-7 syndromes; all 35960 full-code four-punctures have covering radius 4.');