// Exact Li 21D local barrier through sixteen arbitrary old-hole removals.

const lineLists=[[4,5,6,8,10],[2,3,6,7,15],[1,4,9,14,15],[3,9,10,11,16],[1,6,12,13,16],[7,8,9,12,17],[4,7,11,13,18],[2,10,12,14,18],[1,3,5,17,18],[2,5,9,13,19],[6,11,14,17,19],[8,15,16,18,19],[1,2,8,11,20],[5,7,14,16,20],[10,13,15,17,20],[3,4,12,19,20],[3,8,13,14,21],[5,11,12,15,21],[2,4,16,17,21],[1,7,10,19,21],[6,9,18,20,21]];
const tags=[11,13,8,15,3,1,12,5,10,4,7,1,2,3,2,1,9,6,14,3,2];
const lineMasks=lineLists.map(x=>x.reduce((s,j)=>s|(1<<(j-1)),0));
const H=[0], eta=[0], M=new Map([[0,0]]);
for(let i=0;i<21;i++){let x=lineMasks[i];if(!M.has(x)){let n=H.length;for(let j=0;j<n;j++){M.set(H[j]^x,j+n);H.push(H[j]^x);eta.push(eta[j]^tags[i]);}}}
if(H.length!==1024||lineMasks.some((x,i)=>eta[M.get(x)]!==tags[i]))throw Error("geometry rank/tag");
const L=lineMasks.map(x=>M.get(x)), holes=H.flatMap((x,i)=>eta[i]===0?[i]:[]);
const LC=L.filter((x,i)=>tags.filter(y=>y===tags[i]).length===1);
if(LC.length!==12||holes.length!==64)throw Error("phase");
const state={};
const pc=x=>{x-=(x>>>1)&0x55555555;x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0xf0f0f0f)*0x1010101)>>>24};

const ven=holes.filter(i=>pc(H[i])&1),VI=new Map(ven.map((x,i)=>[x,i]));
if(ven.length!==32)throw Error('half');
const C=[];
for(let x=0;x<1024;x++)if(eta[x]&&(pc(H[x])&1)){
 let c=new Uint8Array(32);
 for(let p of LC)for(let t of L){let i=VI.get(x^p^t);if(i!==undefined)c[i]++;}
 C.push({typ:eta[x]<=3?1:0,x,c});
}
const P=[Array(33).fill(0),Array(33).fill(0),Array(33).fill(0)];
for(let i=0;i<C.length;i++)for(let j=0;j<i;j++){
 let x=C[i],y=C[j],w=P[x.typ+y.typ],cnt=[0,0,0,0,0], k=0,z=0;
 for(let t=0;t<32;t++)cnt[x.c[t]+y.c[t]]++;
 for(let v=4;v>=0;v--)for(let t=0;t<cnt[v];t++){z+=v;k++;if(z>w[k])w[k]=z;}
}

const T3=Array.from({length:4},()=>new Uint8Array(33));
const repr=[];for(let tag=1;tag<16;tag++)repr.push(C.findIndex(x=>eta[x.x]===tag));
if(repr.some(x=>x<0))throw Error('missing tag');
let triples=0;
const counts=new Uint8Array(7);
for(let ix of repr){const x=C[ix];for(let j=0;j<C.length;j++)if(j!==ix){const y=C[j];for(let k=0;k<j;k++)if(k!==ix){
 const z=C[k], w=T3[x.typ+y.typ+z.typ];counts.fill(0);triples++;
 for(let a=0;a<32;a++)counts[x.c[a]+y.c[a]+z.c[a]]++;
 for(let v=6,z=0,i=0;v>=0;v--)for(let t=0;t<counts[v];t++){z+=v;i++;if(z>w[i])w[i]=z;}
}}}

const smax=k=>Math.min(2*k,k+4,17), tmax=k=>Math.min(2*k,12);
function possible(k){let all=0,first=null,distribution={},arr=[];
 for(let j=1;j<k;j++)for(let kr=0;kr<=k;kr++){
  let ks=[kr,k-kr];
  for(let a=0;a<=j;a++)for(let b=0;b<=j-a;b++)for(let c=0;c<=j-a-b;c++){
   let d=j-a-b-c;
   let sa=[a,c], ta=[b,d];
   if(!ks.every((p,i)=>p<=sa[i^1]+3*ta[i^1]&&12*p<=smax(p)*sa[i]+tmax(p)*ta[i]+sa[i^1]))continue;
   all++;arr.push({k,j,ks,sa,ta});first??={k,j,ks,sa,ta};let f=ks.join();distribution[f]=(distribution[f]||0)+1;
 }}
 return {all,first,distribution,arr};
}

function pcaps(k,s,t){
 const ps=P[0][k],pm=P[1][k],pt=P[2][k];
 let v=0;
 for(let a=0;a<=(s?smax(k):0);a++)for(let b=0;b<=(t?tmax(k):0);b++){
  if(s&&t&&a+b>pm)continue;
  v=Math.max(v,(s?a+(s-1)*Math.min(a,ps-a):0)+(t?b+(t-1)*Math.min(b,pt-b):0));
 }return v;
}
function refin(k){let aa=possible(k).arr.filter(({ks,sa,ta})=>ks.every((p,i)=>12*p<=pcaps(p,sa[i],ta[i])+sa[i^1]));return aa;}

const tpc=new Map(), seqMemo=new Map();
function seqs(k,s,ty){
 let kw=[k,s,ty].join(),old=seqMemo.get(kw);if(old)return old;
 let max=ty?tmax(k):smax(k), n=Math.min(s,3), r=[],a=[];
 function F(j){
  if(j===n){let z=a.reduce((u,v)=>u+v,0)+(s>3?(s-3)*a[2]:0);r.push([[0,...a.reduce((u,v)=>{u.push((u.at(-1)??0)+v);return u;},[])],z]);return;}
  for(let i=j?a[j-1]:max;i>=0;i--){a.push(i);if((j<1||a[0]+a[1]<=P[2*ty][k])&&(j<2||a[0]+a[1]+a[2]<=T[3*ty][k]))F(j+1);a.pop();}
 }
 F(0);seqMemo.set(kw,r);return r;
}
function tcap(k,s,t){
 const key=[k,s,t].join();if(tpc.has(key))return tpc.get(key);
 let b=0;
 for(let [A,v] of seqs(k,s,0))for(let [B,w] of seqs(k,t,1)){
 if(v+w<=b||(s&&t&&A[1]+B[1]>P[1][k])||(s>=2&&t&&A[2]+B[1]>T[1][k])||(s&&t>=2&&A[1]+B[2]>T[2][k]))continue;b=v+w;
 }
 tpc.set(key,b);return b;
}
function tr(k){return refin(k).filter(({ks,sa,ta})=>ks.every((p,i)=>12*p<=tcap(p,sa[i],ta[i])+sa[i^1]));}

const T=T3.map(x=>Array.from(x));
function check(z,m){if(!z)throw Error(m);}
for(const v of C){const cnt=[0,0,0];for(let u of v.c){check(u>=0&&u<=2,'coverage multiplicity');cnt[u]++;}
 check(v.typ===0?(cnt[1]===9&&cnt[2]===4):(cnt[1]===0&&cnt[2]===6),'literal coverage profiles');}
check(C.length===480&&C.filter(x=>x.typ===0).length===384&&triples===1717215,'coverage/rank counts');
const base=Array.from({length:18},(_,k)=>possible(k).arr.length);
const pair=Array.from({length:18},(_,k)=>refin(k).length);
const mixed=Array.from({length:18},(_,k)=>tr(k).length);
check(base.slice(1,14).every(i=>i===0)&&base[14]===2&&base[15]===18&&base[16]===110,'unary cases');
check(pair.slice(1,16).every(i=>i===0)&&pair[16]===36,'pair cases');
check(mixed.slice(1,17).every(i=>i===0)&&mixed[17]===12,'triple cases');
console.log(JSON.stringify({unary:base,pair,mixed,tripleScan:triples,P,T}));
