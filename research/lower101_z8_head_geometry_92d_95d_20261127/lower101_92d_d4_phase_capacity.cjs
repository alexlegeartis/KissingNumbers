'use strict';
// Exact local combinatorics for the 92D retained-D4 conservative envelope.
const assert=require('node:assert/strict');
function setup101(){
 const ax=Array.from({length:4},(_,i)=>i),sg=[1,-1],P=[];
 for(let i=0;i<4;i++)for(let j=i+1;j<4;j++)for(let a of sg)for(let b of sg)P.push(ax.map(k=>k===i?2*a:k===j?2*b:0));
 const Q=[];for(let i=0;i<4;i++)for(let a of sg)Q.push(ax.map(k=>i===k?2*a:0));for(let b=0;b<16;b++)Q.push(ax.map(k=>b>>>k&1?-1:1)); // small numerator /2
 const word=v=>v.join(','),lookup=new Map(Q.map((q,i)=>[word(q),i]));
 const id=ax.slice(),perm=[];
 function f(S,V){if(!S.length){perm.push(V);return;}for(const i of S)f(S.filter(x=>i!==x),[...V,i]);}f(ax,[]);
 const had=[[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]],inp=s=>Q.map(q=>lookup.get(word(s(q))));
 const aut=new Map(),G=[];
 for(const p of perm)for(let s=0;s<16;s++){
   const v=q=>p.map((j,k)=>q[j]*((s>>>k&1)?-1:1));
   const t=inp(v);if(t.some(x=>x===undefined))throw Error('aut');aut.set(t.join(','),t);
 }
 const h=inp(q=>had.map(r=>r.reduce((a,x,i)=>a+x*q[i],0)/2));
 const norm=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
 const seen=new Set(),out=[id.concat(Array.from({length:20},(_,i)=>4+i))];
 seen.add(out[0].join(','));
 const generators=[...aut.values(),h];
 for(let i=0;i<out.length;i++)for(const g of generators){const p=out[i].map(x=>g[x]),k=p.join(',');if(!seen.has(k)){seen.add(k);out.push(p);}}
 const PB=P.map(p=>P.map(q=>norm(p,q))),N=P.map(u=>Q.map(v=>norm(u,v)));
 const Pd=(u,v)=>u===v?4:(u^4)===v?-4:0,d=(u,v)=>Pd(u>>3,v>>3)+Pd(u&7,v&7);
 function orbits(k){
  const seen=new Set(),ret=[],S=[];
  function g(u){if(S.length===k){const m=word(S);if(seen.has(m))return;for(const p of out)seen.add(word(S.map(i=>p[i]).sort((a,b)=>a-b)));ret.push([...S]);return;}
   for(;u<24;u++){S.push(u);g(u);S.pop();}}
  g(0);return {ret,n:seen.size};
 }
 function graph(S,allow=0xFFFFFF){
  const X=[],sm=S.map((s,i)=>i);
  for(let i=0;i<24;i++)if(allow>>>i&1)for(let p=0;p<64;p++)if(S.every((s,j)=>{const a=d(p,sm[j]);return a===4||a===0&&N[i][s]===-4;}))X.push({i,p});
  const L=X.map(({i,p},j)=>X.flatMap(({i:f,p:q},k)=>j!==k&&i!==f&&p!==q&&PB[i][f]<=d(p,q)?[k]:[])),M=L.map(v=>new Set(v));
  return {X,L,M};
 }
 return {P,Q,out,PB,N,d,orbits,graph};
}
function solve101(O,limit=1e6,target=999,logAll=false,accept=()=>true){
 const {X,L,M}=O,n=X.length,A=[],all=[];
 let best=[],calls=0,stop=false;
 function serial(U){const B=[];for(const i of U){let S=B.find(Q=>Q.every(j=>!M[i].has(j)));if(S)S.push(i);else B.push([i]);}return B;}
 function go(U){
  if(++calls>limit){stop=true;return;}
  if(A.length>best.length)best=[...A];
  if(A.length>=target){if(logAll){all.push(A.map(i=>X[i]));if(all.length>=4096)stop=true;}
    else stop=true; return;
  }
  if(!U.length)return;
  const cs=serial(U);if(A.length+cs.length<(logAll?target:best.length+1))return;
  const Z=[],B=[];for(let i=0;i<cs.length;i++)for(const j of cs[i]){Z.push(j);B.push(i+1);}
  for(let i=Z.length-1;i>=0;i--){if(B[i]+A.length<(logAll?target:best.length+1))break;
   const j=Z[i];A.push(j);
   if(accept(A.map(i=>X[i])))go(Z.slice(0,i).filter(k=>M[j].has(k)));
   A.pop();if(stop)return;
  }
 }
 const root=Array.from({length:n},(_,i)=>i).sort((a,b)=>L[b].length-L[a].length);
 go(root);return {best:best.map(i=>X[i]),calls,complete:!stop,all};
}
const S=setup101();
assert.equal(S.out.length,1152);
for(const f of S.out){
 const inv=new Int8Array(24);f.forEach((x,i)=>inv[x]=i);
 const p=S.N.map(r=>S.N.findIndex(q=>q.every((z,k)=>z===r[inv[k]])));
 assert(!p.includes(-1));assert.equal(new Set(p).size,24);
 for(let i=0;i<24;i++)for(let j=0;j<24;j++)assert.equal(S.PB[i][j],S.PB[p[i]][p[j]]);
}
assert(S.Q.every((x,i)=>S.N.filter(r=>r[i]===-4).length===6));
function induce(g,V){const X=V.map(i=>g.X[i]),M=V.map(i=>new Set(V.flatMap((j,k)=>g.M[i].has(j)?[k]:[])));return{X,L:M.map(s=>[...s]),M};}
const res={};
const g=S.graph([0]);
const h=[-4,0,4].map(q=>{
 const i=S.N.findIndex(r=>r[0]===q), j=g.X.findIndex(z=>z.i===i&&z.p===1);
 assert(j>=0);const v=solve101(induce(g,g.L[j]),2000000,13);
 assert(v.complete&&v.best.length===11);
 return v.best.length+1;
});
res[1]=Math.max(...h);
for(let k=2;k<=4;k++){
 const or=S.orbits(k),H={},todo={2:[5,300],3:[15,2600],4:[59,17550]}[k];
 assert.deepEqual([or.ret.length,or.n],todo);
 let max=0;
 for(const small of or.ret){
  const repeats=small.map((x,i)=>[i,small.flatMap((y,j)=>j>i&&x===y?[j]:[])]);
  const accept=A=>k===2||repeats.every(([i,L])=>!L.length||A.filter(v=>S.d(v.p,i)===0).length<4||L.every(j=>A.filter(v=>S.d(v.p,j)===0).length<4));
  const v=solve101(S.graph(small),2000000,999,false,accept);
  assert(v.complete);max=Math.max(max,v.best.length);H[v.best.length]=(H[v.best.length]||0)+1;
 }
 assert.deepEqual(H,{2:{10:3,11:1,12:1},3:{5:1,6:14},4:{4:47,5:12}}[k]);res[k]=max;
}
assert.deepEqual(res,{1:12,2:12,3:6,4:5});
let tot=0;
for(let i=1;i<=4;i++)for(let j=i;j<=4;j++)for(let k=j;k<=4;k++)if(i+j+k>=9)tot=Math.max(tot,i+j+k+res[i]+res[j]+res[k]);
assert.equal(tot,32);
console.log(JSON.stringify({largeByNumberOfSmall:res,totalCapAlsoForAllReal:tot}));