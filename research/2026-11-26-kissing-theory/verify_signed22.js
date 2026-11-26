'use strict';
const assert = require('node:assert/strict');
const pc=n=>{let q=0;for(;n;n&=n-1)q++;return q;};
const g=[0,1,5,6,7,9,11].reduce((a,b)=>a+(1<<b),0);
let G=[0];for(let i=0;i<12;i++){const p=g<<i,r=p|((pc(p)&1)<<23);G.push(...G.map(x=>x^r));}
const mask=(1<<21)-1, contains=new Set(G);
assert.equal(contains.size,4096);
const O=G.filter(x=>(x>>>21)===0&&pc(x)===8),L=G.filter(x=>(x>>>21)===7&&pc(x)===8).map(x=>x&mask);
assert.equal(O.length,210);assert.equal(L.length,21);
assert(O.every((a,i)=>O.slice(0,i).every(b=>pc(a&b)<=4)));
assert(L.every((a,i)=>L.slice(0,i).every(b=>pc(a&b)===1)));
assert(O.every(a=>L.every(b=>pc(a&b)<=4)));
const rank=v=>{let B=Array(24).fill(0);for(let x of v)for(let j=23;j>=0;j--)if(x>>j&1){if(B[j])x^=B[j];else{B[j]=x;break;}}return B.filter(Boolean).length;};
for(let i=0;i<21;i++){
 const S=O.filter(x=>x>>i&1);assert.equal(S.length,80);assert.equal(rank(S),9);
}
const triples=[[0x2b2230,0x532026,0x820336],[0x2ad060,0x564460,0x9e8048],[0x216f00,0x542e04,0x850d14],[0x2b2230,0x42a294,0x8120f4]];
for(let j=0;j<21;j++)assert(triples.some(S=>S.every(x=>(x>>j&1)===0)));
for(const S of triples){
 assert(S.every(x=>contains.has(x)));assert.equal(rank(S.map(x=>x>>>21)),3);
 let H=[0];for(const x of S)H.push(...H.map(y=>x^y));
 assert(H.slice(1).every(x=>pc(x)===8&&pc(x&mask)>=5&&pc(x&mask)<=7));
}
const row=(t,v)=>[1,t,t*(t-1)/2,v,(v*v-t)/2,(v*v*v-3*v*t+2*v)/6,t*(v*v-t)/2,t*v];
const f4=(t,v)=>(v-1)*(v+2)*(2*v+t+1);
for(let t=0;t<=5;t++)for(let v=-t;v<=t;v+=2){
 const [n,s,k,w,a,b,c,d]=row(t,v);
 assert.equal(f4(t,v),12*b+7*(d-w)+2*c+6*a+t*t+t-2);
 if(v<=1)assert(f4(t,v)<=0);
}
assert.equal(3*f4(5,5),4*336);
const T=[[1,1],[1,-1],[2,-2],[3,1],[4,-2],[5,1],[5,-3]],mul=[160,160,0,0,0,10,5],goal=[336,400,160,0,0,0,0,0];
for(let j=0;j<8;j++)assert.equal(T.reduce((a,v,i)=>a+mul[i]*row(...v)[j],row(5,5)[j]),goal[j]);
function det(M){return M.length===1?M[0][0]:M[0].reduce((z,v,j)=>z+v*((j&1)?-1:1)*det(M.slice(1).map(r=>r.filter((_,i)=>i!==j))),0);}
const A=Array.from({length:7},(_,j)=>T.map(v=>row(...v)[j]));assert.equal(Math.abs(det(A)),672);
console.log('PASS: Golay support counts, 21 rank tests, four explicit incompatibility subspaces, polynomial and equality row system.');
