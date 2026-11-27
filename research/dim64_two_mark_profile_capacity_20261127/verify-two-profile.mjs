#!/usr/bin/env node
// Node 18+; only exact signed 32-bit arithmetic, no data files.
import assert from "node:assert/strict";
const z=(()=>{const M=1<<22;const E=new Int32Array(M);const forbidden=[], permitted=[], all=[];const ij=Array.from({length:6},()=>new Uint32Array(6));let t=7;for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)ij[i][j]=1<<t++;
const qs=new Set();
for(let a=1;a<64;a++)for(let b=a+1;b<64;b++){let c=0;for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)if((((a>>i)&1)*((b>>j)&1))^(((a>>j)&1)*((b>>i)&1)))c^=ij[i][j];if(qs.has(c))continue;qs.add(c);
for(let A=0;A<2;A++)for(let B=0;B<2;B++){let L=((a&b)^(A?b:0)^(B?a:0))<<1;let x=c^L^(A&B);all.push(x);if((A&B)&&(((a&1)^A)&((b&1)^B)))permitted.push(x);else{forbidden.push(x);E[x]=1;}}}
for(let l=1;l<M;l*=2)for(let j=0;j<M;j+=2*l)for(let k=j;k<j+l;k++){let u=E[k],v=E[k+l];E[k]=u+v;E[k+l]=u-v;}
let hist={}, min=0,minchars=[];for(let j=0;j<M;j++){let a=E[j];hist[a]=(hist[a]||0)+1;if(a<min){min=a;minchars=[];}if(a===min&&minchars.length<10)minchars.push(j);}
return {forbidden,permitted,all,hist,min,minchars};})();
assert.equal(z.all.length,2604);
assert.equal(z.permitted.length,155);
assert.equal(z.forbidden.length,2449);
assert.equal(new Set(z.all).size,2604);
// Independently evaluate each generator at both marks, 0 and e_1.
const at=(x,v)=>(x&1)^((x>>1)&(v&1));
assert(z.permitted.every(x=>at(x,0)&&at(x,1)));
assert(z.forbidden.every(x=>!(at(x,0)&&at(x,1))));
assert.equal(z.min,-79);
assert.equal(z.hist[-79],13888);
assert.equal(Object.values(z.hist).reduce((a,b)=>a+b),1<<22);
assert.equal((1<<22)*79/(2449+79),131072);
console.log(JSON.stringify({status:"PASS",full:z.all.length,permitted:z.permitted.length,degree:z.forbidden.length,min:z.min,multiplicity:z.hist[-79],alpha:131072}));