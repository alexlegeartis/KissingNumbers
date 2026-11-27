// Node >=20; canonical QR prefix compatibility replay, no data input.
import assert from 'node:assert/strict';
const need=(x,s)=>assert(x,s);
const pc=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return Math.imul((x+(x>>>4))&0xf0f0f0f,0x01010101)>>>24};
const mask=2**47-1;
const Q=[...new Set(Array.from({length:46},(_,i)=>(i+1)**2%47))].reduce((s,j)=>s+2**j,0);
const B=Array(48).fill(0n);
for(let j=0;j<47;j++){let w=((BigInt(Q)<<BigInt(j))|(BigInt(Q)>>BigInt(47-j)))&BigInt(mask);w|=1n<<47n;for(let t=0;t<48;t++)if((w>>BigInt(t))&1n){if(B[t])w^=B[t];else{B[t]=w;break}}}
need(B.filter(Boolean).length===24,'mother rank');
const G=B.filter(Boolean);
const L=G.map(x=>Number(x&0xffffffffn)),H=G.map(x=>Number(x>>32n));
const min=new Uint8Array(256).fill(40),pop=new Int32Array(256);
let x=0,y=0; ++pop[0]; min[0]=0;
for(let i=1;i<1<<24;i++){
 const k=31-Math.clz32(i&-i); x^=L[k];y^=H[k];
 const a=x&255,d=pc(x>>>8)+pc(y);
 need(pc(x)+pc(y)>=12,'mother min12');
 ++pop[a];if(d<min[a])min[a]=d;
}
need(pop.every(n=>n===65536),'full rank-eight mark map');
const good=a=>a!==0&&min[a]>=10;
const X=Array.from({length:256},(_,a)=>a).filter(good);
assert.deepEqual(X.filter(a=>pc(a)>2),[95,125,189,190,215,219,231,250]);
for(let a=1;a<256;a++)if(pc(a)<=2)need(good(a),'ordinary radius-one safe differences');
const centers=[0,1,2,4,8,16,32,64,128];
const sizes=new Int32Array(10);let max=0;
function go(A,M){
 const k=A.length;max=Math.max(max,k); need(k<=9,'maximum nine whole sheets');
 if(k>=5){
  ++sizes[k];
  const C=centers.filter(c=>A.every(a=>pc(a^c)<=1));
  assert.equal(C.length,1);
  const U=X.filter(x=>A.every(a=>good(x^a)));
  assert.equal(U.length,9-k);
  need(U.every(a=>pc(a^C[0])<=1),'all partial-compatible marks lie in the same ball');
 }
 while(k+M.length>=5&&M.length){
  const x=M.shift();go([...A,x],M.filter(y=>good(x^y)));
 }
}
go([0],X.slice());
assert.equal(max,9);
assert.deepEqual(Array.from(sizes.slice(5)),[630,504,252,72,9]);
console.log('AUDIT PASS',JSON.stringify({max,exoticDifferences:X.filter(a=>pc(a)>2),normalizedCliqueCounts:Array.from(sizes.slice(5))}));
