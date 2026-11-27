// Node >=20. Run: node verify.mjs path/to/a40.8.8.4267H
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
const raw=readFileSync(process.argv[2]);
assert.equal(createHash('sha256').update(raw).digest('hex'),'0566ebf4da8b4997c9a782ea537ba3f35125d9c53201fa0cd94551a803ac47e2');
const smart=raw.toString('ascii'),need=(x,s)=>assert(x,s);
const pc=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return Math.imul((x+(x>>>4))&0xf0f0f0f,0x01010101)>>>24};
const mask=2**47-1;
const Q=[...new Set(Array.from({length:46},(_,i)=>(i+1)**2%47))].reduce((s,j)=>s+2**j,0);
const X=[], B=Array(48).fill(0n);
for(let j=0;j<47;j++){let w=((BigInt(Q)<<BigInt(j))|(BigInt(Q)>>BigInt(47-j)))&BigInt(mask);w|=1n<<47n;for(let t=0;t<48;t++)if((w>>BigInt(t))&1n){if(B[t])w^=B[t];else{B[t]=w;break}}}
need(B.filter(Boolean).length===24,'mother rank');
const G=B.slice(8).filter(Boolean).map(x=>x>>8n);
need(G.length===16,'shortened rank');
const L=G.map(w=>Number(w&0xffffffffn)),H=G.map(w=>Number(w>>32n)),ls=new Int32Array(65536),hs=new Int32Array(65536),wt=new Uint8Array(65536);
for(let i=1;i<65536;i++){const k=31-Math.clz32(i&-i);ls[i]=ls[i&(i-1)]^L[k];hs[i]=hs[i&(i-1)]^H[k];wt[i]=pc(ls[i])+pc(hs[i]);need(wt[i]>=12,'QR min12');}
const at=Array.from({length:40},(_,i)=>G.reduce((s,w,k)=>s+(((w>>BigInt(i))&1n)?1<<k:0),0));
const oct=smart.split(/\r?\n/).filter(x=>/^[a-f0-9]{10}$/i.test(x.trim())).map(x=>{const w=BigInt('0x'+x.trim()), ix=[];for(let i=0;i<40;i++)if((w>>BigInt(i))&1n)ix.push(i);return ix;});
need(oct.length===4267&&oct.every(x=>x.length===8),'Smart stock');
function fwht(A){for(let d=1;d<A.length;d*=2)for(let i=0;i<A.length;i+=2*d)for(let j=i;j<i+d;j++){const a=A[j],b=A[j+d];A[j]=a+b;A[j+d]=a-b;}return A;}

const order=oct.map(S=>S.join(','));
need(new Set(order).size===4267,'unique Smart octads');
// Constant-weight distance 8: no covered quintuple occurs twice.
const quin=new Set();
for(const S of oct)for(let a=0;a<8;a++)for(let b=a+1;b<8;b++)for(let c=b+1;c<8;c++)for(let d=c+1;d<8;d++)for(let e=d+1;e<8;e++){
 const k=((((40*S[a]+S[b])*40+S[c])*40+S[d])*40+S[e]);
 need(!quin.has(k),'duplicate Smart quintuple');quin.add(k);
}
const col=Array.from({length:40},()=>new Uint8Array(65536)), parity=[], spec={};
let maxAlpha=0;
for(const S of oct){
 let K=Array.from({length:16},(_,i)=>1<<i);
 for(const i of S){const p=K.find(x=>pc(x&at[i])&1);if(p!==undefined)K=K.filter(x=>x!==p).map(x=>(pc(x&at[i])&1)?x^p:x);}
 const n=1<<K.length, A=new Int32Array(n), V=new Int32Array(n);let degree=0;
 for(let i=1;i<n;i++){
  const u=V[i]=V[i&(i-1)]^K[31-Math.clz32(i&-i)];
  A[i]=(wt[u]===12);degree+=A[i];
  if(wt[u]>12)for(const j of S){need(col[j][u]<255,'agreement counter overflow');++col[j][u];}
 }
 fwht(A);const low=A.reduce((s,x)=>Math.min(s,x),0);
 const alpha=degree?Math.floor(n*(-low)/(degree-low)):n;
 maxAlpha=Math.max(maxAlpha,alpha);spec[alpha]=(spec[alpha]||0)+1;
 parity.push(S.reduce((a,j)=>a^at[j],0));
}
assert.equal(maxAlpha,53);
const mins=[],covers=[],excess=[],degrees=[],cut=[];
for(let i=0;i<40;i++){
 const f=new Int32Array(65536);let n=0;
 for(let j=0;j<oct.length;j++)if(oct[j].includes(i)){++f[parity[j]];++n;}
 assert.equal(f[0],0);
 const rep=n-f.reduce((a,j)=>a+Number(j>0),0);
 need(rep<=11&&n>=852&&n<=855,'nonzero distinct incident parity characters');
 fwht(f);let t=n;
 for(let v=1;v<65536;v++)if(wt[v]!==12||(pc(v&at[i])&1)){
  need((n-f[v])%2===0,'integer Walsh disagreement');
  t=Math.min(t,(n-f[v])/2);
 }
 const r=col[i].reduce((a,b)=>Math.max(a,b),0);
 mins.push(t);covers.push(r);excess.push(rep);degrees.push(n);
 need(t>=340&&r<=33,'allowed pair pressure');
 for(let m=2;m<=85;m++)need(t*m>=16*Math.floor(m*m/4),'mixed-sign arbitrary bank payment');
 // Minimum possible total A of chosen-parity parent multiplicities, and collision lower bound.
 function charge(m){
  let lo=0,hi=Math.floor(m*n/2);
  const F=A=>2*m*n*A-2*A*A-t*m*(m-1)*n;
  need(F(hi)>=0,'nonempty quadratic relaxation');
  while(lo<hi){const k=(lo+hi)>>1;if(F(k)>=0)hi=k;else lo=k+1;}
  const den=lo+r*m*(m-1), num=lo*lo;
  need(Number.isSafeInteger(num)&&Number.isSafeInteger(8*(m-1)*den),'integer collision bound');
  return num>=8*(m-1)*den;
 }
 let m=2;while(m<=512&&charge(m))++m;
 cut.push(m-1);need(m-1>=212,'fixed-parity arbitrary bank payment');
 // Rank k>=8: total active S count exceeds 8 times 64, after inactive-character budget.
 for(let k=8;k<=16;k++)need(n-((1<<(16-k))-1)-rep>=512,'large affine bank payment');
}
assert.equal(Math.min(...mins),340);
assert.equal(Math.max(...covers),33);
assert.equal(Math.max(...excess),11);
assert.equal(Math.min(...cut),212);
console.log('AUDIT PASS',JSON.stringify({maxAlpha,parityDisagreementMin:mins,fullRestrictionCollisionMax:covers,repeatExcess:excess,fixedParityBankCaps:cut,degrees,uniformMixedSignArbitraryCap:85}));
