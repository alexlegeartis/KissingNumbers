/* Portable node replay: exact finite data and inequalities.
   Source of the 19D E code, 19D body, and 21D G code: Kevin Li,
   https://hexagonmath.org/pdf/2610.00248v1, Sections 3-4.
   Usage: node verify.mjs. No packages or external files. */
import assert from 'node:assert/strict';

const M19=(()=>{
const TS=[0,3,6,8],T=TS.reduce((x,j)=>x|(1<<j),0);
const G0=[[1,8,9,12,16,17,18,19],[2,10,11,14,15,17,18,20],[3,7,9,13,15,16,17,19],[4,7,8,10,12,15,16,19],[5,10,12,13,15,16,17,18],[6,7,8,9,10,13,16,18],[1,4,7,9],[1,5,6,18,20],[1,3,12,15,20],[1,10,13,19,20],[1,3,5,6,7,13,14,15,18,20],[2,4,6,7,8,13,14,16,17,18]].map(v=>v.reduce((x,j)=>x|1<<(j-1),0));
const pc=x=>x.toString(2).replaceAll('0','').length;
let E=[0],dual=[0],piv=Array(20).fill(0);
for(let g of G0){E.push(...E.map(x=>x^g));for(let j=19;j>=0;j--)if(g>>>j&1){if(piv[j])g^=piv[j];else{piv[j]=g;break;}}}
for(let j=0;j<20;j++)if(!piv[j]){let x=1<<j;for(let k=0;k<20;k++)if(pc(x&piv[k])%2)x|=1<<k;dual.push(...dual.map(a=>a^x));}
if(E.length!==4096||dual.length!==256||dual.some(w=>G0.some(g=>pc(g&w)%2)))throw Error("binary sources");
const O=dual.filter(w=>pc(w)===8),Z=E.filter(w=>(w&T)===0);
const B=O.filter(w=>pc(w&T)===2),F={};let hh={};
for(let w of B){
 const f=new Set(Z.map(z=>z&w)),r=Math.log2(f.size),key=w&T;F[key]??=[];
 F[key].push([w&~T,r]);hh[r]=(hh[r]||0)+1;
}
return {T,TS,Z,F,hh,O,B,pc};
})();
const {Z,F,hh,O,B,TS,pc}=M19;
assert.equal(Z.length,256);
assert.equal(new Set(Z).size,256);
assert.equal(O.length,130);
assert.equal(B.length,96);
assert.deepEqual(hh,{'5':96});
assert.equal(Object.keys(F).length,6);
assert(Object.values(F).every(rows=>rows.length===16 && new Set(rows.map(row=>row[0])).size===16));
for(const o of B) {
  const freq=new Map();
  for(const z of Z) {
    const t=z&o;
    assert.equal(pc(t)%2,0);
    freq.set(t,(freq.get(t)||0)+1);
  }
  assert.equal(freq.size,32);
  assert([...freq.values()].every(n=>n===8));
}
const R=[[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]];
const vv=[];
for(let i=0;i<4;i++) for(let j=0;j<4;j++) if(i!==j)
  vv.push(R.map(r=>(r[i]-r[j])/2).join(','));
const want=[];
for(let i=0;i<3;i++) for(let j=i+1;j<3;j++)
  for(const x of [-1,1]) for(const y of [-1,1]) {
    const v=[0,0,0];v[i]=x;v[j]=y;want.push(v.join(','));
  }
assert.deepEqual(vv.sort(),want.sort());
assert.deepEqual(TS,[0,3,6,8]);
assert.equal(12*16*32,6144);
assert.equal(480+12+3840+6144+192+1792,12460);
// A=117/20, B=82/20; the three-direction cut is beta=47/82.
assert(117n**2n>2n*82n**2n); // A > sqrt(2) B
assert(2n*47n**2n<82n**2n);
assert(11n*47n**2n<4n*82n**2n);
const kappaN=760n,kappaD=761n;
assert(2n*kappaN**2n>kappaD**2n);
assert(500n*(kappaD-kappaN)*kappaD<2n*kappaN**2n-kappaD**2n);
assert(500n*47n**2n>2n*82n**2n); // c^2 > 1/500
assert.equal(3*16,2*8*3);
console.log('19D: 256 zero-tetrad heads; 96 rank-five supports; 12 D3 orientations; charge 48/24 (two for one).');

const M21=(()=>{
const pc=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0x0f0f0f0f)*0x01010101)>>>24;};
const g=[0,1,5,6,7,9,11].reduce((a,x)=>a|1<<x,0);
const gen=Array.from({length:12},(_,i)=>(g<<i)|((pc(g)&1)<<23));
let G=[0];for(let v of gen)G.push(...G.map(x=>v^x));
const D=G.map(x=>x&((1<<21)-1)).sort((a,b)=>a-b);
const O=G.filter(x=>x<1<<21&&pc(x)===8);
const L=D.filter(x=>pc(x)===5);
if(new Set(D).size!==4096||O.length!==210||L.length!==21)throw Error("Golay counts");
const os=O.map(x=>Array.from({length:21},(_,i)=>i).filter(i=>x>>>i&1));
let cd=(x,i)=>(pc(x)&1)^(x>>>i&1);
const H=L[0],D0=Array.from({length:21},(_,i)=>i).filter(i=>H>>>i&1),A0=Array.from({length:21},(_,i)=>i).filter(i=>!(H>>>i&1));
const pmask=(p)=>D0.reduce((r,i,j)=>r|(cd(p,i)<<j),0);
const Cs=[...new Set(D.map(pmask))];
if(Cs.length!==16||Cs.some(p=>pc(p)%2))throw Error("pivots");
const pairs=[];
for(let i=0;i<21;i++)for(let j=i+1;j<21;j++){
 const m=L.find(x=>(x>>>i&1)&&(x>>>j&1));
 for(let k=0;k<D0.length;k++)if(D0[k]!==i&&D0[k]!==j&&(m>>>D0[k]&1)){
  let I=D0[k],a=L.map(P=>+(P>>>i&1)+ +(P>>>j&1)),f=L.map(P=>!(P>>>I&1)?1:0);
  if(a.some((x,p)=>!f[p]&&x===1))throw Error("bad cancellation");
  pairs.push({i,j,k,m,g:L.map((P,p)=>1-(f[p]?-1:1)*(2-2*a[p]))});
 }
}
const near=D.filter(x=>pc(x)<=7&&x>0);
return {D,O,L,os,pc,gen,pmask,cd,seen:Cs,A0,D0,pairs,near};
})();
const {D,L,os,pairs,pmask,gen}=M21;
assert.equal(D.length,4096);
assert(gen.every((g,i)=>gen.every((h,j)=>M21.pc(g&h)%2===0)));
assert.equal(L.length,21);
assert(D.every(w=>M21.O.every(o=>(M21.pc(w&o)&1)===0)));
for(let i=0;i<21;i++) {
  assert.equal(L.filter(P=>P>>>i&1).length,5);
  assert.equal(os.filter(O=>O.includes(i)).length,80);
  for(let j=i+1;j<21;j++){
    assert.equal(L.filter(P=>(P>>>i&1)&&(P>>>j&1)).length,1);
    assert.equal(os.filter(O=>O.includes(i)&&O.includes(j)).length,28);
  }
}
assert.equal(new Set(D.map(pmask)).size,16);
assert(D.every(w=>M21.pc(pmask(w))%2===0));
assert(pairs.length>0);
for(const {g} of pairs) {
  assert.equal(g.filter(x=>x===-1).length,4);
  assert.equal(g.filter(x=>x===1).length,8);
  assert.equal(g.filter(x=>x===3).length,9);
  assert.equal(g.reduce((x,y)=>x+y,0),31);
}
for(const k of [9,10]) {
  assert(52*2*k>4*210); // positive RMS > 2 after maximizing over phase
  assert(21*50**2>62*k*9**2); // q > 9/50
  assert(60+2*(25-4*k)>0); // polynomial derivative > 0 on [0,1]
  assert(-6*50**2+60*9*50+(25-4*k)*9**2>0);
}
assert.equal(128*210+840+3840,31560);
console.log('21D: 16 even phases; (r,lambda)=(80,28); every legal pair has g census 4,8,9; uniform k=9,10 impossible.');

const H=(()=>{
function mul(a,b){let z=0;for(;b;b>>=1){if(b&1)z^=a;a<<=1;if(a&4)a^=7;}return z;}
const conj=a=>mul(a,a),n=a=>a===0?0:1;
const dot=x=>x.toString(2).replaceAll("0","").length&1;
const pack=(a,b,c)=>a|(b<<1)|(c<<2);
function B(a,b){
 const [x,y]=[a&3,a>>>2], [u,v]=[b&3,b>>>2];
 return pack(mul(x,conj(u))^mul(u,conj(x)),mul(y,conj(v))^mul(v,conj(y)),mul(x,conj(v))^mul(u,conj(y)));
}
const C=Array.from({length:16},(_,z)=>Array.from({length:16},(_,u)=>B(z,u)));
const F=C.map((A,z)=>{
 const S=Array.from({length:16},(_,v)=>v).filter(v=>A.every(a=>!dot(a&v)));
 if(z&&!(S.length===2&&S[0]===0))throw Error("annihilator");
 return S[1];
});
function fitting(Z){
 const R=new Map();let bad=false,cert=null;
 for(let a of Z)for(let b of Z)if(a<b){
  let u=BigInt(F[a^b])<<BigInt(4*a)|BigInt(F[a^b])<<BigInt(4*b),c=1;
  while(u){let k=u.toString(2).length-1;if(R.has(k)){let [v,t]=R.get(k);u^=v;c^=t;}else{R.set(k,[u,c]);break;}}
  if(!u&&c){bad=true;break;}
 }
 if(bad)return null;
 let s=0n;
 for(let k=0;k<64;k++)if(R.has(k)){let [u,t]=R.get(k);
  let ct=0;for(let a=u&s;a;a&=a-1n)ct^=1;
  if(ct!==t)s|=1n<<BigInt(k);
 }
 return {r:R.size,H:Z.map(z=>Number(s>>BigInt(4*z)&15n))};
}
return {C,F,fitting};
})();
const fs=[...new Set(H.F.slice(1))].sort((a,b)=>a-b);
assert.deepEqual(fs,[1,2,7,11,15]);
for(let j=0;j<fs.length;j++){
  const basis=fs.filter((_,i)=>i!==j);
  let sums=[0];for(const a of basis)sums.push(...sums.map(t=>t^a));
  assert.equal(new Set(sums).size,16);
}
assert.deepEqual(fs.map(f=>H.F.slice(1).filter(x=>x===f).length),[3,3,3,3,3]);
for(let x=1;x<16;x++)for(let y=1;y<16;y++)if(H.F[x]===H.F[y]&&x!==y)
  assert.equal(H.F[x^y],H.F[x]);
console.log('21D Hermitian split: five rank-three polar directions; any four annihilators independent.');

