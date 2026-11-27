// Node 20+: node verify.mjs witness.npz
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {inflateRawSync} from 'node:zlib';
const decode=(BASE,inflate)=>{
const raw=Uint8Array.from(BASE);
const view=new DataView(raw.buffer,raw.byteOffset,raw.byteLength);
const u16=(p)=>view.getUint16(p,true),u32=(p)=>view.getUint32(p,true);
let e=raw.length-22;while(e>=0&&u32(e)!==0x06054b50)e--;if(e<0)throw Error("EOCD");
let n=u16(e+10),p=u32(e+16), out={};
for(let i=0;i<n;i++){
 if(u32(p)!==0x02014b50)throw Error("central");
 const z=u32(p+20),l=u32(p+24),k=u16(p+28),q=u32(p+42);
 const name=String.fromCharCode(...raw.slice(p+46,p+46+k));
 const off=q+30+u16(q+26)+u16(q+28), chunk=raw.slice(off,off+z);
 const data=u16(p+10)===8?inflate(chunk,l):chunk;
 if(data.length!==l||data[0]!==147)throw Error("npy");
 const v=data[6],d=new DataView(data.buffer,data.byteOffset,data.byteLength);
 const start=v===1?10:12;
 const end=start+(v===1?d.getUint16(8,true):d.getUint32(8,true));
 const head=String.fromCharCode(...data.slice(start,end));
 const dtype=head.match(/'descr': '([^']+)/)?.[1],shape=head.match(/'shape': \(([^)]*)/)?.[1]?.split(",").map(s=>s.trim()).filter(Boolean).map(Number);
 let C={"<i2":Int16Array,"<i4":Int32Array,"<i8":BigInt64Array,"|i1":Int8Array,"<u2":Uint16Array,"<u4":Uint32Array,"<u8":BigUint64Array,"|u1":Uint8Array}[dtype];
 if(!C)throw Error("missing dtype "+dtype);
 const a=new C(data.buffer,data.byteOffset+end,(data.length-end)/C.BYTES_PER_ELEMENT);
 if(a.length!==shape.reduce((s,k)=>s*k,1))throw Error("shape "+[name,head]);
 out[name.replace(".npy","")]={shape,dtype,a};
 p+=46+k+u16(p+30)+u16(p+32);
}
return out;
};

const raw=readFileSync(process.argv[2]);assert.equal(createHash('sha256').update(raw).digest('hex'),'90d130ba1bb96172fae7b7ae10d1169747b76f2e720946bc2350c99a08d427c9');
const X=decode(raw,inflateRawSync);
const TRADES=[[1074431249,285933841,1074955537,336],[4580344841,4580183809,4580245519,612],[4295055681,4312861953,4563472593,744],[4295055681,4564520193,4563472593,768],[4429259523,4514185473,4563470097,584],[4429259523,4546691329,4563470097,584],[4568654401,4581257217,4564453137,516],[2170621961,4578152449,285282571,664],[2170621961,4578152449,2969637121,600]];
const maskbits=x=>Array.from({length:33},(_,i)=>i).filter(i=>Math.floor(x/2**i)&1);
const B=Array.from(X.supports.a,Number), oldMasks=new Set(B), par=Array.from(X.parity.a,Number),v=Array.from({length:4320},(_,i)=>Array.from(X.sparse.a.slice(33*i,33*(i+1))));
assert.equal(B.length,1837);assert.equal(oldMasks.size,1837);assert.deepEqual(Array.from(X.dense_positions.a),Array.from({length:32},(_,i)=>i));
assert.equal(X.dense_generator.a.length,17);
for(const R of v){
 assert.equal(R.reduce((s,x)=>s+x*x,0),32);assert(R.every(x=>x%2===0));
}
const ix=B.map(maskbits),N=33,C=Array.from({length:34},()=>[1,0,0,0,0,0]);for(let i=1;i<34;i++)for(let j=1;j<=5;j++)C[i][j]=C[i-1][j]+C[i-1][j-1];
const f=new Uint16Array(C[N][5]);
for(const [i,A] of ix.entries()){
 assert.equal(A.length,8);
 for(let a=0;a<4;a++)for(let b=a+1;b<5;b++)for(let c=b+1;c<6;c++)for(let d=c+1;d<7;d++)for(let e=d+1;e<8;e++){
  const k=A[a]+C[A[b]][2]+C[A[c]][3]+C[A[d]][4]+C[A[e]][5];assert.equal(f[k],0);f[k]=i+1;
 }
}
const h=new Int32Array(34),H=[],p=[];
function go(j){
 ++h[p.length];if(p.length>=9)H.push(p.slice());
 for(;j<N;j++){const x=C[j][5],k=p.length;let bad=false;
 outer:for(let a=0;a<k;a++)for(let b=a+1;b<k;b++)for(let c=b+1;c<k;c++)for(let d=c+1;d<k;d++)
  if(f[p[a]+C[p[b]][2]+C[p[c]][3]+C[p[d]][4]+x]){bad=true;break outer;}
 if(!bad){p.push(j);go(j+1);p.pop();}
 }
}
go(0);
assert.deepEqual(Array.from(h.slice(5,11)),[134464,79706,3044,197,18,0]);
// Bounds on sqrt(32/7), sqrt(24/7) are strict integer-square comparisons.
for(const [a,b,q] of [[21380,21381,32],[18516,18517,24]])assert(7*a*a<q*10**8&&7*b*b>q*10**8);
assert(21381+18517<40000&&21381+3*18517<80000);
let total=0,min=100,max=0,keys=new Set(),surv=0;
for(const A of H){
 assert.equal(A.length,9);assert.equal(A[8],32);
 const r=v.flatMap((z,id)=>{
  const t=A.map(i=>z[i]);return t.reduce((s,x)=>s+Math.abs(x),0)>=8?[[id,t]]:[];
 });
 for(let s=0;s<512;s++){
  const S=A.map((_,i)=>s>>i&1?-1:1);
  let k=0;
  for(const [id,t] of r){
   const a=t[8]*S[8],b=t.slice(0,8).reduce((sum,z,i)=>sum+z*S[i],0);
   const l=a*(a>=0?21380:21381)+b*(b>=0?18516:18517);
   const u=a*(a>=0?21381:21380)+b*(b>=0?18517:18516);
   const isotropic=(a+b>8);
   assert(isotropic?l>160000:u<=160000);
   if(isotropic){assert(v[id].filter(x=>Math.abs(x)===2).length===8);keys.add(id);++k;}
  }
  min=Math.min(min,k);max=Math.max(max,k);++total;if(!k)++surv;
 }
}
assert.deepEqual([total,min,max,keys.size,surv],[9216,8,32,1864,0]);
const overlap=(A,B)=>A.filter(i=>B.includes(i)).length;
for(const [old,a,b,expected] of TRADES){
 assert(oldMasks.has(old)&&!oldMasks.has(a)&&!oldMasks.has(b));
 const A=maskbits(a),B=maskbits(b);assert(A.length===8&&B.length===8&&overlap(A,B)<=4);
 for(const Z of [A,B])assert.deepEqual(ix.flatMap((Q,i)=>overlap(Q,Z)>4?[i]:[]),[Array.from(X.supports.a,Number).indexOf(old)]);
 function price(c){return v.flatMap((R,i)=>{const signs=R.map(Math.abs);const bit=c[1];const t=c[0].map(x=>R[x]);const z=t.reduce((s,x)=>s+2*Math.abs(x),0)-4*Math.min(...t.map(Math.abs))*((t.filter(x=>x<0).length%2)!==bit);return z>16?[i]:[];});}
 for(let p=0;p<=1;p++)for(let q=0;q<=1;q++)assert.equal(new Set([...price([A,p]),...price([B,q])]).size,expected);
}
console.log('AUDIT PASS',JSON.stringify({fullOctads:ix.length,oldExplicit:v.length,uncoveredNine:H.length,donors:total,minBlockers:min,maxBlockers:max,uniqueOldBlockers:keys.size,trades:TRADES}));
