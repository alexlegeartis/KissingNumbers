'use strict';
// Exact two-heavy signed carrier census in the literal 370556 scaffold.
const assert=require('node:assert/strict'),fs=require('node:fs'),crypto=require('node:crypto');
const pc=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return Math.imul((x+(x>>>4))&0xf0f0f0f,0x01010101)>>>24};
const alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
function base64(s){let x=0,k=0,out=[];for(const c of s){const v=alphabet.indexOf(c);if(v<0)continue;x=x<<6|v;k+=6;if(k>=8){k-=8;out.push(x>>>k&255);}}return Uint8Array.from(out);}
function inf(by,n){let at=0,k=0,b=0,w=0;const O=new Uint8Array(n);
const rd=(j)=>{while(k<j){if(at>=by.length)throw Error('EOF');b|=by[at++]<<k;k+=8;}const x=b&((1<<j)-1);b>>>=j;k-=j;return x};
const R=(L)=>{const N=new Int32Array(16),F=new Int32Array(16),T=Array.from({length:16},()=>new Map());for(const x of L)if(x)N[x]++;for(let i=1;i<16;i++)F[i]=(F[i-1]+N[i-1])<<1;for(let i=0;i<L.length;i++)if(L[i])T[L[i]].set(F[L[i]]++,i);return T;};
const sym=T=>{let x=0;for(let i=1;i<16;i++){x=x<<1|rd(1);if(T[i].has(x))return T[i].get(x);}throw Error('code');};
const LB=[3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258],LE=LB.map((x,i)=>i<8||i===28?0:(i>>2)-1);
const DB=[1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577],DE=DB.map((x,i)=>i<4?0:(i>>1)-1);
let done=0;while(!done){done=rd(1);const flag=rd(2);let L,D;
 if(flag===0){rd(k&7);const m=rd(16);if((m^rd(16))!==65535)throw Error('stored');for(let i=0;i<m;i++)O[w++]=rd(8);continue;}
 if(flag===1){L=R(Array.from({length:288},(_,i)=>i<144?8:i<256?9:i<280?7:8));D=R(Array(32).fill(5));}
 if(flag===2){const m=rd(5)+257,n=rd(5)+1,g=rd(4)+4,ord=[16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15],u=Array(19).fill(0);for(let i=0;i<g;i++)u[ord[i]]=rd(3);
  const T=R(u),S=[];while(S.length<m+n){const x=sym(T);if(x<16)S.push(x);else if(x===16)S.push(...Array(rd(2)+3).fill(S.at(-1)));else if(x===17)S.push(...Array(rd(3)+3).fill(0));else S.push(...Array(rd(7)+11).fill(0));}
  L=R(S.slice(0,m));D=R(S.slice(m));}
 if(!L)throw Error('invalid block');
 for(;;){const x=sym(L);if(x<256)O[w++]=x;else if(x===256)break;else {const z=LB[x-257]+rd(LE[x-257]),q=sym(D),s=DB[q]+rd(DE[q]);if(w<s)throw Error('bad distance');for(let y=0;y<z;y++){O[w]=O[w-s];w++;}}}
}
if(w!==n)throw Error('length');return O;}
function decode(arr){
 const by=Uint8Array.from(arr),V=new DataView(by.buffer);const rd=(j)=>V.getUint32(j,true),rh=(j)=>V.getUint16(j,true);
 let e=by.length-22;while(e>=0&&rd(e)!==0x06054b50)e--;if(e<0)throw Error('missing dir');
 const n=rh(e+10);let j=rd(e+16),Z={};for(let i=0;i<n;i++){
 if(rd(j)!==0x02014b50)throw Error('central dir');
 const sz=rd(j+20),out=rd(j+24),m=rh(j+28),x=rd(j+42),name=String.fromCharCode(...by.slice(j+46,j+46+m)),off=x+30+rh(x+26)+rh(x+28);
 const T=by.slice(off,off+sz),U=rh(j+10)===8?inf(T,out):T;
 if(U[0]!==147||U.length!==out)throw Error('NPY');
 const v=new DataView(U.buffer,U.byteOffset,U.byteLength),q=U[6]===1?10:12,end=q+(q===10?v.getUint16(8,true):v.getUint32(8,true));const head=String.fromCharCode(...U.slice(q,end)),D=head.match(/'descr': '([^']+)/)[1],shape=head.match(/'shape': \(([^)]*)/)[1].split(',').map(x=>x.trim()).filter(Boolean).map(Number);
 const cls={'<i2':Int16Array,'<i4':Int32Array,'<i8':BigInt64Array,'|i1':Int8Array,'<u2':Uint16Array,'<u4':Uint32Array,'<u8':BigUint64Array,'|u1':Uint8Array}[D];
 if(!cls)throw Error(D);
 const data=Array.from(new cls(U.buffer,U.byteOffset+end,(U.length-end)/cls.BYTES_PER_ELEMENT),Number);
 if(data.length!==shape.reduce((a,b)=>a*b,1))throw Error('shape');
 Z[name.replace('.npy','')]={shape,data};j+=46+m+rh(j+30)+rh(j+32);
 }
 return Z;
}


function scope(D){ const j={full_even_supports:D.supports.data};
const B=j.full_even_supports.map(v=>[v|0,(v>=2**32)?1:0]);const M=B.map(([b,h])=>BigInt(b>>>0)|(BigInt(h)<<32n));
const p=Array.from({length:33},(_,i)=>1n<<BigInt(i));const A=Array.from({length:33},(_,i)=>i);const cov=new Set();for(const X of M){const R=A.filter(i=>X&p[i]);for(let a=0;a<8;a++)for(let b=a+1;b<8;b++)for(let c=b+1;c<8;c++)for(let d=c+1;d<8;d++)cov.add(p[R[a]]|p[R[b]]|p[R[c]]|p[R[d]]);}
function comb(L,k,f){const P=[];function rec(i){if(P.length===k){f(P);return;}for(;i<=L.length-k+P.length;i++){P.push(L[i]);rec(i+1);P.pop()}}rec(0)}

let coarse=0,one=0,two=0;
const E=[];
for(const i of A)for(const j of A)if(i<j){
 const L=A.filter(k=>k!==i&&k!==j);
 const U=[];for(const a of L)for(const b of L)if(a<b&&!cov.has(p[i]|p[j]|p[a]|p[b]))U.push([a,b]);
 for(const [a,b] of U)for(const c of L)if(c>b&&!cov.has(p[i]|p[j]|p[a]|p[c])&&!cov.has(p[i]|p[j]|p[b]|p[c])){
  coarse++;
  const T=p[a]|p[b]|p[c];
  if(M.some(m=>(m&(p[i]|p[j]))&&(m&T)===T))continue;
  one++;
  const H=p[i]|p[j]; let bad=H|T;
  for(const m of M)if((m&H)===H&&(m&T))bad|=m;
  const V=A.filter(k=>!(bad&p[k]));if(V.length<2)continue;two++;
  for(const x of V)for(const y of V)if(x<y){
   const q=p[x]|p[y];
   if(M.some(m=>(m&H)&&(m&q)===q&&[a,b,c].filter(k=>m&p[k]).length>=2))continue;
   E.push([i,j,a,b,c,x,y]);
  }
 }
}
assert.deepEqual([coarse,one,two],[1417,40,3]);
assert.deepEqual(E,[[0,4,8,16,28,20,32],[0,4,8,20,28,16,32]]);
const weights=[3,3,2,2,2,1,1],edges=[];
for(const I of E){
 assert(M.every(m=>I.reduce((a,i,k)=>a+((m&p[i])?weights[k]:0),0)<=8));
 for(let b=0;b<128;b++){
  const X=weights.map((w,i)=>(b>>>i&1)?-w:w);
  assert.equal(X.reduce((a,b)=>a+b*b,0),32);
  const P=[];
  for(let i=0;i<2216;i++){
   const v=D.sparse.data.slice(i*33,(i+1)*33);
   assert(v.every(x=>x===0||Math.abs(x)===2));assert.equal(v.filter(x=>x).length,8);
   if(I.reduce((a,t,k)=>a+v[t]*X[k],0)>16)P.push(i);
  }
  assert.equal(P.length,14);edges.push(P);
 }
}
assert.equal(new Set(edges.flat()).size,448);
const pair=new Int16Array(2216).fill(-1),left=new Int16Array(256).fill(-1);
function aug(i,seen){
 for(const j of edges[i])if(!seen[j]){
  seen[j]=1;
  if(pair[j]===-1||aug(pair[j],seen)){pair[j]=i;left[i]=j;return true;}
 }
 return false;
}
for(let i=0;i<256;i++)assert(aug(i,new Uint8Array(2216)));
assert.equal(new Set(left).size,256);
for(let i=0;i<256;i++)assert(edges[i].includes(left[i]));
return {carriers:E,sparseBlockers:448,match:Array.from(left)};
}

const b=fs.readFileSync(process.argv[2]);assert.equal(crypto.createHash('sha1').update('blob '+b.length+'\0').update(b).digest('hex'),'a282e197b71dc04ee5295ebfbeafa223fd51cac0');console.log(JSON.stringify(scope(decode(b))));