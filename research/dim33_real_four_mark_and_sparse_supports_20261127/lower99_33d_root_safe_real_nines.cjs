'use strict';
// All exact for the included proof; Node >=20.
// Source: verifications/improved/dim33-root-star-exchange/witness.npz
// Pass the literal witness from its published 370556 revision as sole argument.
const assert=require('node:assert/strict'), fs=require('node:fs'), crypto=require('node:crypto');
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


function makeNineData(D){
const P=Array.from({length:33},(_,i)=>1n<<BigInt(i)), B=D.supports.data.map(BigInt),A=P.map((p,i)=>i), C=new Set();
assert.equal(B.length,1837);assert(B.every(m=>pcB(m)===8));assert(D.parity.data.every(x=>x===0));
for(const m of B){const V=P.filter(p=>m&p);for(let a=0;a<8;a++)for(let b=a+1;b<8;b++)
 for(let c=b+1;c<8;c++)for(let d=c+1;d<8;d++)for(let e=d+1;e<8;e++)C.add(V[a]|V[b]|V[c]|V[d]|V[e]);}
let U=[];for(let a=0;a<29;a++)for(let b=a+1;b<30;b++)for(let c=b+1;c<31;c++)for(let d=c+1;d<32;d++)for(let e=d+1;e<33;e++){
 const m=P[a]|P[b]|P[c]|P[d]|P[e];if(!C.has(m))U.push(m);}
const counts=[U.length];
for(let k=6;k<=9;k++){
 const prev=new Set(U),V=new Set();
 for(const m of U)for(const p of P)if(p>m&&P.every(q=>!(m&q)||prev.has((m^q)|p)))V.add(m|p);
 U=[...V];counts.push(U.length);
}
assert.deepEqual(counts,[134464,79706,3044,197,18]);
assert(U.every(m=>m&P[32]));
assert.deepEqual(D.sparse.shape,[4348,33]);
const S=U.map(m=>A.filter(i=>m&P[i])),oct=[];
for(let i=0;i<2216;i++){
 const row=D.sparse.data.slice(33*i,33*i+33);assert(row.every(v=>v===0||Math.abs(v)===2));assert.equal(row.filter(v=>v).length,8);
 oct.push({pos:row.reduce((m,v,k)=>v===2?m|P[k]:m,0n),neg:row.reduce((m,v,k)=>v===-2?m|P[k]:m,0n)});
}
const R=[],rev=Array.from({length:2216},()=>[]),dist={};
for(const I of S)for(let s=0;s<512;s++){
 const pos=I.reduce((m,i,k)=>!(s>>>k&1)?m|P[i]:m,0n),neg=I.reduce((m,i,k)=>(s>>>k&1)?m|P[i]:m,0n),bad=[];
 for(let j=0;j<oct.length;j++){
  const o=oct[j];
  if(pcB((o.pos|o.neg)&(pos|neg))-2*pcB((o.pos&neg)|(o.neg&pos))>=5){bad.push(j);rev[j].push(R.length);}
 }
 dist[bad.length]=(dist[bad.length]||0)+1;R.push({pos:Number(pos),neg:Number(neg),B:bad});
}
assert.equal(R.length,9216);assert.equal(Math.min(...R.map(x=>x.B.length)),8);
return {S,R,rev,priced:rev.filter(z=>z.length).length};
}
function pcB(x){let n=0;for(;x;x&=x-1n)n++;return n}
function coreNine(T,p,UCheck=true){
const R=T.R.flatMap((z,i)=>z.B.length<=p?[{...z,id:i}]:[]);const rev=T.rev.map(u=>u.filter(i=>T.R[i].B.length<=p));const n=R.length,back=new Int32Array(T.R.length).fill(-1),edges=R.map(()=>[]);
R.forEach(({id},i)=>back[id]=i);const pc=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return Math.imul((x+(x>>>4))&0xf0f0f0f,0x01010101)>>>24};
const q=(a,b)=>pc((a.pos|a.neg)&(b.pos|b.neg))-2*pc((a.pos&b.neg)|(a.neg&b.pos))+(((a.pos/2**32|0)&(b.pos/2**32|0))|((a.neg/2**32|0)&(b.neg/2**32|0)))-(((a.pos/2**32|0)&(b.neg/2**32|0))|((a.neg/2**32|0)&(b.pos/2**32|0)));
let ne=0;
for(let i=0;i<n;i++){
 const a=R[i],C=new Int8Array(n);
 for(const j of a.B)for(const id of rev[j])C[back[id]]++;
 for(let j=i+1;j<n;j++)if(a.B.length+R[j].B.length-C[j]<=p&&q(a,R[j])<=5){edges[i].push(j);edges[j].push(i);ne++;}
}
const D=Int32Array.from(edges,v=>v.length),alive=new Uint8Array(n).fill(1),Q=[];for(let i=0;i<n;i++)if(D[i]<p){Q.push(i);alive[i]=0;}
for(let k=0;k<Q.length;k++)for(const j of edges[Q[k]])if(alive[j]&&--D[j]<p){Q.push(j);alive[j]=0;}
const V=R.flatMap((z,i)=>alive[i]?[i]:[]); if(!UCheck)return {p,n,ne,core:V.length,valid:false,calls:0};
const ind=new Uint32Array(2216),Vset=new Set(V);
const e=edges.map(a=>new Set(a.filter(i=>Vset.has(i))));
V.sort((i,j)=>D[j]-D[i]||i-j);
let calls=0;
function run(P,total,k){
 calls++;
 if(k>p)return true;
 if(k+P.length<=p)return false;
 for(let q=0;q<P.length;q++){
  const i=P[q],add=R[i].B.filter(j=>!ind[j]),z=total+add.length;
  if(z>p)continue;
  add.forEach(a=>ind[a]=1);
  const rest=P.slice(q+1).filter(j=>e[i].has(j)&&z+R[j].B.filter(t=>!ind[t]).length<=p);
  const ret=run(rest,z,k+1);
  add.forEach(a=>ind[a]=0);
  if(ret)return true;
  if(k+P.length-q-1<=p)break;
 }
 return false;
}
return {p,n,ne,core:V.length,valid:run(V,0,0),calls};
}
function main(buffer){
 assert(Buffer.isBuffer(buffer));
 const sha=crypto.createHash('sha1').update('blob '+buffer.length+'\0').update(buffer).digest('hex');
 assert.equal(sha,'a282e197b71dc04ee5295ebfbeafa223fd51cac0');
 const D=decode(buffer),T=makeNineData(D);
 assert.equal(T.priced,1864);
 const extensions=new Map();
 for(const I of T.S)for(let i=0;i<33;i++)if(!I.includes(i)){
   const key=[...I,i].sort((a,b)=>a-b).join(',');
   extensions.set(key,(extensions.get(key)||0)+1);
 }
 const exHist={};for(const d of extensions.values())exHist[d]=(exHist[d]||0)+1;
 assert.equal(extensions.size,420);assert.deepEqual(exHist,{1:414,2:2,3:2,4:2});
 assert(16*199<57*57);assert(64*199<117*117); // exact inequalities for constructive 10-profile
 const output=[];
 for(let p=8;p<=23;p++){
  const result=coreNine(T,p);
  assert(!result.valid);
  output.push(result);
 }
 console.log(JSON.stringify({supports:T.S.length,candidates:T.R.length,prices:T.priced,tenSupports:extensions.size,tenExtensions:exHist,budgets:output}));
}
main(fs.readFileSync(process.argv[2]));
