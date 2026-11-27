'use strict';
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
const S=U.map(m=>A.filter(i=>m&P[i]));return S;
const oct=[];
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

function makeCharts(D,S){
 const A=Array.from({length:33},(_,i)=>i),P=A.map(i=>1n<<BigInt(i)),B=D.supports.data.map(BigInt),VV=new Map();
 for(const s of S)for(const j of A)if(!s.includes(j)){const X=[...s,j].sort((a,b)=>a-b);VV.set(X.join(','),X);}
 assert.equal(VV.size,420);
 const R=[],H={};
 for(const X of VV.values()){const support=X.reduce((m,i)=>m|P[i],0n),F=[...new Set(B.map(s=>s&support).filter(s=>pcB(s)===5))];
  for(let a=0;a<10;a++)for(let b=a+1;b<10;b++)for(let c=b+1;c<10;c++){
   const V=[X[a],X[b],X[c]],m=P[X[a]]|P[X[b]]|P[X[c]],Q=new Set();
   if(F.some(f=>pcB(f&m)<2))continue;
   for(const f of F){const q=f&m;if(pcB(q)===2)Q.add(Number(q));}
   if(V.includes(32))Q.add(Number(m^P[32]));
   if(Q.size>=3)continue;
   const e=Q.size;if(!H[e])H[e]=[0,0];H[e][V.includes(32)?1:0]++;
   R.push({H:X.filter(i=>!V.includes(i)),L:V,con:[...Q],support:X}); 
  }
 }
 assert.deepEqual(H,{0:[10,0],1:[110,20],2:[23,8]});
 return R;
}
function rows(D,k){return Array.from({length:k},(_,i)=>D.sparse.data.slice(33*i,33*i+33));}
function phaseCertificate(D,R){
 const charts=R.filter(f=>f.con.length===0),old=rows(D,2348);
 assert.equal(charts.length,10);
 const poss=Array.from({length:2216},()=>[]),meta=[];
 for(const f of charts){
  const U=f.H.reduce((s,x)=>s|(1n<<BigInt(x)),0n),V=f.L.reduce((s,x)=>s|(1n<<BigInt(x)),0n);
  for(let h=0;h<128;h++){
   const n=f.H.reduce((s,x,j)=>s|((h>>>j&1)?1n<<BigInt(x):0n),0n),id=meta.length;
   const E=[];
   for(let j=0;j<2216;j++){
    const r=old[j],a=f.H.reduce((s,x,i)=>s+((h>>>i&1)?-r[x]:r[x]),0),k=f.L.filter(x=>r[x]!==0).length;
    if(a>=12||a===10&&k<=1){E.push(j);poss[j].push(id);}
   }
   assert([2,3,4].includes(E.length));meta.push({U,V,n,w:12/E.length});
  }
 }
 function forced(a,b){
  const HH=a.U&b.U,h=pcB(HH)-2*pcB(HH&(a.n^b.n)),k1=pcB(a.U&b.V),k2=pcB(a.V&b.U),l=pcB(a.V&b.V),s=h-4;
  if(s<=0)return false;
  if(l===0||k1&&k2){
   const strict=!!l||k1>0&&k1<3||k2>0&&k2<3,d=s*s-k1-k2;
   return d>=0&&(d*d>4*k1*k2||(d*d===4*k1*k2&&strict));
  } else if(k1||k2) return s*s>=Math.max(k1,k2)+1;
  else return s>1||s===1&&l<3;
 }
 const counts={};
 for(const L of poss)if(L.length){
  let best=0;
  for(let m=1;m<1<<L.length;m++){
   const I=L.filter((x,j)=>m>>>j&1);if(I.some((i,t)=>I.slice(t+1).some(j=>forced(meta[i],meta[j]))))continue;
   best=Math.max(best,I.reduce((s,i)=>s+meta[i].w,0));
  }
  assert(best<=12);const key=[L.length,best];counts[key]=(counts[key]||0)+1;
 }
 assert.deepEqual(counts,{'2,8':384,'6,12':320,'6,8':72,'8,12':64});
 return {counts,charts};
}

function main(buffer){
 assert(Buffer.isBuffer(buffer));
 const sha=crypto.createHash('sha1').update('blob '+buffer.length+'\0').update(buffer).digest('hex');
 assert.equal(sha,'a282e197b71dc04ee5295ebfbeafa223fd51cac0');

 const D=decode(buffer),R=makeCharts(D,makeNineData(D));
 const {counts,charts}=phaseCertificate(D,R);
 console.log(JSON.stringify({counts,charts,allCharts:R.length}));
}
main(fs.readFileSync(process.argv[2]));