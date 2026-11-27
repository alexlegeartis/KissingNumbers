// Exact union matching witness, 29D. No solver or search.
// Usage: node verify-union.mjs official-v2.npz official-v1-alt.npz own-exchange.npz own-grid-trade.npz matching-edits.json [actual-own-232157.npz]
import {readFileSync} from "node:fs"; import {createHash} from "node:crypto"; import {inflateRawSync} from "node:zlib";
const [L0,L1,L2,L3,EDITFILE,AACT]=process.argv.slice(2); if(!EDITFILE)throw Error("usage in source header");
const hash=x=>createHash("sha256").update(x).digest("hex");
const wanted=["b27ec26a1d186ed6e9b6885921abf3b2e6bbde132206fde2689dfdaa7e39f824","ded2da5a8141af3606d375df464f2ce3bfa0c45bb02204610d94d19ccdcf8cd5","130c0c144b014e1afa2284242718918fa7ccc2a9af8dd026a10d1a6a5aa9bd3e","d85b5cdf65eb128ad55104711cfc57621554de1d31005ad181757b46b2f3f9cd"];
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
 let C={"<i2":Int16Array,"<i4":Int32Array,"<i8":BigInt64Array,"|i1":Int8Array,"<u2":Uint16Array,"<u4":Uint32Array}[dtype];
 if(!C)throw Error("missing dtype "+dtype);
 const a=new C(data.buffer,data.byteOffset+end,(data.length-end)/C.BYTES_PER_ELEMENT);
 if(a.length!==shape.reduce((s,k)=>s*k,1))throw Error("shape "+[name,head]);
 out[name.replace(".npy","")]={shape,dtype,a};
 p+=46+k+u16(p+30)+u16(p+32);
}
return out;
};

const L=[L0,L1,L2,L3].map((p,i)=>{let b=readFileSync(p);if(hash(b)!==wanted[i])throw Error("pinned input hash "+i);const q=decode(b,inflateRawSync);for(const d of Object.values(q))if(d.dtype==='<'+'i8')d.a=Array.from(d.a,Number);return q;});
const makeStock=(D,check=true)=>{
const I=[0,2,4,5,6,10,11,23],J=Array.from({length:24},(_,i)=>i).filter(x=>!I.includes(x));
const pop=x=>{let s=0;for(;x;x&=x-1)s++;return s;};
const norm=a=>a.reduce((s,x)=>s+x*x,0),dot=(a,b)=>{let s=0;for(let i=0;i<a.length;i++)s+=a[i]*b[i];return s;};
const zz=()=>Array(24).fill(0),ls=[];
const lex=(a,b)=>{for(let j=0;j<a.length;j++)if(a[j]!==b[j])return a[j]-b[j];return 0;};
const G=[0],g=(1<<11)|(1<<10)|(1<<6)|(1<<5)|(1<<4)|(1<<2)|1;
for(let i=0;i<12;i++){
 let c=(g<<i)|(1<<23);
 G.push(...G.slice().map(v=>v^c));
}
for(let a=0;a<24;a++)for(let b=a+1;b<24;b++)for(let s of [-8,8])for(let t of [-8,8]){let v=zz();v[a]=s;v[b]=t;ls.push(v);}
for(let w of G)if(pop(w)===8){
 let ss=[];for(let i=0;i<24;i++)if(w>>i&1)ss.push(i);
 for(let b=0;b<256;b++)if(!(pop(b)&1)){
  const v=zz();for(let i=0;i<8;i++)v[ss[i]]=b>>i&1?-4:4;
  ls.push(v);
 }
}
for(let w of G){
 let x=Array.from({length:24},(_,i)=>w>>i&1?2:-2);
 for(let j=0;j<24;j++){
  let v=x.slice();v[j]*=-3;ls.push(v);
 }
}
ls.sort(lex);if(ls.length!==196560)throw Error("Leech count");
const added=[];
for(let i=0;i<D.added_integer.shape[0];i++){
 const x=Array.from(D.added_integer.a.slice(i*29,(i+1)*29));
 if(norm(x)!==128)throw Error("off norm");
 added.push(x);
}
const deleted=new Set(D.deleted_leech_indices.a);
if(check&&deleted.size!==32198)throw Error("deletion");
return {ls,added,deleted};};

const norm=x=>x.reduce((t,z)=>t+z*z,0), dot=(x,y)=>x.reduce((s,q,i)=>s+q*y[i],0);
const D=makeStock(L[0]);
const ar=(l,n=29)=>Array.from({length:l.length/n},(_,i)=>l.slice(i*n,(i+1)*n));
const kept=(N,del)=>ls.flatMap((v,i)=>del.has(i)?[]:[v.concat([0,0,0,0,0])]).concat(N);
const {ls}=D;
const B=kept(D.added,D.deleted), AH=kept(ar(L[1].added_integer.a),new Set(L[1].deleted_leech_indices.a));
const P=L[2],rs=new Set(ar(P.removed.a).map(x=>x.join()));
let A0=AH.filter(x=>!rs.has(x.join())).concat(ar(P.added.a));
if(A0.length!==232154)throw Error('trade1 size');
const Q=L[3],r2=ar(Q.rays.a,13).map((a,i)=>{let v=Array(29).fill(0),t=Math.sqrt(128/Q.squared_norms.a[i]);a.forEach((x,j)=>v[Q.coordinates.a[j]]=x*t);return v;});
const removed2=[], A=A0.filter(x=>{if(r2.some(a=>dot(x,a)>64)){removed2.push(x);return false;}return true;}).concat(r2);
if(A.length!==232157||B.length!==232230||removed2.length!==11)throw Error('count');
const akeys=new Set(A.map(x=>x.join())),bkeys=new Set(B.map(x=>x.join())),AX=A.filter(x=>!bkeys.has(x.join())),BX=B.filter(x=>!akeys.has(x.join()));
if(akeys.size!==A.length||bkeys.size!==B.length||BX.length-AX.length!==73)throw Error('unique');

const keybody=x=>x.slice(0,24).join(), X=new Map();
for(let j=0;j<BX.length;j++){let k=keybody(BX[j]);if(!X.has(k))X.set(k,[]);X.get(k).push(j);}
const ma=AX.map(()=>-1),mb=BX.map(()=>-1);
for(let i=0;i<AX.length;i++)for(let j of X.get(keybody(AX[i]))??[])if(mb[j]<0&&dot(AX[i],BX[j])>64){ma[i]=j;mb[j]=i;break;}
const seeds=ma.filter(x=>x>=0).length;
const edits=JSON.parse(readFileSync(EDITFILE,'utf-8'));
if(edits.length!==2675||new Set(edits.map(p=>p[0])).size!==2675)throw Error("edited pair count");
for(let [i,j] of edits){if(i<0||j<0||i>=AX.length||j>=BX.length)throw Error("pair bounds");ma[i]=j;}
const R=new Set(),hist={};
for(let i=0;i<AX.length;i++){
  let j=ma[i],u=dot(AX[i],BX[j]??[]);if(j<0||R.has(j)||!(u>64))throw Error("conflict matching at "+i+" "+j);
  R.add(j);hist[u]=(hist[u]??0)+1;
}
if(AX.length!==59911||BX.length!==59984||A.length-AX.length!==172246||seeds!==57385||R.size!==59911)throw Error("matching counts");
const lex29=(a,b)=>{for(let i=0;i<29;i++)if(a[i]!==b[i])return a[i]-b[i];return 0};
const canon=V=>{const a=Buffer.alloc(V.length*58);let k=0;
 for(let x of V.slice().sort(lex29))for(let y of x){a.writeInt16LE(y,k);k+=2;}return hash(a);};
const hA=canon(A), hB=canon(B);
if(canon(A0)!=="f62df9eb2b3cdf7442828722aeb8b7e66441ffc6c1a16043429f88bec0904a3a"||hA!=="85d48a37bccd4aa5950c3b5279b9a673511a71ebc1eef7f0df9647fa1d32495d")throw Error("own chain canonical hashes");
if(hB!=="3675c465c219c0e542f881356774f85dc9d67cb7ad0f83cfb805d0e808138233")throw Error("v0.2 official canonical set");
if(AACT){
 const exact=decode(readFileSync(AACT),inflateRawSync);
 const rawA=exact.points??exact.integer_points;
 if(rawA.shape.join()!=="232157,29")throw Error("reference shape");
 const publishedRows=Buffer.alloc(rawA.a.length*2);
 for(let i=0;i<rawA.a.length;i++)publishedRows.writeInt16LE(Number(rawA.a[i]),2*i);
 if(hash(publishedRows)!=="80aa09bd7b78170e2693cf221980eec312a7ac54e619218d3d5eb90b1fd7eecf")throw Error("published row order hash");
 if(canon(Array.from({length:232157},(_,i)=>Array.from(rawA.a.slice(i*29,(i+1)*29),Number)))!==hA)throw Error("published A canonical set");
}
console.log(JSON.stringify({common:172246,Aonly:59911,Bonly:59984,matching:R.size,unusedB:59984-R.size,hA,hB,hist,unionBest:232230}));
