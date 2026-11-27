// Usage: node verify-radial.mjs <official v0.2 d29_certificate.npz>
import {readFileSync} from "node:fs"; import {createHash} from "node:crypto"; import {inflateRawSync} from "node:zlib";
const INPUT=readFileSync(process.argv[2]);if(createHash("sha256").update(INPUT).digest("hex")!=="b27ec26a1d186ed6e9b6885921abf3b2e6bbde132206fde2689dfdaa7e39f824")throw Error("input SHA");
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

const X=[[0,2,0,2,0,0,0,2,2,2,0,0,-2,2,2,-2,-2,2,2,2,2,2,-2,0,-8,0,0,0,0],[0,4,0,0,0,0,0,0,0,0,0,0,0,-4,0,0,4,0,0,0,0,0,-4,0,-8,0,0,0,0],[0,4,0,0,0,0,0,0,0,0,0,0,4,0,4,0,0,0,0,0,0,0,4,0,-8,0,0,0,0],[4,2,4,0,0,0,0,2,2,-2,-4,0,0,0,-2,2,-2,0,0,0,-2,0,0,-4,0,0,0,-4,-4]], B=makeStock(decode(INPUT,inflateRawSync)), keep=B.ls.filter((_,i)=>!B.deleted.has(i)).map(z=>z.concat([0,0,0,0,0])).concat(B.added);
function U(n){let a=[];for(let b=1;b<32;b++){const v=Array.from({length:5},(_,j)=>(b>>j&1)?j:-1).filter(x=>x>=0);if(v.length!==n)continue;
for(let k=0;k<(1<<n);k++){let z=Array(5).fill(0);v.forEach((i,j)=>z[i]=k>>j&1?-1:1);a.push(z);}
}return a;}
const T=[1,2,4,5].flatMap(U);
const dp=(x,y,k)=>{let d=0;for(let j=0;j<k;j++)d+=x[j]*y[j];return d;};

const gcd=(a,b)=>{a=Math.abs(a);b=Math.abs(b);while(b)[a,b]=[b,a%b];return a};
function vertices(v){
 let V=[[0,0,1],[16,0,1],[16,16,1],[0,16,1]];
 const check=(u)=>{if(u.some(x=>!Number.isSafeInteger(x)))throw Error('exact overflow')};
 for(let [a,b] of v){
  let W=[];const D=V.map(([x,y,s])=>a*x+b*y-64*s);
  for(let j=0;j<V.length;j++){
   let k=(j+1)%V.length, r=D[j],s=D[k];if(r<=0)W.push(V[j]);
   if((r<0&&s>0)||(r>0&&s<0)){
    let t=V[j].map((x,i)=>s*x-r*V[k][i]);check(t);let q=t.reduce(gcd,0)*Math.sign(t[2]);if(!q)throw Error('zero ray');W.push(t.map(x=>x/q));
   }
  }
  if(!W.length)throw Error('no origin'); V=W;
 }
 return V;
}
const exactHist=[],global=[];
for(let i=0;i<X.length;i++){
 const h=X[i].slice(0,24), n=dp(h,h,24), P=new Map();
 for(let y of keep){let t=y.slice(24),key=t.join(),p=dp(h,y,24);if(!P.has(key)||P.get(key)[0]<p)P.set(key,[p,t]);}
 let no=0,on=0,maxn=0,maxd=1,traces=[], opts=[];
 for(let v of T){
  const m=dp(v,v,5),R=[];for(let [p,w] of P.values())R.push([p,dp(v,w,5)]);
  const Z=vertices(R);let best=0,den=1,tight=[];
  for(let u of Z){if(R.some(([p,q])=>p*u[0]+q*u[1]>64*u[2]))throw Error('invalid vertex');
   const q=n*u[0]**2+m*u[1]**2,d=u[2]**2;
   if(!Number.isSafeInteger(q)||!Number.isSafeInteger(d))throw Error('rational norm overflow');
   if(BigInt(q)*BigInt(den)>BigInt(best)*BigInt(d)){best=q;den=d;}
   if(q>=128*d){
    if(q>128*d||u[0]*u[1]!==0)opts.push([v,u,q,d]);
    else tight.push(u);
   }
  }
  if(best<128*den)no++;else if(best===128*den&&!opts.length)on++;
  if(BigInt(best)*BigInt(maxd)>BigInt(maxn)*BigInt(den)){maxn=best;maxd=den;}
  traces.push([v,best,den,tight.length]);
 }
 exactHist.push({scanned:T.length,norm:n,strict:no,axis:on,max:[maxn,maxd],opts});
}
if(exactHist.length!==4||exactHist.some((e,i)=>e.scanned!==162||e.strict!==162||e.opts.length||e.max[0]!==[80,80,80,2560][i]||e.max[1]!==[1,1,1,25][i]))throw Error("radial exact census");
console.log(JSON.stringify(exactHist));
