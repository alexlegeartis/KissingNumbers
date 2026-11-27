import {readFileSync} from "node:fs";
import {createHash} from "node:crypto";
import {inflateRawSync} from "node:zlib";
const file=process.argv[2];if(!file)throw Error("usage: node verify-bw29.mjs <pinned Takhanov d29_certificate.npz> [--hadamard]");
const bytes=readFileSync(file),want="b27ec26a1d186ed6e9b6885921abf3b2e6bbde132206fde2689dfdaa7e39f824";
if(createHash("sha256").update(bytes).digest("hex")!==want)throw Error("source hash");
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
const g=(D)=>{
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
if(deleted.size!==32198)throw Error("deletion");
const map=new Map(),rest=[];
function g1(b){
 const k=b.join();
 if(!map.has(k))map.set(k,{b,heads:[],l:[],ph:[],mid:[],orig:[]});
 return map.get(k);
}
for(let i=0;i<ls.length;i++){
 const x=ls[i],b=J.map(i=>x[i]);if(norm(b)!==64)continue;
 const a=I.map(i=>x[i]),r=g1(b);r.heads.push(a);
 if(!deleted.has(i)){r.l.push(i);r.orig.push(a.concat([0,0,0,0,0]));}
}
const short=[],aux=[],poles=[];
for(let i=0;i<added.length;i++){
 const x=added[i],b=J.map(i=>x[i]),h=I.map(i=>x[i]),t=x.slice(24),v=h.concat(t);
 let n=norm(b),s=norm(t);
 if(n===64){
   const gg=g1(b);
   gg.orig.push(v);
   if(s===32)gg.mid.push(i);else if(s===64){gg.ph.push(i);aux.push(i);}else throw Error("body 64 shape");
 }else if(n===32&&s===32)short.push(i);
 else if(n===0&&s===128)poles.push(i);
 else if(n!==0)rest.push(i);
}
if(map.size!==4368||aux.length!==2416||short.length!==580||poles.length!==40||[...map.values()].some(x=>![0,16].includes(x.heads.length)))throw Error("projection inventory");
const S=[...map.values()].filter(x=>x.heads.length===16).sort((a,b)=>lex(a.b,b.b));
const Lx=new Map(S.map((v,i)=>[v.b.join(),i]));
const hist={},by={};
for(let i=0;i<S.length;i++){
 const s=S[i], n=s.orig.length;
 if(n>26)throw Error("orthoplex");
 const P=new Set(s.orig.map(u=>u.join()));
 const ant=s.orig.filter(v=>P.has(v.map(x=>-x).join())).length;
 hist[n]=(hist[n]??0)+1;
 const h=[s.l.length,s.mid.length,s.ph.length,ant].join();
 by[h]=(by[h]??0)+1;
 for(let j=0;j<n;j++)for(let k=0;k<j;k++)if(dot(s.orig[j],s.orig[k])>0)throw Error("base local");
 s.axis=[];for(let a=0;a<5;a++)if(s.orig.every(v=>v[8+a]===0))s.axis.push(a);
 s.ant=ant;
 for(let j=0;j<i;j++)if(dot(s.b,S[j].b)>32)throw Error("BW bound");
}
const axes={},axisA={};
for(let i=0;i<S.length;i++){
 const u=S[i];const h=[u.orig.length,u.axis.length,u.ant].join();axes[h]=(axes[h]??0)+1;
}
const auxShape={};for(let i of aux){let t=added[i].slice(24),k=t.map(Math.abs).sort((a,b)=>a-b).join();auxShape[k]=(auxShape[k]??0)+1;}
return {I,J,ls,deleted,added,S,short,aux,poles,rest,Lx,hist,by,axes,auxShape};
};
const ax=(g)=>{
const G=g.added, frame=g.S, norm=x=>x.reduce((s,a)=>s+a*a,0), dot=(x,y)=>x.reduce((s,a,i)=>s+a*y[i],0), T=G.map(x=>x.slice(24)),J=G.map(x=>g.J.map(i=>x[i]));
const axial=Array.from({length:10},()=>[]);
for(let j of g.aux)if(T[j].some(x=>Math.abs(x)===8)){
 let i=T[j].findIndex(x=>x!==0);if(norm(T[j])!==64)throw Error("tail");
 axial[2*i+(T[j][i]>0?1:0)].push(j);
}
const cdeg=[];
for(const L of axial){
 if(L.length!==32)throw Error("axis cap");
 const s=new Set(L.map(i=>J[i].join()));
 if(s.size!==32||L.some(i=>!s.has(J[i].map(x=>-x).join())))throw Error("axis antipodes");
 for(let i of L)for(let j of L)if(i<j&&![0,-64].includes(dot(J[i],J[j])))throw Error("axis frame");
}
const hist={}, pass=[],forms={};
for(let i=0;i<frame.length;i++){
 const f=frame[i]; if(f.orig.length!==24)continue;
 if(f.axis.length!==1||f.ant!==24)throw Error("fibre24");
 const a=f.axis[0];
 const p=axial.slice(2*a,2*a+2).map(z=>Math.max(...z.map(j=>dot(J[j],f.b))));
 const key=p.slice().sort((x,y)=>x-y).join();
 hist[key]=(hist[key]??0)+1;
 if((64-p[0])**2+(64-p[1])**2>=4096)pass.push({i,a,b:f.b,p});
 forms[a]=(forms[a]??0)+1;
}
return {hist,pass,forms};
};
const diag=(g,ax)=>{
const R=g.added.map(x=>g.I.map(i=>x[i]).concat(x.slice(24))),J=g.added.map(x=>g.J.map(i=>x[i]));
const rd=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
let count=0,by=Array.from({length:5},(_,j)=>g.added.map((_,i)=>i).filter(i=>Math.abs(g.added[i][24+j])===4));
for(const e of ax.pass){
 const f=g.S[e.i],A=by[e.a].filter(i=>rd(f.b,J[i])===32);
 for(const r of f.orig.filter(x=>x.find(k=>k)>0)){
  let mask=0;
  for(let i of A){
   const q=rd(r,R[i]);if(Math.abs(q)!==32)continue;
   const t=(q>0?1:0)+(g.added[i][24+e.a]>0?2:0);mask|=1<<t;
   if(mask===15)break;
  }
  if(mask!==15)throw Error("missing exact diagonal wall");
  count++;
 }
}
if(count!==4020)throw Error("angular count");
return count;
};
const polish=(g)=>{
const tRows=[],J=[];
const dot=(a,b)=>{let s=0;for(let i=0;i<a.length;i++)s+=a[i]*b[i];return s;};
for(const x of g.added){tRows.push(x.slice(24));J.push(g.J.map(i=>x[i]));}
const groups=[];let key=new Map();
for(let i=0;i<tRows.length;i++){
 const t=tRows[i];if(!t.some(Boolean))continue;
 const k=t.join();
 if(!key.has(k)){key.set(k,groups.length);groups.push({t,L:[]});}
 groups[key.get(k)].L.push(i);
}
const F=[],hist={};
for(let i=0;i<g.S.length;i++){
 const f=g.S[i];if(f.orig.length>=24)continue;
 const W=[];
 for(const {t,L} of groups){
  let q=-Infinity,owner=-1;
  for(let j of L){
   const at=dot(f.b,J[j]);if(at>q){q=at;owner=j;}
  }
  const rhs=64-q;if(rhs<0)throw Error("body dot");
  // If rhs >= max on norm8 sphere the constraint is unnecessary.
  if(rhs*rhs<64*dot(t,t))W.push([t,rhs,owner]);
 }
 let a=W.filter(w=>w[1]===0).length,k=[f.orig.length,f.ant,a].join();
 hist[k]=(hist[k]??0)+1;F.push({i,b:f.b,orig:f.orig.length,ant:f.ant,W,l:f.l,mid:f.mid,ph:f.ph});
}
return {F,hist,rows:groups.length};
};
const maxpol=(all)=>{
const gcd=(a,b)=>{a=Math.abs(a);b=Math.abs(b);while(b){let c=a%b;a=b;b=c;}return a;};
const dot=(a,v)=>a.reduce((s,x,i)=>s+x*v[i],0);
function canonical(a){let D=a.reduce(gcd,0);if(!D||a[5]<=0)throw Error("homogeneous");return a.map(x=>x/D);}
function verts(W){
 let P=Array.from({length:32},(_,k)=>{
  const a=Array.from({length:5},(_,i)=>k>>i&1?8:-8).concat(1);
  let M=0n;for(let i=0;i<5;i++)M|=1n<<BigInt(2*i+(k>>i&1));
  return {a,M};
 });
 const ww=W.slice().sort((u,v)=>u[1]/Math.sqrt(dot(u[0],u[0]))-v[1]/Math.sqrt(dot(v[0],v[0])));
 let maximum=P.length,made=0;
 for(let s=0;s<ww.length;s++){
  const [b,r0] = ww[s], h=1n<<BigInt(10+s);
  const Q=P.map(({a})=>dot(b,a)-r0*a[5]);
  let low=[],high=[];
  for(let j=0;j<P.length;j++){if(Q[j]<0)low.push(j);if(Q[j]>0)high.push(j);}
  if(!high.length){for(let j=0;j<P.length;j++)if(!Q[j])P[j]={a:P[j].a,M:P[j].M|h};continue;}
  const B=new Map();
  for(let j=0;j<P.length;j++)if(Q[j]<=0){
   const w=P[j];B.set(w.a.join(),{a:w.a,M:w.M|(!Q[j]?h:0n)});
  }
  let C=new Map();
  for(let j of low)for(let k of high){
   const x=P[j], y=P[k], M=x.M&y.M;
   let c=C.get(M);
   if(c===undefined){c=0;for(const p of P)if((p.M&M)===M){c++;if(c>=3)break;}C.set(M,c);}
   if(c!==2)continue;
   const raw=x.a.map((a,i)=>Q[k]*a-Q[j]*y.a[i]);
   if(raw.some(a=>!Number.isSafeInteger(a)))throw Error("int precision");
   const a=canonical(raw),k0=a.join(),old=B.get(k0);
   B.set(k0,{a,M:(old?.M??0n)|M|h});made++;
  }
  P=[...B.values()];
  maximum=Math.max(maximum,P.length);
  if(!P.length)throw Error("origin always feasible");
 }
 for(const {a} of P)if(W.some(([q,t])=>dot(q,a)>t*a[5]))throw Error("invalid vertices");
 return {P:P.map(w=>w.a),maximum,made};
}
const rows=[], hist={},wit=[], near=[];
for(const f of all.F){
 const V=verts(f.W),mm=V.P.map(a=>[dot(a.slice(0,5),a),a[5]**2,a]);
 let max=mm.reduce((z,x)=>x[0]*z[1]>z[0]*x[1]?x:z);
 const k=max.slice(0,2).map(x=>x/gcd(max[0],max[1])).join("/");
 hist[k]=(hist[k]??0)+1;
 const row={...f,count:V.P.length,peak:V.maximum,made:V.made,norm:k,max:max[2],P:V.P};
 rows.push(row);
 if(64*max[1]<=max[0])wit.push(row);
 else near.push(row);
}
return {rows,hist,wit,near};
};
const mix=(g)=>{
const gcd=(a,b)=>{a=Math.abs(a);b=Math.abs(b);while(b){let c=a%b;a=b;b=c;}return a;};
const dot=(a,v)=>a.reduce((s,x,i)=>s+x*v[i],0);
function canonical(a,d){let D=a.reduce(gcd,0);if(!D||a[d]<=0)throw Error("homogeneous");return a.map(x=>x/D);}
function verts(W,d){
 let P=Array.from({length:1<<d},(_,k)=>{
  const a=Array.from({length:d},(_,i)=>k>>i&1?8:-8).concat(1);
  let M=0n;for(let i=0;i<d;i++)M|=1n<<BigInt(2*i+(k>>i&1));
  return {a,M};
 });
 const ww=W.slice().sort((u,v)=>u[1]/Math.sqrt(dot(u[0],u[0]))-v[1]/Math.sqrt(dot(v[0],v[0])));
 let maximum=P.length,made=0;
 for(let s=0;s<ww.length;s++){
  const [b,r0] = ww[s], h=1n<<BigInt(2*d+s);
  const Q=P.map(({a})=>dot(b,a)-r0*a[d]);
  let low=[],high=[];
  for(let j=0;j<P.length;j++){if(Q[j]<0)low.push(j);if(Q[j]>0)high.push(j);}
  if(!high.length){for(let j=0;j<P.length;j++)if(!Q[j])P[j]={a:P[j].a,M:P[j].M|h};continue;}
  const B=new Map();
  for(let j=0;j<P.length;j++)if(Q[j]<=0){
   const w=P[j];B.set(w.a.join(),{a:w.a,M:w.M|(!Q[j]?h:0n)});
  }
  let C=new Map();
  for(let j of low)for(let k of high){
   const x=P[j], y=P[k], M=x.M&y.M;
   let c=C.get(M);
   if(c===undefined){c=0;for(const p of P)if((p.M&M)===M){c++;if(c>=3)break;}C.set(M,c);}
   if(c!==2)continue;
   const raw=x.a.map((a,i)=>Q[k]*a-Q[j]*y.a[i]);
   if(raw.some(a=>!Number.isSafeInteger(a)))throw Error("int precision");
   const a=canonical(raw,d),k0=a.join(),old=B.get(k0);
   B.set(k0,{a,M:(old?.M??0n)|M|h});made++;
  }
  P=[...B.values()];
  maximum=Math.max(maximum,P.length);
  if(!P.length)throw Error("origin always feasible");
 }
 for(const {a} of P)if(W.some(([q,t])=>dot(q,a)>t*a[d]))throw Error("invalid vertices");
 return {P:P.map(w=>w.a),maximum,made};
}

const J=g.added.map(x=>g.J.map(i=>x[i])),R=g.added.map(x=>g.I.map(i=>x[i]).concat(x.slice(24)));
const rd=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
function Null(A,n){
 A=A.map(x=>x.slice());const piv=[];let r=0;
 for(let j=0;j<n;j++){
  let k=A.findIndex((x,i)=>i>=r&&x[j]);if(k<0)continue;
  [A[r],A[k]]=[A[k],A[r]];
  let p=A[r][j];if(p<0){A[r]=A[r].map(x=>-x);p=-p;}
  for(let i=0;i<A.length;i++)if(i!==r&&A[i][j]){
   const t=A[i][j];A[i]=A[i].map((x,k)=>x*p-t*A[r][k]);
   const d=A[i].reduce(gcd,0);if(d)A[i]=A[i].map(x=>x/d);
  }
  piv.push(j);r++;
 }
 const free=Array.from({length:n},(_,i)=>i).filter(i=>!piv.includes(i));
 const D=piv.reduce((x,j,i)=>x*A[i][j]/gcd(x,A[i][j]),1);
 const B=Array.from({length:n},(_,i)=>free.map((j,k)=>{
  let t=piv.indexOf(i);return t<0?(i===j?D:0):-D*A[t][j]/A[t][i];
 }));
 if(B.flat().some(x=>!Number.isSafeInteger(x)))throw Error("null basis");
 for(const v of A)for(let i=0;i<free.length;i++)if(rd(v,B.map(x=>x[i])))throw Error("null correction");
 return {B,D,free};
}
const rows=[],hist={},wit=[];
for(const f of g.S.map((v,i)=>({...v,i}))){
 if(f.orig.length>=24||f.l.length===16)continue;
 const ant=new Set(f.orig.map(a=>a.join()));
 const normals=f.orig.filter(a=>a.find(x=>x)>0&&ant.has(a.map(x=>-x).join()));
 const {B,D,free}=Null(normals,13),d=free.length;
 if(d<3||d>5)throw Error("unexpected kernel "+d);
 const proj=B.map((_,j)=>j);let W=new Map();
 for(let j=0;j<R.length;j++){
  const a=free.map((_,k)=>rd(R[j],B.map(x=>x[k]))),b=D*(64-rd(f.b,J[j]));
  if(!a.some(Boolean))continue;
  if(8*a.reduce((s,x)=>s+Math.abs(x),0)<=b)continue;
  const k=a.join(),p=W.get(k);
  if(!p||b<p[1])W.set(k,[a,b,j]);
 }
 let W0=[...W.values()];
 const P=verts(W0,d);
 const NN=P.P.map(a=>{const t=B.map(x=>rd(x,a)),n=rd(t,t),den=(D*a[d])**2;return [n,den,t,a[d]];});
 let X=NN.reduce((z,q)=>q[0]*z[1]>z[0]*q[1]?q:z),m=X.slice(0,2).map(c=>c/gcd(X[0],X[1])).join("/");
 hist[m]=(hist[m]??0)+1;
 if(X[0]>=64*X[1])wit.push({i:f.i,d,f:f.b,v:X[2],norm:m,B,D,W0,verts:P.P,newoptions:NN.filter(v=>64*v[1]<=v[0]).map(([n,d,t,z])=>[n,d,t])});
 rows.push({i:f.i,d,norm:m,w:W0.length,vs:P.P.length,max:Math.round(1000*X[0]/X[1])/1000});
}
return {rows,hist,wit};};
const fold=(g,base)=>{
const gcd=(a,b)=>{a=Math.abs(a);b=Math.abs(b);while(b){let c=a%b;a=b;b=c;}return a;};
const dot=(a,v)=>a.reduce((s,x,i)=>s+x*v[i],0);
function canonical(a,d){let D=a.reduce(gcd,0);if(!D||a[d]<=0)throw Error("homogeneous");return a.map(x=>x/D);}
function verts(W,d){
 let P=Array.from({length:1<<d},(_,k)=>{
  const a=Array.from({length:d},(_,i)=>k>>i&1?8:-8).concat(1);
  let M=0n;for(let i=0;i<d;i++)M|=1n<<BigInt(2*i+(k>>i&1));
  return {a,M};
 });
 const ww=W.slice().sort((u,v)=>u[1]/Math.sqrt(dot(u[0],u[0]))-v[1]/Math.sqrt(dot(v[0],v[0])));
 let maximum=P.length,made=0;
 for(let s=0;s<ww.length;s++){
  const [b,r0] = ww[s], h=1n<<BigInt(2*d+s);
  const Q=P.map(({a})=>dot(b,a)-r0*a[d]);
  let low=[],high=[];
  for(let j=0;j<P.length;j++){if(Q[j]<0)low.push(j);if(Q[j]>0)high.push(j);}
  if(!high.length){for(let j=0;j<P.length;j++)if(!Q[j])P[j]={a:P[j].a,M:P[j].M|h};continue;}
  const B=new Map();
  for(let j=0;j<P.length;j++)if(Q[j]<=0){
   const w=P[j];B.set(w.a.join(),{a:w.a,M:w.M|(!Q[j]?h:0n)});
  }
  let C=new Map();
  for(let j of low)for(let k of high){
   const x=P[j], y=P[k], M=x.M&y.M;
   let c=C.get(M);
   if(c===undefined){c=0;for(const p of P)if((p.M&M)===M){c++;if(c>=3)break;}C.set(M,c);}
   if(c!==2)continue;
   const raw=x.a.map((a,i)=>Q[k]*a-Q[j]*y.a[i]);
   if(raw.some(a=>!Number.isSafeInteger(a)))throw Error("int precision");
   const a=canonical(raw,d),k0=a.join(),old=B.get(k0);
   B.set(k0,{a,M:(old?.M??0n)|M|h});made++;
  }
  P=[...B.values()];
  maximum=Math.max(maximum,P.length);
  if(!P.length)throw Error("origin always feasible");
 }
 for(const {a} of P)if(W.some(([q,t])=>dot(q,a)>t*a[d]))throw Error("invalid vertices");
 return {P:P.map(w=>w.a),maximum,made};
}


const Jold=[],Hold=[], keep=[];
const norm=x=>x.reduce((s,a)=>s+a*a,0),rd=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
for(let i=0;i<g.ls.length;i++)if(!g.deleted.has(i)){const x=g.ls[i];const h=g.I.map(k=>x[k]);if(!h.some(Boolean))continue;Hold.push(h);Jold.push(g.J.map(k=>x[k]));keep.push(i);}
const hist={},rows=[],win=[],edge=[];
for(const o of base.wit){
 const W=new Map(o.W0.map(w=>[w[0].join(),w]));
 let ns=[...o.newoptions],aa=[];
 for(const [n,d,t] of ns){
  const blockers=[];for(let j=0;j<Hold.length;j++){
   const p=rd(o.f,Jold[j]),A=rd(t.slice(0,8),Hold[j]);
   if(A>0&&64*A*A>(64-p)**2*n){
    blockers.push(keep[j]);if(blockers.length>=10)break;
   }
  }
  aa.push({n,d,t,blockers});
 }
 aa.sort((a,b)=>a.blockers.length-b.blockers.length);
 const best=aa[0];
 if(!best.blockers.length)win.push({i:o.i,b:o.f,...best});
 if(best.blockers.length<=2)edge.push({i:o.i,b:o.f,found:aa.filter(v=>v.blockers.length<=2)});
 for(let j=0;j<Hold.length;j++){
  const a=Array.from({length:o.d},(_,k)=>rd(Hold[j],o.B.slice(0,8).map(x=>x[k]))),b=o.D*(64-rd(o.f,Jold[j]));
  if(b<0)throw Error("Leech equatorial midpoint");
  if(!a.some(Boolean)||8*a.reduce((s,x)=>s+Math.abs(x),0)<=b)continue;
  const k=a.join(),p=W.get(k);if(!p||b<p[1])W.set(k,[a,b,"eq"+keep[j]]);
 }
 const P=verts([...W.values()],o.d);
 const NN=P.P.map(a=>{const t=o.B.map(x=>rd(x,a));return [norm(t),(o.D*a[o.d])**2,t];});
 const X=NN.reduce((z,q)=>q[0]*z[1]>z[0]*q[1]?q:z),M=X.slice(0,2).map(c=>c/gcd(X[0],X[1])).join("/");
 hist[M]=(hist[M]??0)+1;
 for(const q of NN)if(q[0]>=64*q[1])win.push({i:o.i,b:o.f,n:q[0],d:q[1],t:q[2],full:true});
 rows.push({i:o.i,d:o.d,n:M,options:ns.length,minBlock:best.blockers.length,H:W.size,nd:P.P.length,max:X});
}
return {hist,rows,win,edge,kept:keep.length};
};
const paid=(g,base,near)=>{
const gcd=(a,b)=>{a=Math.abs(a);b=Math.abs(b);while(b){let c=a%b;a=b;b=c;}return a;};
const dot=(a,v)=>a.reduce((s,x,i)=>s+x*v[i],0);
function canonical(a,d){let D=a.reduce(gcd,0);if(!D||a[d]<=0)throw Error("homogeneous");return a.map(x=>x/D);}
function verts(W,d){
 let P=Array.from({length:1<<d},(_,k)=>{
  const a=Array.from({length:d},(_,i)=>k>>i&1?8:-8).concat(1);
  let M=0n;for(let i=0;i<d;i++)M|=1n<<BigInt(2*i+(k>>i&1));
  return {a,M};
 });
 const ww=W.slice().sort((u,v)=>u[1]/Math.sqrt(dot(u[0],u[0]))-v[1]/Math.sqrt(dot(v[0],v[0])));
 let maximum=P.length,made=0;
 for(let s=0;s<ww.length;s++){
  const [b,r0] = ww[s], h=1n<<BigInt(2*d+s);
  const Q=P.map(({a})=>dot(b,a)-r0*a[d]);
  let low=[],high=[];
  for(let j=0;j<P.length;j++){if(Q[j]<0)low.push(j);if(Q[j]>0)high.push(j);}
  if(!high.length){for(let j=0;j<P.length;j++)if(!Q[j])P[j]={a:P[j].a,M:P[j].M|h};continue;}
  const B=new Map();
  for(let j=0;j<P.length;j++)if(Q[j]<=0){
   const w=P[j];B.set(w.a.join(),{a:w.a,M:w.M|(!Q[j]?h:0n)});
  }
  let C=new Map();
  for(let j of low)for(let k of high){
   const x=P[j], y=P[k], M=x.M&y.M;
   let c=C.get(M);
   if(c===undefined){c=0;for(const p of P)if((p.M&M)===M){c++;if(c>=3)break;}C.set(M,c);}
   if(c!==2)continue;
   const raw=x.a.map((a,i)=>Q[k]*a-Q[j]*y.a[i]);
   if(raw.some(a=>!Number.isSafeInteger(a)))throw Error("int precision");
   const a=canonical(raw,d),k0=a.join(),old=B.get(k0);
   B.set(k0,{a,M:(old?.M??0n)|M|h});made++;
  }
  P=[...B.values()];
  maximum=Math.max(maximum,P.length);
  if(!P.length)throw Error("origin always feasible");
 }
 for(const {a} of P)if(W.some(([q,t])=>dot(q,a)>t*a[d]))throw Error("invalid vertices");
 return {P:P.map(w=>w.a),maximum,made};
}


const rd=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0),norm=x=>rd(x,x);
const Jold=[],Hold=[],idx=[];
for(let i=0;i<g.ls.length;i++)if(!g.deleted.has(i)){
 const x=g.ls[i],h=g.I.map(k=>x[k]);if(!h.some(Boolean))continue;
 Hold.push(h);Jold.push(g.J.map(k=>x[k]));idx.push(i);
}
const out=[],all=[];
for(const o of near.edge){
 const d=base.wit.find(x=>x.i===o.i), B=d.B, D=d.D;
 const S=new Set(o.found[0].blockers);
 if(S.size!==2||!o.found.every(x=>x.blockers.every(y=>S.has(y))))throw Error("owner pair");
 const owners=[...S].map(i=>g.ls[i]);
 const head=o.found[0].t.slice(0,8);
 const mid=norm(owners[0].map((x,j)=>x+owners[1][j]))===4*96;
 const aa=g.I.every((k,j)=>2*head[j]===owners[0][k]+owners[1][k])&&g.J.every((k,j)=>2*o.b[j]===owners[0][k]+owners[1][k]);
 const W=new Map(d.W0.map(w=>[w[0].join(),w]));
 for(let j=0;j<Hold.length;j++){
  if(S.has(idx[j]))continue;
  const a=Array.from({length:d.d},(_,k)=>rd(Hold[j],B.slice(0,8).map(x=>x[k]))),b=D*(64-rd(o.b,Jold[j]));
  if(!a.some(Boolean)||8*a.reduce((s,x)=>s+Math.abs(x),0)<=b)continue;
  const key=a.join(),p=W.get(key);
  if(!p||b<p[1])W.set(key,[a,b,"eq"+idx[j]]);
 }
 const P=verts([...W.values()],d.d);
 const Q=P.P.map(a=>{const v=B.map(x=>rd(x,a));return {a,v,n:norm(v),h:(D*a[d.d])**2};});
 const u=Q.reduce((p,q)=>q.n*p.h>p.n*q.h?q:p),ft=Q.filter(x=>x.n>=64*x.h);
 out.push({i:o.i,pair:[...S],site:o.b,owners,mid,aa,ownProducts:owners.map(x=>rd(head,g.I.map(i=>x[i]))+rd(o.b,g.J.map(i=>x[i]))),oldOwnerNormLabels:owners.map(x=>norm(g.I.map(j=>x[j]))),sn:[u.n,u.h],candidates:ft.length,vertices:P.P.length,B:d.B,D:d.D,W:[...W.values()]});
 for(const q of o.found){all.push({i:o.i,b:o.b,v:q.t,n:q.n,c:q.blockers});}
}
const bad=[],dup=new Set(),shared=[];
for(let i=0;i<all.length;i++)for(let j=0;j<i;j++){
 const x=all[i],y=all[j];const a=rd(x.b,y.b),v=rd(x.v,y.v);
 if(v>0&&4096*v*v>(64-a)**2*x.n*y.n)bad.push([i,j,"mutual"]);
 if(x.i!==y.i&&x.c.some(w=>y.c.includes(w)))shared.push([i,j]);
}
return {out,bad,shared,all};
};
const singles=(g,base)=>{
const gcd=(a,b)=>{a=Math.abs(a);b=Math.abs(b);while(b){let c=a%b;a=b;b=c;}return a;};
const dot=(a,v)=>a.reduce((s,x,i)=>s+x*v[i],0);
function canonical(a,d){let D=a.reduce(gcd,0);if(!D||a[d]<=0)throw Error("homogeneous");return a.map(x=>x/D);}
function verts(W,d){
 let P=Array.from({length:1<<d},(_,k)=>{
  const a=Array.from({length:d},(_,i)=>k>>i&1?8:-8).concat(1);
  let M=0n;for(let i=0;i<d;i++)M|=1n<<BigInt(2*i+(k>>i&1));
  return {a,M};
 });
 const ww=W.slice().sort((u,v)=>u[1]/Math.sqrt(dot(u[0],u[0]))-v[1]/Math.sqrt(dot(v[0],v[0])));
 let maximum=P.length,made=0;
 for(let s=0;s<ww.length;s++){
  const [b,r0] = ww[s], h=1n<<BigInt(2*d+s);
  const Q=P.map(({a})=>dot(b,a)-r0*a[d]);
  let low=[],high=[];
  for(let j=0;j<P.length;j++){if(Q[j]<0)low.push(j);if(Q[j]>0)high.push(j);}
  if(!high.length){for(let j=0;j<P.length;j++)if(!Q[j])P[j]={a:P[j].a,M:P[j].M|h};continue;}
  const B=new Map();
  for(let j=0;j<P.length;j++)if(Q[j]<=0){
   const w=P[j];B.set(w.a.join(),{a:w.a,M:w.M|(!Q[j]?h:0n)});
  }
  let C=new Map();
  for(let j of low)for(let k of high){
   const x=P[j], y=P[k], M=x.M&y.M;
   let c=C.get(M);
   if(c===undefined){c=0;for(const p of P)if((p.M&M)===M){c++;if(c>=3)break;}C.set(M,c);}
   if(c!==2)continue;
   const raw=x.a.map((a,i)=>Q[k]*a-Q[j]*y.a[i]);
   if(raw.some(a=>!Number.isSafeInteger(a)))throw Error("int precision");
   const a=canonical(raw,d),k0=a.join(),old=B.get(k0);
   B.set(k0,{a,M:(old?.M??0n)|M|h});made++;
  }
  P=[...B.values()];
  maximum=Math.max(maximum,P.length);
  if(!P.length)throw Error("origin always feasible");
 }
 for(const {a} of P)if(W.some(([q,t])=>dot(q,a)>t*a[d]))throw Error("invalid vertices");
 return {P:P.map(w=>w.a),maximum,made};
}


const dotn=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0), norm=x=>dotn(x,x);
const R=g.added.map(x=>g.I.map(i=>x[i]).concat(x.slice(24))), J=g.added.map(x=>g.J.map(i=>x[i]));
const EJ=[],EH=[],EL=[];
for(let j=0;j<g.ls.length;j++)if(!g.deleted.has(j)){
 const x=g.ls[j],h=g.I.map(i=>x[i]);if(!h.some(Boolean))continue;
 EH.push(h);EJ.push(g.J.map(i=>x[i]));EL.push(j);
}
const out=[],win=[],maxhist={};
function codeMax(C){
 const d=C.length, ad=C.map((p,j)=>C.map((q,i)=>i!==j&&dotn(p,q)<=0));
 let best=[];
 function rec(S,X){
  if(X.length>best.length)best=X;
  if(best.length>=5||X.length+S.length<=best.length)return;
  while(S.length){
   if(X.length+S.length<=best.length)return;
   const i=S.pop();rec(S.filter(j=>ad[i][j]),X.concat(i));
   if(best.length>=5)return;
  }
 }
 rec(Array.from({length:d},(_,i)=>i),[]);return best;
}
for(let ii=0;ii<g.S.length;ii++){
 const f=g.S[ii];if(f.orig.length!==20||f.ant!==16)continue;
 const known=new Set(f.orig.map(w=>w.join()));
 const U=f.orig.filter(w=>!known.has(w.map(x=>-x).join()));
 const rm=new Set(U.map(x=>x.join())), eq=f.l.filter(j=>rm.has(g.I.map(i=>g.ls[j][i]).concat([0,0,0,0,0]).join())), 
   off=[...f.mid,...f.ph].filter(j=>rm.has(R[j].join()));
 if(rm.size!==4||eq.length+off.length!==4)throw Error("four singles");
 const d0=base.wit.find(x=>x.i===ii), d=d0.d,B=d0.B,D=d0.D;
 if(d!==5||!d0)throw Error("kernel");
 const W=new Map(), Bs=Array.from({length:d},(_,k)=>B.map(a=>a[k]));
 const put=(r,j,owner)=>{
  const a=Bs.map(u=>dotn(r,u)),b=D*(64-dotn(f.b,j));
  if(8*a.reduce((s,x)=>s+Math.abs(x),0)<=b)return;
  const k=a.join(),w=W.get(k);if(!w||b<w[1])W.set(k,[a,b,owner]);
 };
 for(let j=0;j<R.length;j++)if(!off.includes(j))put(R[j],J[j],"o"+j);
 for(let j=0;j<EH.length;j++)if(!eq.includes(EL[j]))put(EH[j].concat([0,0,0,0,0]),EJ[j],"e"+EL[j]);
 const P=verts([...W.values()],d), 
     Q=P.P.map(a=>{const x=B.map(z=>dotn(z,a)), n=norm(x), e=D*a[d];return {a,x,n,e};});
 const max=Q.reduce((z,q)=>q.n*z.e*z.e>z.n*q.e*q.e?q:z), k=max.n/max.e/max.e;
 const C=new Map(), putC=x=>{const s=x.reduce(gcd,0),v=x.map(t=>t/s);C.set(v.join(),v);};
 for(let q of Q)if(q.n>=64*q.e*q.e)putC(q.x);
 const near=Q.filter(q=>q.n>=49*q.e*q.e);
 for(let i=0;i<near.length;i++)for(let j=0;j<i;j++){
  const p=near[i],q=near[j],s=p.x.map((x,k)=>x*q.e+q.x[k]*p.e),den=2*p.e*q.e;
  if(norm(s)>=64*den*den)putC(s);
 }
 const rays=[...C.values()],cl=codeMax(rays), info={ii,k,wall:W.size,vertices:P.P.length,rays:rays.length,cl:cl.length,eq,off,b:f.b,chosen:cl.map(i=>rays[i])};
 if(cl.length>4)win.push(info);
 out.push(info);
}
return {out,win};
};
const had=(g,mode,CAP)=>{
const R=g.added.map(x=>g.I.map(i=>x[i]).concat(x.slice(24)));
const J=g.added.map(x=>g.J.map(i=>x[i]));
const norm=x=>x.reduce((s,v)=>s+v*v,0),dot=(x,y)=>{let s=0;for(let j=0;j<x.length;j++)s+=x[j]*y[j];return s;};
const J2=J.map(norm), keys=J.map(j=>j.join()), out=[],stats={};
const comb=(n,k)=>{let a=[];function rec(s,P){if(P.length===k){a.push(P);return;}for(let j=s;j<n;j++)rec(j+1,P.concat(j));}rec(0,[]);return a;};
const tris=comb(12,3);
const PM=(x,y,z,t,p)=>{
 const w=Math.abs(x)+Math.abs(y)+Math.abs(z)+Math.abs(t);
 if(!x||!y||!z||!t||((x<0)+(y<0)+(z<0)+(t<0))%2===p)return w;
 return w-2*Math.min(Math.abs(x),Math.abs(y),Math.abs(z),Math.abs(t));
};
const axis=Array.from({length:10},()=>[]);
for(let j of g.aux){
 const z=R[j].slice(8);
 if(z.some(x=>Math.abs(x)===8)){
  const u=z.findIndex(x=>x!==0);
  axis[2*u+(z[u]>0?1:0)].push(J[j]);
 }
}
const frameInfo=[];
for(let ii=0;ii<g.S.length;ii++){
 const f=g.S[ii];if(f.orig.length!==24||f.ph.length!==(mode==0?8:0))continue;
 if(mode==0&&f.l.length!==16)continue;
 if(mode==1&&f.mid.length!==16)continue;
 const j=f.axis[0];
 const p=axis.slice(2*j,2*j+2).map(L=>Math.max(...L.map(x=>dot(x,f.b))));
 if(Math.max(...p)>32)continue;
 frameInfo.push(ii);
}
let passed0=0,total=0,scoped=0,full=0,wallCounts={};
for(const ii of frameInfo){
 const f=g.S[ii],j=f.axis[0];
 const ors=f.orig.filter(x=>x.find(y=>y)!==undefined&&x.find(y=>y)>0);
 const F=f.heads.filter(x=>x.find(a=>a)>0), Q=ors.map(x=>F.map(a=>dot(a,x.slice(0,8))));
 if(F.length!==8||ors.length!==12)throw Error("frames");
 const T=[];
 for(const [a,b,c] of tris){
  if(mode==0&& (norm(ors[a].slice(8))>0)+(norm(ors[b].slice(8))>0)+(norm(ors[c].slice(8))>0)!==1)continue;
  // all eight signed head sums over those three axes clear the entire Leech frame convex hull
  let ok=true;
  for(let sig of [[1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1]]){
   let l=0;for(let i=0;i<8;i++)l+=Math.abs(sig[0]*Q[a][i]+sig[1]*Q[b][i]+sig[2]*Q[c][i]);
   if(l>128){ok=false;break;}
  }
  if(!ok)continue;
  // all pole tests, including the newly active vacancy, are exact
  if([0,1,2,3,4].some(k=>(k!==j)&&[a,b,c].filter(r=>ors[r][8+k]!==0).length>2))continue;
  T.push([a,b,c]);
 }
 if(!T.length)continue;
 total+=T.length*2;scoped++;

 const key=f.b.join(),FQ=[];
 // Begin with strongly restricted J64 rows from nearby fibres, then handle short/zero heads.
 for(const h of g.S){
  if(h===f)continue;
  const px=dot(f.b,h.b);
  if(px<=0)continue;
  FQ.push([px,h]);
 }
 FQ.sort((a,b)=>b[0]-a[0]);
 let CC=T.flatMap((x)=>[x.concat(0,[[]]),x.concat(1,[[]])]),rowCount=0;
 const filter=(k,G)=>{
   if(!CC.length)return;const c=8*R[k][8+j];
   const A=ors.map(r=>dot(r,R[k])),m=Math.max(...A.map(Math.abs));
   if(3*m+Math.abs(c)<=2*G)return;
   rowCount++;
   for(let u=CC.length-1;u>=0;u--){
    const [a,b,c0,p,ee]=CC[u];
    if(PM(A[a],A[b],A[c0],c,p)>2*G){
      ee.push(k);if(ee.length>CAP){CC[u]=CC.at(-1);CC.pop();}
    }
   }
 };
 for(const [px,h] of FQ){
  if(!CC.length)break;
  const G=64-px;for(const k of h.mid)filter(k,G);
  for(const k of h.ph)filter(k,G);
 }
 if(CC.length){
  for(let k=0;k<J.length&&CC.length;k++){
   // Everything with standard J64 was already checked when it could be relevant.
   if(J2[k]===64&&g.Lx.has(keys[k]))continue;
   if(keys[k]===key)continue;
   const G=64-dot(f.b,J[k]);
   if(G*G>=64*(128-J2[k]))continue;
   filter(k,G);
  }
 }
 wallCounts[Math.floor(rowCount/10)*10]=(wallCounts[Math.floor(rowCount/10)*10]??0)+1;
 for(const [a,b,c,p,blockers] of CC){
  passed0++;
  const r=[ors[a],ors[b],ors[c]];
  let vv=[];
  for(let ss=0;ss<16;ss++){
   const sg=Array.from({length:4},(_,k)=>ss>>k&1?-1:1);
   if(sg.filter(x=>x<0).length%2!==p)continue;
   const v=Array.from({length:13},(_,k)=>(sg[0]*r[0][k]+sg[1]*r[1][k]+sg[2]*r[2][k]+(k===8+j?8*sg[3]:0))/2);
   if(norm(v)!==64||v.some(x=>!Number.isInteger(x)))throw Error("daughter");
   const pt=Array(29).fill(0);g.I.forEach((k,u)=>pt[k]=v[u]);g.J.forEach((k,u)=>pt[k]=f.b[u]);for(let k=0;k<5;k++)pt[24+k]=v[8+k];
   vv.push(pt);
  }
  out.push({ii,j,p,old:r,body:f.b,ids:[a,b,c],newrows:vv,blockers});
 }
}
return {candidateFrames:frameInfo.length,total,scoped,passed0,out,wallCounts,summ:mode};
};
const move=(g,near)=>{
const gcd=(a,b)=>{a=Math.abs(a);b=Math.abs(b);while(b){let c=a%b;a=b;b=c;}return a;};
const dot=(a,v)=>a.reduce((s,x,i)=>s+x*v[i],0);
function canonical(a,d){let D=a.reduce(gcd,0);if(!D||a[d]<=0)throw Error("homogeneous");return a.map(x=>x/D);}
function verts(W,d){
 let P=Array.from({length:1<<d},(_,k)=>{
  const a=Array.from({length:d},(_,i)=>k>>i&1?8:-8).concat(1);
  let M=0n;for(let i=0;i<d;i++)M|=1n<<BigInt(2*i+(k>>i&1));
  return {a,M};
 });
 const ww=W.slice().sort((u,v)=>u[1]/Math.sqrt(dot(u[0],u[0]))-v[1]/Math.sqrt(dot(v[0],v[0])));
 let maximum=P.length,made=0;
 for(let s=0;s<ww.length;s++){
  const [b,r0] = ww[s], h=1n<<BigInt(2*d+s);
  const Q=P.map(({a})=>dot(b,a)-r0*a[d]);
  let low=[],high=[];
  for(let j=0;j<P.length;j++){if(Q[j]<0)low.push(j);if(Q[j]>0)high.push(j);}
  if(!high.length){for(let j=0;j<P.length;j++)if(!Q[j])P[j]={a:P[j].a,M:P[j].M|h};continue;}
  const B=new Map();
  for(let j=0;j<P.length;j++)if(Q[j]<=0){
   const w=P[j];B.set(w.a.join(),{a:w.a,M:w.M|(!Q[j]?h:0n)});
  }
  let C=new Map();
  for(let j of low)for(let k of high){
   const x=P[j], y=P[k], M=x.M&y.M;
   let c=C.get(M);
   if(c===undefined){c=0;for(const p of P)if((p.M&M)===M){c++;if(c>=3)break;}C.set(M,c);}
   if(c!==2)continue;
   const raw=x.a.map((a,i)=>Q[k]*a-Q[j]*y.a[i]);
   if(raw.some(a=>!Number.isSafeInteger(a)))throw Error("int precision");
   const a=canonical(raw,d),k0=a.join(),old=B.get(k0);
   B.set(k0,{a,M:(old?.M??0n)|M|h});made++;
  }
  P=[...B.values()];
  maximum=Math.max(maximum,P.length);
  if(!P.length)throw Error("origin always feasible");
 }
 for(const {a} of P)if(W.some(([q,t])=>dot(q,a)>t*a[d]))throw Error("invalid vertices");
 return {P:P.map(w=>w.a),maximum,made};
}


const dotn=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0), GM=[64,32,32,32];
const norm=(x)=>x.slice(0,4).reduce((s,w,i)=>s+GM[i]*w*w,0);
const out=[],win=[];
function anchorTriangle(C){
 const d=C.length;let edges=Array.from({length:d},()=>[]);
 const comp=(a,b)=>{
  const u=a.reduce((s,x,i)=>s+GM[i]*x*b[i],0);
  return u<=0||4*u*u<=norm(a)*norm(b);
 };
 for(let k=0;k<d;k++)for(let j=0;j<k;j++)if(comp(C[k],C[j])){
  for(let i of edges[j])if(comp(C[k],C[i]))return [i,j,k];
  edges[k].push(j);
 }
 return [];
}
for(const o of near.edge){
 const q=o.found[0], pay=q.blockers,s=new Set(pay), X0=g.ls[pay[0]],X1=g.ls[pay[1]];
 const B=Array.from({length:4},()=>Array(29).fill(0));
 g.J.forEach((j,i)=>{B[0][j]=o.b[i];if(X0[j]!==o.b[i]||X1[j]!==o.b[i])throw Error("control J");});
 g.I.forEach(k=>{B[1][k]=(X0[k]+X1[k])/2; B[2][k]=(X0[k]-X1[k])/2;});
 for(let k=0;k<5;k++)B[3][24+k]=q.t[8+k];
 if(B.some((r,k)=>dotn(r,r)!==GM[k]||B.some((v,j)=>k!==j&&dotn(r,v))))throw Error("axes");
 const W=new Map();const aeq=(x)=>B.map(u=>dotn(x,u));
 function ins(a,owner){
  if(!a.some(Boolean))return;
  const k=a.join();
  if(!W.has(k))W.set(k,[a,64,owner]);
 }
 for(let i=0;i<g.ls.length;i++)if(!g.deleted.has(i)&&!s.has(i))ins(aeq(g.ls[i]),"e"+i);
 for(let i=0;i<g.added.length;i++)ins(aeq(g.added[i]),"o"+i);
 // Radius-128 sphere fits these ambient coordinate boxes; they do not cut any candidate.
 for(let i=0;i<4;i++)for(let p of [-1,1]){
   const w=Array.from({length:4},(_,k)=>i===k?p:0);
   W.set(w.join(),[w,2,"box"]);
 }
 const P=verts([...W.values()],4);
 const qn=P.P.slice().map(a=>[a,norm(a),a[4]**2]);
 const M=qn.reduce((a,b)=>b[1]*a[2]>a[1]*b[2]?b:a);
 const C=new Map();
 const put=x=>{const s=x.reduce(gcd,0),v=x.map(y=>y/s);C.set(v.join(),v);};
 for(const [a,n,d] of qn)if(n>=128*d)put(a.slice(0,4));
 const base=C.size, L=qn.filter(([a,n,d])=>n>=110*d);
 for(let i=0;i<L.length;i++)for(let j=0;j<i;j++){
  const a=L[i][0],b=L[j][0],c=a.slice(0,4).map((u,k)=>u*b[4]+b[k]*a[4]),m=(2*a[4]*b[4])**2;
  if(norm(c)>=128*m)put(c);
 }
 const ray=[...C.values()],t=anchorTriangle(ray);
 const info={i:o.i,body:o.b,three:t.map(j=>ray[j]),paid:pay,B,max:[M[1],M[2]],wall:W.size,vertices:P.P.length,base,rays:ray.length};
 if(t.length)win.push(info);out.push(info);
}
return {win,out};
};
const S=g(decode(bytes,x=>Uint8Array.from(inflateRawSync(x))));
const A=ax(S),diagCount=diag(S,A),P=maxpol(polish(S)),M=mix(S),E=fold(S,M),B=paid(S,M,E),U=singles(S,M);
if(A.pass.length!==335||P.rows.length!==64||P.wit.length!==0||M.rows.length!==32||M.wit.length!==23||E.win.length!==0||E.edge.length!==11||
B.bad.length||B.shared.length||B.out.some(o=>!o.aa||o.sn[0]!==64*o.sn[1]||o.candidates!==4)||
U.out.length!==11||U.win.length||U.out.some(o=>o.k!==64||o.cl!==4))throw Error("final proof count");
const moved=move(S,E);
if(moved.out.length!==11||moved.win.length||moved.out.some(x=>x.max[0]!==128*x.max[1]||x.base!==4||x.rays!==4))throw Error("moving-body polar maxima");
const summary={movingBodyPolars:moved.out.map(x=>[x.i,128,x.wall,x.vertices,x.rays]),headProfiles:S.hist,axialEnvelope:A.hist,onePairExactDiagonals:diagCount,pureTailPolar:P.hist,fullProjectedPolar:E.hist,
neutralPackets:B.out.map(o=>[o.i,o.pair]),allSingles:U.out.map(o=>[o.ii,o.k,o.rays,o.cl,o.eq,o.off])};
if(process.argv.includes("--hadamard")){
 const h=[had(S,0,4),had(S,1,4)];
 if(h[0].total!==46592||h[1].total!==453248||h.some(q=>q.out.length))throw Error("Hadamard counts");
 summary.hadamard={candidates:h.map(q=>q.total),everyPacketAtLeastFiveOffEquatorBlockers:true};
}
console.log(JSON.stringify(summary,null,2));
