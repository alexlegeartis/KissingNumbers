import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {resolve} from 'node:path';

const base=resolve(process.argv[2]??'.');
function input(path,sha,git=false){
 const buffer=readFileSync(resolve(base,path));
 const h=createHash(git?'sha1':'sha256');
 if(git)h.update('blob '+buffer.length+'\0');h.update(buffer);
 assert.equal(h.digest('hex'),sha,path);
 return JSON.parse(buffer.toString('utf8'));
}
const dep='verifications/improved/dim73-82-gamma72-depth6-lattice/';
const cert=input(dep+'certificate.json','e5e656d728ef6ca1a7242e2863c8c41447247a9165a7b5aa862ea15d11a1584e');
const gr=input(dep+'catalog_gram.json','1791fab8c152302db2652ec66f227eb473ce4007f934431ef6bbb2ab2f2d87ab');
const live=input('verifications/improved/dim71-gamma72-single-block-refill/pattern.json','605d91b1a19e189379f52254bb8144b2495e0f50',true);
assert.equal(live.prior+live.added-live.removed,2623656288);
const D=((cert,gr,lam)=>{

const N=24,I=Array.from({length:N},(_,i)=>Array.from({length:N},(_,j)=>i===j?1:0));
const tr=M=>M[0].map((_,i)=>M.map(r=>r[i]));
const mm=(M,T)=>M.map(r=>T[0].map((_,j)=>r.reduce((s,v,k)=>s+v*T[k][j],0)));

const mod2=M=>M.map(r=>r.map(x=>x&1));
function inv2(M){
 const Z=mod2(M).map((r,i)=>r.concat(I[i]));
 for(let j=0;j<N;j++){
  let p=j;while(p<N&&!Z[p][j])p++;if(p===N)throw Error("singular gf2");
  [Z[j],Z[p]]=[Z[p],Z[j]];
  for(let i=0;i<N;i++)if(i!==j&&Z[i][j])for(let k=j;k<2*N;k++)Z[i][k]^=Z[j][k];
 }
 return Z.map(r=>r.slice(N));
}
const equal=(A,B)=>JSON.stringify(A)===JSON.stringify(B);
const {leech_basis_change:U,leech_frame:J}=cert;
const G=Array.from({length:72},(_,i)=>gr.slice(i*72,(i+1)*72));
const F0=G.slice(0,24).map(r=>r.slice(0,24).map(x=>x/2));

const CX=G.slice(0,24).map(r=>r.slice(48,72));
const A0=mod2(mm(CX,inv2(F0)));
if(!equal(mod2(mm(A0,F0)),mod2(CX)))throw Error("A0 exact mod2");
const Ui=inv2(U);
if(!equal(mod2(mm(U,Ui)),I)||!equal(mod2(mm(Ui,U)),I))throw Error("inverse exact mod2");
const F=mm(mm(U,F0),tr(U)),A=mod2(mm(mm(U,A0),Ui));
const FI=mm(F,tr(J)),B=I.map((r,i)=>r.map((x,j)=>x-A[i][j]));
if(!equal(mod2(mm(A,A)),A))throw Error('idempotent A');
for(const P of [A,B]){
 const H=mm(mm(P,F),tr(P));
 if(H.some((r,i)=>r.some((x,j)=>x%(i===j?4:2))))throw Error('non-singular residue half');
}
if(!equal(mm(J,FI),I.map(r=>r.map(x=>8*x))))throw Error('frame exact');
const p24=r=>r.reduce((acc,v,i)=>acc|((v&1)<<i),0);
const pa=A.map(p24),pb=B.map(p24);
function lab(x,proj){let z=0;for(;x;x&=x-1)z^=proj[31-Math.clz32(x&-x)];return z;}
const pop=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0x0f0f0f0f)*0x01010101)>>>24;};
const Jt=tr(J);
function c24(q){let s=0;for(let j=0;j<N;j++){
  let d=0;const t=Jt[j];for(let i=0;i<N;i++)d+=q[i]*t[i];
  if(d&7)return -1;s|=(((d>>3)&1)<<j);
 }return s;}
const pack4=r=>p24(r)+0x1000000*p24(r.map(x=>x>>1));
function xor4(x,y){const a=x&0xffffff,b=y&0xffffff,c=Math.floor(x/0x1000000),d=Math.floor(y/0x1000000);
return (a^b)+0x1000000*(c^d^(a&b));}
let R=new Set([0]);
for(const r of FI){const a=pack4(r),aa=xor4(a,a),aaa=xor4(a,aa);const S=new Set();for(const x of R){S.add(x);S.add(xor4(x,a));S.add(xor4(x,aa));S.add(xor4(x,aaa));}R=S;}
if(R.size!==8192)throw Error("R");
const four=[];
const push=(q)=>{const t=c24(q);if(t>=0){
 let j=0;while(!q[j])j++;if(q[j]>0)four.push({q:Int8Array.from(q),a:lab(t,pa),b:lab(t,pb)});
}};
for(let i=0;i<N;i++)for(let j=i+1;j<N;j++)for(const x of [4,-4]){const q=new Int8Array(N);q[i]=4;q[j]=x;push(q);}
for(const a of R){
 const lo=a&0xffffff,hi=Math.floor(a/0x1000000);
 if(lo===0&&pop(hi)===8){
  const idx=Array.from({length:N},(_,j)=>j).filter(j=>(hi>>>j)&1),q=new Int8Array(N);
  for(let u=0;u<128;u++){
   idx.forEach((j,k)=>q[j]=((u<<(1))>>>k)&1?-2:2); push(q);
  }
 } else if(lo===0xffffff){
  const q=Int8Array.from({length:N},(_,j)=>(hi>>>j)&1?-1:1);
  for(let j=0;j<N;j++){q[j]*=-3;push(q);q[j]/=-3;}
 }
}
const weights=(u,v)=>{let d=0;for(let i=0;i<N;i++)d+=u[i]*v[i];return d;};
const qlam=mm([lam],FI)[0],ca=lab(p24(lam),pa),cb=lab(p24(lam),pb);
if(four.length!==98280||weights(qlam,qlam)!==32||!ca)throw Error([four.length,weights(qlam,qlam),ca]);
const g4=new Map();for(const z of four){if(!g4.has(z.a))g4.set(z.a,[]);g4.get(z.a).push(z);}
if(g4.size!==4095||g4.has(0)||[...g4.values()].some(L=>L.length!==24||new Set(L.map(z=>z.b)).size!==24))throw Error('norm4 group');
return {I,J,A,B,F,FI,p24,pa,pb,lab,pop,c24,R,four,g4,weights,qlam,ca,cb};

})(cert,gr,live.lambda);
const E=((D)=>{

const {four,g4,weights,qlam,ca,cb,lab,pa,pb,p24,pop}=D;
const K=new Map();
function qkey(q){return String.fromCharCode(...q.map(x=>x+20));}
function canon(x,y,t){
 const q=Array.from(x,(a,i)=>a+t*y[i]);if (q.find(x=>x!==0)<0) for(let i=0;i<q.length;i++)q[i]=-q[i];return q;
}
for(const L of g4.values())for(let i=0;i<24;i++)for(let j=i+1;j<24;j++){
 const u=L[i],v=L[j],b=u.b^v.b;
 if(!K.has(b))K.set(b,new Map());
 const m=K.get(b);if(m.size===24)continue;
 if(weights(u.q,v.q)!==0)throw Error('A4 frame');
 for(const t of [-1,1]){
  const q=canon(u.q,v.q,t),s=qkey(q);if(!m.has(s))m.set(s,Int8Array.from(q));
 }
 if(m.size>24)throw Error('oversize');
}
if(K.size!==4095||[...K.values()].some(X=>X.size!==24))throw Error("B8 frame "+[...K.values()].reduce((H,X)=>(H[X.size]=(H[X.size]||0)+1,H),{}));
const frame8=new Map([...K].map(([b,W])=>[b,[...W.values()]]));
const survive=[];
const census={};
for(const [b,W] of frame8){
 let z=[0,0,0,0,0];
 for(const q of W){
  const a=Math.abs(weights(qlam,q))/8;if(!Number.isInteger(a)||a>4)throw Error('dot');
  z[a]++;
 }
 if(z.some((_,i)=>z[i]&&z[(i+1)%2]))throw Error('parity');
 census[z]=(census[z]||0)+1;if(!z[3]&&!z[4])survive.push(b);
}
const W6a=new Map();
for(const [a,L] of g4)if(a<(a^ca)&&g4.has(a^ca)){
 const M=g4.get(a^ca);
 for(const u of L) for(const v of M){
  const b=u.b^v.b,d=weights(u.q,v.q);
  if (d!==8&&d!==-8 || W6a.has(b))continue;
  const q=canon(u.q,v.q,d>0?-1:1);
  if(weights(q,q)!==48)throw Error('6 norm');
  W6a.set(b,Int8Array.from(q));
 }
}
if(W6a.size!==2048)throw Error('W6a '+W6a.size);
for(const b of W6a.keys())if(!W6a.has(b^cb))throw Error("W6a parity");
return {frame8,survive,census,W6a,canon,qkey};

})(D);
E.passSet=new Set(E.survive);
assert.deepEqual(E.census,{'0,23,0,1,0':2048,'16,0,8,0,0':2024,'22,0,0,0,2':23});
assert.equal(E.frame8.size*48**2,9434880);
assert.equal(E.passSet.size,2024);
console.log('B-only raw stock 9434880; whole banks hit by (8,±3)/(8,4): 2071.');

const Q=((D,E)=>{

const {weights,cb,ca,pop}=D,{frame8,W6a}=E;
const axes=g4=>g4.get(ca).map(z=>z.q);
const L=axes(D.g4);
const w6m=new Map([...W6a].map(([b,w])=>{
  let m=0; for(let i=0;i<24;i++){
    const h=weights(L[i],w);if(h!==8&&h!==-8)throw Error('all raw norm-six phase magnitudes one');
    if(h<0)m|=1<<i;
  }
  return [b,m];
}));
const pairs=[...w6m].filter(([b])=>b<(b^cb)).map(([b,r])=>[r,w6m.get(b^cb)]);
if(pairs.length!==1024)throw Error('pair6');
const phase=new Map();
for(const [b,U] of frame8){
 const code=[];
 for(const u of U){
   let m=0,s=0,tri=-1,c2=0,c4=0;
   for(let i=0;i<24;i++){
     const h=weights(L[i],u)/8,j=1<<i;
     if(h<0)m|=j;if(h)s|=j;
     const a=Math.abs(h);
     if(a===3){if(tri!==-1)throw Error('two threes');tri=i;}
     else if(a===2)c2++;
     else if(a===4)c4++;
     else if(a!==0&&a!==1)throw Error('basis magnitude');
   }
   if(tri>=0&&(s!==0xffffff||c2||c4)||tri<0&&!(c2===8&&c4===0&&pop(s)===8||c2===0&&c4===2&&pop(s)===2))throw Error('basis 24-8-2 shapes');
   code.push({m,s,tri,c4});
 }
 phase.set(b,code);
}
function h(q,w){
   const d=(q.m^w)&q.s,n=pop(d);
   return q.tri>=0?Math.abs((26-2*n-4*(d>>>q.tri&1))/4):q.c4?Math.abs(2-2*n):Math.abs(4-n);
}
// Integer identity above is also spot-checked against a physical dot and no assumption about Golay equivalence is used.
const sample=[...W6a].slice(0,32);
for(const [b,U] of [...frame8].slice(0,20))for(let i=0;i<U.length;i++)for(const [c,r] of sample)if(h(phase.get(b)[i],w6m.get(c))!==Math.abs(weights(U[i],r))/8)throw Error('compressed norm-six dot');
return {phase,pairs,h,w6m};

})(D,E);
const R=((D,E,Q)=>{

const {weights,pop,qlam}=D,{h,pairs,phase}=Q;
let min6=Infinity,min8=Infinity,mn=Infinity,mx=0;
const H={},res=[],cats={};
const cap=0xffffff;
const stage2={histo:{},rows:[]};
const rr=new Int32Array(6),ss=new Int32Array(6);
for(const [b,U] of E.frame8){
 let e6=0,e8=0;
 const Uq=phase.get(b),odd=Uq[0].tri>=0;
 if(Uq.some(x=>(x.tri>=0)!==odd))throw Error('B parity');
 const active=new Uint32Array(24);active.fill(cap);
 let liveMask=E.passSet.has(b)?cap:0;
 for(let k=0;k<pairs.length;k++){
  let a2=0,a4=0,b2=0,b4=0;
  const [r,s]=pairs[k];rr.fill(0);ss.fill(0);
  for(let i=0;i<24;i++){
   const x=h(Uq[i],r),y=h(Uq[i],s);
   if((x&1)!==+odd||(y&1)!==+odd||x>5||y>5||!Number.isInteger(x)||!Number.isInteger(y))throw Error('fast local range');
   rr[x]++;ss[y]++;
   if(liveMask){
    if(x===2)a2|=1<<i;else if(x===4)a4|=1<<i;
    if(y===2)b2|=1<<i;else if(y===4)b4|=1<<i;
   }
  }
  if(odd){e6+=rr[1]*ss[5]+rr[5]*ss[1]+rr[3]*ss[3];e8+=rr[3]*ss[5]+rr[5]*ss[3];}
  else{e6+=rr[2]*ss[4]+rr[4]*ss[2];e8+=rr[4]*ss[4];}
  if(liveMask){
    liveMask=0;
    for(let i=0;i<24;i++){
      if(!active[i])continue;const m=1<<i;
      if(a4&m)active[i]&=~(b2|b4);
      if(a2&m)active[i]&=~b4;
      if(b4&m)active[i]&=~(a2|a4);
      if(b2&m)active[i]&=~a4;
      liveMask|=active[i];
    }
  }
 }
 const k6=e6*8,k8=e8*8,k=k6+k8;
 min6=Math.min(min6,k6);min8=Math.min(min8,k8);mn=Math.min(mn,k);mx=Math.max(mx,k);H[k]=(H[k]||0)+1;
 const key=odd?'odd':'even';
 cats[key]??={min6:Infinity,max6:0,min8:Infinity,max8:0,minSum:Infinity,maxSum:0};
 const C=cats[key];C.min6=Math.min(C.min6,k6);C.max6=Math.max(C.max6,k6);C.min8=Math.min(C.min8,k8);C.max8=Math.max(C.max8,k8);C.minSum=Math.min(C.minSum,k);C.maxSum=Math.max(C.maxSum,k);
 res.push([b,k6,k8]);
 if(E.passSet.has(b)){
   const n=active.reduce((s,m)=>s+pop(m),0),H2=stage2.histo;
   H2[n]=(H2[n]||0)+1;if(n)stage2.rows.push({b,active:[...active],ct:n});
 }
}
return {min6,min8,mn,mx,H,res,cats,stage2};

})(D,E,Q);
assert.deepEqual(R.stage2.histo,{0:1789,4:220,16:15});
assert.equal(R.stage2.rows.reduce((s,r)=>s+4*r.ct,0),4480);
assert.equal(R.min6,49152);
assert.equal(R.mn,52224);
assert.equal(R.mx,76608);
assert.deepEqual(R.cats.odd,{min6:76032,max6:76032,min8:576,max8:576,minSum:76608,maxSum:76608});
assert(22*2304<R.mn && 21*2304<R.min6 && 33*2304<R.cats.odd.minSum);
console.log('6+6 masks leave 4480 signed rays in 235 banks after the first cut.');
console.log('Trade: min 52224 distinct old blockers/bank; min 49152 norm-ten alone; every odd bank 76608.');

const C=((D,E,old)=>{

const {four,g4,weights,qlam,ca,cb,pop}=D,{frame8,canon,qkey}=E;
const basis=g4.get(ca),slots=new Map(basis.map(z=>[cb^z.b,z.q]));
if(slots.size!==24)throw Error("slots8");
const G=new Map(Array.from(slots.keys(),b=>[b,new Map()]));
for(const [a,L] of g4) if(a<(a^ca)&&g4.has(a^ca)){
 const M=g4.get(a^ca);
 for(const u of L)for(const v of M){
  const b=u.b^v.b,m=G.get(b);if(m===undefined||m.size===24)continue;
  if(weights(u.q,v.q)!==0)continue;
  for(const t of [-1,1]){const q=canon(u.q,v.q,t);m.set(qkey(q),Int8Array.from(q));}
  if(m.size>24)throw Error('oversize');
 }
}
if([...G.values()].some(x=>x.size!==24))throw Error('A+B frame '+[...G.values()].map(x=>x.size));
const contacts=[...G].flatMap(([b,W])=>[...W.values()].map(z=>[slots.get(b),z]));
if(contacts.length!==576)throw Error("contacts");

const hist={},res=[];
for(const row of old.rows){
  const U=frame8.get(row.b);let n=0;
  for(let i=0;i<24;i++)for(let j=0;j<24;j++)if(row.active[i]>>>j&1){
    let seen=false;
    for(const [r,s] of contacts){
      const q=weights(U[i],r)/8,t=weights(U[j],s)/8;
      const qr=weights(U[j],r)/8,tr=weights(U[i],s)/8;
      if(![q,t,qr,tr].every(Number.isInteger))throw Error("local integral range");
      const v=Math.abs(q)+Math.abs(t),vr=Math.abs(qr)+Math.abs(tr);
      if(v===6||v===8||vr===6||vr===8){seen=true;break;}
    }
    if(!seen){n++;res.push([row.b,i,j]);}
  }
  hist[n]=(hist[n]||0)+1;
}
return {hist,res,Bs:G.size};

})(D,E,R.stage2);
assert.deepEqual(C.hist,{0:235});
assert.equal(C.res.length,0);
assert.equal(C.Bs*24,576);
console.log('4+8 and 8+4 masks exclude all 4480 residual rays; no gain at fixed stock.');
