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

const {h,pairs,phase}=Q;
const masks=[0,0,0,0,0,0],nasks=masks.slice();
const stage2={histo:{},rows:[]},histDegree={};
let minimum=10000,minEven=10000,minOdd=10000;
const kMax=4095;
const bit=x=>31-Math.clz32(x&-x);
const mm=(m,h)=>h===1?m[5]:h===2?m[4]:h===3?m[3]|m[5]:h===4?m[2]|m[4]:h===5?m[1]|m[3]:0;
for(const [b,Uq] of phase){
 const odd=Uq[0].tri>=0;
 const degrees=new Uint16Array(24*24);
 const ha=new Uint8Array(24),hb=new Uint8Array(24);
 for(const [r,s] of pairs){
  masks.fill(0);nasks.fill(0);
  for(let i=0;i<24;i++){
    const a=h(Uq[i],r),b=h(Uq[i],s);
    ha[i]=a;hb[i]=b;masks[a]|=1<<i;nasks[b]|=1<<i;
  }
  for(let i=0;i<24;i++){
    const v=i*24;
    for(let j=mm(nasks,ha[i]);j;j&=j-1) degrees[v+bit(j)]++;
    for(let j=mm(masks,hb[i]);j;j&=j-1) degrees[v+bit(j)]++;
  }
 }
 const active=new Uint32Array(24),mh=odd?'odd':'even',hh=stage2.histo[mh]??={};
 let n=0;
 for(let i=0;i<24;i++)for(let j=0;j<24;j++){
    const d=degrees[24*i+j];
    if(!d){active[i]|=1<<j;n++;}
    else{minimum=Math.min(minimum,d);if(odd)minOdd=Math.min(minOdd,d);else minEven=Math.min(minEven,d);}
    histDegree[d]=(histDegree[d]||0)+1;
 }
 hh[n]=(hh[n]||0)+1;if(n)stage2.rows.push({b,active:[...active],ct:n,odd});
}
return {minimum,minEven,minOdd,stage2,histDegree};

})(D,E,Q);
assert.equal(R.minEven,4);
assert.equal(R.minOdd,13);
assert.deepEqual(R.stage2.histo,{odd:{0:2048},even:{0:1789,4:240,16:18}});
assert.equal(R.stage2.rows.length,258);
assert.equal(R.stage2.rows.reduce((s,r)=>s+r.ct,0),1248);
assert(R.stage2.rows.every(r=>!r.odd));
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


function canonRaw(U,V,sgn){
  const a=Array.from(U,(x,i)=>x+sgn*V[i]);
  const s=a.find(x=>x!==0);
  if(s===undefined)throw Error('unoccupied old tail');
  return s>0?a:a.map(x=>-x);
}
const hist={},summary={};
let rightMax=0; let processed=0,full=0,need=0,minDegree=Infinity,minUnion=Infinity,maxUnion=0,firstFailure=null,allOrbitEdges=0;
for(const row of old.rows){
  const U=frame8.get(row.b),rows=[];
  for(let i=0;i<24;i++)for(let j=0;j<24;j++)if(row.active[i]>>>j&1){
    const candidates=new Set();
    for(const [r,s] of contacts)for(const [a,b] of [[r,s],[s,r]]){
      const ar=weights(U[i],a)/8,bs=weights(U[j],b)/8;
      if(!Number.isInteger(ar)||!Number.isInteger(bs))throw Error("local integral range");
      const n=Math.abs(ar)+Math.abs(bs);
      if(n!==6&&n!==8)continue;
      for(const sigma of ar>0?[-1]:ar<0?[1]:[-1,1])for(const tau of bs>0?[-1]:bs<0?[1]:[-1,1]){
        const Y=canonRaw(U[i],a,sigma),Z=canonRaw(U[j],b,tau);
        if([weights(Y,Y)/8,weights(Z,Z)/8].includes(6))throw Error('C edge collides with B signature');
        if((weights(Y,Y)+weights(Z,Z))/8!==28-2*n)throw Error('raw z norms');
        candidates.add(qkey(Y)+qkey(Z));
      }
    }
    minDegree=Math.min(minDegree,candidates.size);allOrbitEdges+=candidates.size;
    rows.push([...candidates]);
  }
  const seen=new Map();
  function augment(i,visited){
    for(const edge of rows[i]){
      if(visited.has(edge))continue;
      visited.add(edge);
      if(!seen.has(edge)||augment(seen.get(edge),visited)){
        seen.set(edge,i); return true;
      }
    }
    return false;
  }
  const N=rows.filter((_,i)=>augment(i,new Set())).length;
  const model=new Map();for(const E of rows)for(const v of E)model.set(v,(model.get(v)||0)+1);
  const total=model.size;for(const t of model.values())rightMax=Math.max(rightMax,t);
  minUnion=Math.min(minUnion,total);maxUnion=Math.max(maxUnion,total);
  processed++;need+=rows.length;full+=N;
  hist[N+','+rows.length]=(hist[N+','+rows.length]||0)+1;
  if(N!==rows.length&&firstFailure===null)firstFailure={bank:row.b,N,rows:rows.length,total};
}
return {processed,need,full,hist,rightMax,minDegree,minUnion,maxUnion,firstFailure,allOrbitEdges};

})(D,E,R.stage2);
assert.equal(C.processed,258);
assert.equal(C.need,1248);
assert.equal(C.full,1248);
assert.equal(C.minDegree,40);
assert.equal(C.rightMax,2);
assert.equal(C.allOrbitEdges,62016);
assert.equal(C.firstFailure,null);
console.log('Hall gate: six/six min positive degrees even 4, odd 13; odd residue none.');
console.log('Even 4+8: 1248 left axes in 258 banks, degree >= 40, per-bank right degree <= 2.');
