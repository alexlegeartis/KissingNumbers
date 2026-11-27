const pc=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return Math.imul((x+(x>>>4))&0xf0f0f0f,0x01010101)>>>24};

const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const tr=A=>A[0].map((_,i)=>A.map(r=>r[i]));
const mm=(A,B)=>{const T=tr(B);return A.map(r=>T.map(t=>dot(r,t)));};
const I=Array.from({length:24},(_,i)=>Array.from({length:24},(_,j)=>+(i===j)));
const p24=(v)=>v.reduce((s,x,i)=>s+((x&1)<<i),0);
const mp=(x,M)=>{let y=0;for(let a=x;a;a&=a-1)y^=M[31-Math.clz32(a&-a)];return y;};
function invert(A){
 const V=A.map(p24), W=I.map(p24);
 for(let i=0;i<24;i++){
  const j=V.findIndex((v,k)=>k>=i&&(v>>i&1));if(j<0)throw Error('singular');
  [V[i],V[j]]=[V[j],V[i]];[W[i],W[j]]=[W[j],W[i]];
  for(let k=0;k<24;k++)if(k!==i&&(V[k]>>i&1)){V[k]^=V[i];W[k]^=W[i];}
 }
 return W;
}
function gf(A,B){const R=B.map(p24);return A.map(a=>Array.from({length:24},(_,j)=>mp(p24(a),R)>>j&1));}
function gamma(c,vals){
 const G=Array.from({length:72},(_,i)=>vals.slice(72*i,72*i+72));
 const U=c.leech_basis_change,J=c.leech_frame, F0=G.slice(0,24).map(r=>r.slice(0,24).map(x=>x/2));
 const X=G.slice(0,24).map(r=>r.slice(48,72));
 const if0=invert(F0),iu=invert(U);
 const A=U.map(r=>mp(mp(mp(p24(r),X.map(p24)),if0),iu));
 const B=A.map((x,i)=>x^(1<<i));
 const F=mm(mm(U,F0),tr(U)),FM=F.map(p24), fwd=mm(F,tr(J));
 if(mm(fwd,J).some((row,i)=>row.some((x,j)=>x!==8*+(i===j))))throw Error('frame');
 const add4=(a,b)=>[(a[0]^b[0]),(a[1]^b[1]^(a[0]&b[0]))];
 let R=[[0,0]];
 for(const r of fwd){
  const a=[p24(r),p24(r.map(x=>x>>1))],a2=add4(a,a);
  const T=new Map();for(const t of R)for(const b of [[0,0],a,a2,add4(a2,a)]){const p=add4(t,b);T.set(p.join(','),p);}
  R=[...T.values()];
 }
 if(R.length!==8192)throw Error('mod 4');
 const M=new Map(),Jcols=tr(J);
 let n=0;
 function reg(w){
  const u=Jcols.map(y=>dot(y,w));
  if(u.some(x=>x%8))return;
  n++;const x=p24(u.map(x=>x/8)), a=mp(x,A), b=mp(x,B);
  if(!M.has(a))M.set(a,new Set());M.get(a).add(b);
 }
 for(let i=0;i<24;i++)for(let j=i+1;j<24;j++)for(const a of [-4,4])for(const b of [-4,4]){
  const w=Array(24).fill(0);w[i]=a;w[j]=b;reg(w);
 }
 for(const [lo,hi] of R){
  if(!lo&&pc(hi)===8){
   const X=[];for(let m=hi;m;m&=m-1)X.push(31-Math.clz32(m&-m));
   for(let z=0;z<256;z++){const w=Array(24).fill(0);X.forEach((j,i)=>w[j]=z>>i&1?-2:2);reg(w);}
  }else if(lo===0xffffff){
   const v=Array.from({length:24},(_,j)=>hi>>j&1?-1:1);
   for(let j=0;j<24;j++){const w=v.slice();w[j]*=-3;reg(w);}
  }
 }
 if(n!==196560||M.size!==4095||M.has(0)||[...M.values()].some(s=>s.size!==24))throw Error('M4');
 return {A,B,FM,M4:[...M].sort((a,b)=>a[0]-b[0]).map(([a,V])=>[a,[...V]])};
}
function planes(d){

function linear(words){const B=new Uint32Array(24);for(let w of words){while(w){const i=31-Math.clz32(w);if(B[i])w^=B[i];else{B[i]=w;break;}}}return Array.from(B).filter(Boolean);}
function group(B){let L=[0];for(const i of B)L.push(...L.map(x=>x^i));return L;}
const MA=new Map(d.M4), A=group(linear(d.A)).sort((a,b)=>a-b);
const ord=new Map(A.map((a,i)=>[a,i]));
const D=new Map(d.M4.map(([a,L])=>[a,L[0]]));const K=new Map(A.map(a=>[a,mp(a,d.FM)]));
const E= new Map(), h={},MAX=[];
for(const [u,L] of d.M4){
 const B=linear(L.map(x=>x^L[0]));
 h[B.length]=(h[B.length]||0)+1;
 const p=A.filter(v=>{let w=K.get(v);return B.every(x=>pc(x&w)%2===0);});
 E.set(u,p);
}
const L=[],L0=[];const deg=Array(A.length).fill(0),Llow=[],hu={};
for(const [a,Ea] of E)for(const b of Ea){
 if(b<=a||(a^b)<=b)continue;
 const c=a^b,X=[a,b,c];
 if(!E.get(b).includes(a)||!E.get(b).includes(c)||!E.get(c).includes(a)||!E.get(c).includes(b))continue;
 if(X.some((v,i)=>pc(K.get(v)&D.get(X[(i+1)%3]))%2!==pc(K.get(v)&D.get(X[(i+2)%3]))%2))continue;
 let t=0;X.forEach((v,i)=>t|=(pc(K.get(v)&D.get(X[(i+1)%3]))%2)<<i);
 if(((t^(t>>1)^(t>>2))&1)!==0)continue; // t is linear on plane: a+b=c.
 L.push([a,b,c,t]);for(const x of X)deg[ord.get(x)]++;
 if(t)continue;
 if(X.some((x,i)=>MA.get(x).some(u=>MA.get(X[(i+1)%3]).includes(u))))continue;
 L0.push(X);
}
for(let i=0;i<A.length;i++)hu[deg[i]]=(hu[deg[i]]||0)+1;
return {h,L,L0,A,E:[...E],ord:[...ord],depthDegree:hu};
}
function invReal(A){
 const n=A.length, T=A.map((a,i)=>[...a,...Array.from({length:n},(_,j)=>+(i===j))]);
 for(let i=0;i<n;i++){
  let k=i;for(let j=i+1;j<n;j++)if(Math.abs(T[j][i])>Math.abs(T[k][i]))k=j;
  if(Math.abs(T[k][i])<1e-9)throw Error('singular');
  [T[i],T[k]]=[T[k],T[i]];const p=T[i][i];for(let j=i;j<2*n;j++)T[i][j]/=p;
  for(let k=0;k<n;k++)if(k!==i){const r=T[k][i];for(let j=i;j<2*n;j++)T[k][j]-=r*T[i][j];}
 }
 return T.map(r=>r.slice(n));
}
function hmat(C,vals){
 const G=Array.from({length:72},(_,i)=>vals.slice(72*i,72*i+72)),U=C.leech_basis_change,J=C.leech_frame,F0=G.slice(0,24).map(r=>r.slice(0,24).map(x=>x/2));
 const X=G.slice(0,24).map(r=>r.slice(48,72));
 const A0=mm(X,invReal(F0)).map(r=>r.map(Math.round)),Ui=invReal(U).map(r=>r.map(Math.round));
 const eq=(a,b)=>a.every((r,i)=>r.every((x,j)=>x===b[i][j]));
 if(!eq(mm(U,Ui),I)||!eq(mm(A0,F0),X))throw Error('A0');
 const A=mm(mm(U,A0),Ui),F=mm(mm(U,F0),tr(U));
 if(!eq(mm(mm(A,F),tr(A)),F.map(r=>r.map(x=>2*x))))throw Error('norm');
 const K=mm(mm(J,A),mm(F,tr(J))).map(r=>r.map(x=>x/8));
 const K2=K.map(r=>r.map(x=>x*2));
 return {A,K,K2,F,J};
}
function phaseGraph(P){

const decode=new Map(P.A.map((a,i)=>[a,i])), planes=P.L0.map(t=>t.map(a=>decode.get(a)));
const Good=new Uint8Array(4096*4096), V=[], global=[];
for(const [u,v,w] of planes){
 for(const a of [u,v,w])for(const b of [u,v,w])if(a!==b)Good[(a<<12)|b]=1;
 const start=V.length;
 for(const a of [0,u,v,w])for(const b of [0,u,v,w])for(const c of [0,u,v,w]){
  const set=new Set([a,b,c]);set.delete(0);if(set.size>=2)V.push([a,b,c]);
 }
 global.push([start,V.length]);
}
function is(a,b,c){
 if(a!==0&&b!==0&&a!==b)return (c===0||c===a||c===b||c===(a^b))&&Good[(a<<12)|b]!==0;
 const e=a||b;return !!Good[(e<<12)|c];
}
function test(i,j){
 const u=V[i],v=V[j];return is(u[0]^v[0],u[1]^v[1],u[2]^v[2]);
}
function ngh(i){return V.flatMap((u,j)=>test(i,j)?[j]:[]);}
function maximal(i,goal=8){
 const Q=ngh(i);
 const E=Q.map((u,p)=>Q.map((v,q)=>test(u,v)?q:-1).filter(x=>x>=0)).map(l=>new Set(l));
 const A=Array.from({length:Q.length},(_,i)=>i).sort((i,j)=>E[j].size-E[i].size);
 let best=[];
 function colors(P){
   const O=[],C=[];let col=0;
   while(P.length){col++;const L=[],remain=[];
    for(const x of P){if(L.every(y=>!E[x].has(y))){L.push(x);O.push(x);C.push(col);}else remain.push(x);}
    P=remain;
   }
   return [O,C];
 }
 function rec(P,X){
   if(X.length>best.length){best=X.slice();if(best.length+1>=goal)return true;}
   const [O,C]=colors(P);
   for(let i=O.length-1;i>=0;i--){
    if(X.length+C[i]<=best.length)break;
    const v=O[i];X.push(v);
    if(rec(O.slice(0,i).filter(x=>E[v].has(x)),X))return true;
    X.pop();
   }
   return false;
 }
 rec(A,[]);
 return {N:Q.length,best:[i,...best.map(q=>Q[q])].map(j=>V[j]),E:E.reduce((s,A)=>s+A.size,0)/2};
}
return {is,V,Good,global,maximal};
}
const assert=require('node:assert/strict'),fs=require('node:fs'),crypto=require('node:crypto');
const cbuf=fs.readFileSync(process.argv[2]),gbuf=fs.readFileSync(process.argv[3]);
const blobHash=x=>crypto.createHash('sha1').update(Buffer.from('blob '+x.length+'\0')).update(x).digest('hex');
assert.equal(blobHash(cbuf),'5cc5386748aa477917fa7fd113a9c54a343800f5');
assert.equal(crypto.createHash('sha256').update(gbuf).digest('hex'),'1791fab8c152302db2652ec66f227eb473ce4007f934431ef6bbb2ab2f2d87ab');


const C=JSON.parse(cbuf),G=JSON.parse(gbuf),d=gamma(C,G),P=planes(d),H=hmat(C,G);
assert.equal(P.L0.length,975);
for(let i=0;i<P.A.length;i++)for(let b=0;b<12;b++)assert.equal(P.A[i]^P.A[1<<b],P.A[i^(1<<b)]);
const k=phaseGraph(P);assert.equal(k.V.length,40950);
for(const [idx,last] of k.global){
 assert.equal(last-idx,42);
 const s=k.maximal(idx,8);assert.equal(s.N,110);assert.equal(s.best.length,7);
}
const U=720,V=8359937;
const triples=Array.from({length:8},(_,i)=>{const x=i&1,y=(i>>1)&1,z=i>>2;return [(x?U:0)^(z?V:0),(y?U:0)^((x^z)?V:0),(z?U:0)^(y?V:0)];});
assert(P.L0.some(S=>S.includes(U)&&S.includes(V)&&S.includes(U^V)));
const gen=[
 x=>[x[0],x[1],x[2]^x[1]],
 x=>[x[0],x[0]^x[1]^x[2],x[2]],
 x=>[x[0]^x[1],x[1],x[2]],
 x=>[x[0],x[1]^x[0],x[2]],
 x=>[x[0],x[1],x[2]^x[0]]
];
const Kt=tr(H.K);assert(H.K.flat().every(x=>Number.isInteger(4*x)));
function iso(s,path){
 let x=s.slice(0,24),y=s.slice(24,48),z=s.slice(48);
 for(const i of path){
  if(i===0)[x,y]=[y,x];else if(i===1)[y,z]=[z,y];
  else if(i===3)x=x.map(v=>-v);else if(i===4)z=z.map(v=>-v);
  else if(i===2){
   const betaz=Kt.map((f,j)=>z[j]-dot(z,f));
   const p=x.map((w,j)=>(w+y[j]+betaz[j])/2);
   const ap=Kt.map(f=>dot(p,f));
   x=x.map((w,j)=>w-p[j]);y=y.map((w,j)=>w-p[j]);z=z.map((w,j)=>w-ap[j]);
  }else throw Error('gen');
 }
 const r=[...x,...y,...z];assert(r.every(Number.isSafeInteger));assert.equal(dot(r,r),384);
 return r;
}
function integralPhase(x){
 const gs=[];for(let j=0;j<3;j++){
  const w=x.slice(24*j,24*j+24);assert.equal(dot(w,w),128);
  const q=tr(H.J).map(t=>dot(t,w)/16);assert(q.every(Number.isSafeInteger));
  const z=p24(q);gs.push([mp(z,d.A),mp(z,d.B)]);
 }
 assert.deepEqual([gs[0][0]^gs[2][0],gs[1][0]^gs[2][0],gs[0][1]^gs[1][1]^gs[2][1]],[U^V,V,0]);
}
const X=[-2,-2,-2,-2,-2,-2,-2,-2,2,-2,-2,-2,2,-2,6,-2,-2,2,-2,-2,-2,-2,-2,-2,4,0,4,0,0,-4,0,0,4,0,0,0,0,0,-4,0,0,0,0,-4,-4,-4,0,0,-4,0,0,-4,0,0,0,0,-4,0,4,0,-4,0,0,0,0,-4,4,0,0,0,4,0];
integralPhase(X);
const literal=[[4,[0,0,0,0,0,4,4,4,0,-4,0,0,-4,4,0,0,0,4,0,0,4,0,0,0,0,4,4,0,0,-4,0,0,0,0,0,-4,0,0,-4,0,0,-4,0,-4,0,0,4,0,-4,-4,-4,0,0,0,0,0,0,-4,0,0,0,0,0,0,0,4,-4,0,0,0,-4,-4],[0,1,0,2],[80,80,32]],[3,[-6,-2,-2,2,2,2,2,2,-2,2,-2,2,2,-2,2,2,2,2,2,2,-2,2,-2,2,-2,-2,-2,-2,2,2,-2,-2,-2,-2,2,-2,-2,2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-6,4,0,0,0,0,0,0,0,4,0,-4,4,0,0,-4,4,0,0,0,4,0,0,-4,0],[2,1,3],[64,64,64]],[6,[2,-2,-2,2,2,2,2,2,2,2,2,2,-2,-2,-2,2,2,-6,2,2,2,-2,2,-2,2,-2,-2,2,-2,2,-2,-6,2,-2,-2,-2,2,-2,-2,-2,2,2,-2,2,-2,2,-2,2,-4,0,0,0,4,0,0,0,0,0,0,0,0,-4,-4,-4,0,0,-4,4,-4,0,0,0],[0,1,2,1],[48,80,64]]];
assert.deepEqual(literal.map(([a])=>a),[4,3,6]);
for(const [a,y,path,byp] of literal){
 integralPhase(y);
 const p=path.reduce((s,i)=>gen[i](s),triples[2]);assert.deepEqual(p,triples[a]);
 const Y=iso(y,path);
 const bp=[0,1,2].map(k=>dot(X.slice(k*24,k*24+24),Y.slice(k*24,k*24+24)));
 assert.deepEqual(bp,byp);assert.equal(dot(X,Y),192);assert.equal(dot(X.map((v,i)=>v-Y[i]),X.map((v,i)=>v-Y[i])),384);
}
console.log('PASS rank-two global palette cap 8; all 3 unordered F8 ratio types attain actual horizontal squared distance 6.');
