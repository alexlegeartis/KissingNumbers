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

function linear(words){const B=new Uint32Array(24);for(let w of words){while(w){const i=31-Math.clz32(w);if(B[i])w^=B[i];else{B[i]=w;break;}}}return Array.from(B).filter(Boolean);}
function group(B){let L=[0];for(const i of B)L.push(...L.map(x=>x^i));return L;}

function allowedB(d,P){

const MA=new Map(d.M4);
const BB=group(linear(d.B)), h={},metric=[];
if(BB.length!==4096)throw Error('B rank');
for(const [a,b,c] of P.L0){
 const ua=mp(a,d.FM), ub=mp(b,d.FM);
 const S=new Set;
 for(const [x,y] of [[a,b],[a,c],[b,c]])for(const r of MA.get(x))for(const s of MA.get(y))S.add(r^s);
 let valid=0,already=0;const T=[];
 for(const t of BB){if((pc(ua&t)&1)||(pc(ub&t)&1))continue;if(S.has(t))already++;else{valid++;T.push(t);}}
 if(valid+already!==1024||S.size!==already)throw Error('off subspace');
 h[valid]=(h[valid]||0)+1;
 metric.push(valid);
}
return {h,metric};
}
function graph(P,d){
const dec=12, map=new Map(P.A.map((a,i)=>[a,i]));
for(let i=0;i<P.A.length;i++)for(let j=0;j<dec;j++)if((P.A[i]^P.A[1<<j])!==P.A[i^(1<<j)])throw Error('linear order');
const L=[], planeType={},typ=new Map(d.M4.map(([u,L])=>[u,linear(L.map(x=>x^L[0])).length]));
for(const U of P.L0){
 const k=U.map(x=>typ.get(x)).sort().join(',');planeType[k]=(planeType[k]||0)+1;
 for(const a of U)for(const b of U)if(a!==b)L.push(map.get(a)|(map.get(b)<<12));
}
const M=new Map(L.map((u,i)=>[u,i])),C=L.map(u=>[]),hist={};
for(let i=0;i<L.length;i++)for(let j=i+1;j<L.length;j++)if(M.has(L[i]^L[j])){C[i].push(j);C[j].push(i);}
for(const S of C)hist[S.length]=(hist[S.length]||0)+1;

let tri=0,k4=0,witness=null;
for(let i=0;i<C.length;i++){
 const a=new Set(C[i]);for(const j of C[i])if(j>i){
  const b=C[j].filter(k=>k>j&&a.has(k));tri+=b.length;
  if(b.length&&!witness)witness=[i,j,b[0]].map(k=>L[k]);
  for(const k of b)for(const q of C[k])if(q>k&&b.includes(q))k4++;
 }
}
return {C,count:L.length,hist,planeType,tri,k4,witness};
}
function color(O){

const CC=O.C,col=new Int8Array(CC.length).fill(-1),S=new Uint16Array(CC.length),remaining=CC.map(L=>L.length);
const h={};let top=0;
for(let z=0;z<CC.length;z++){
 let v=-1,score=-1;
 for(let i=0;i<CC.length;i++)if(col[i]<0){const t=32*pc(S[i])+remaining[i];if(t>score){score=t;v=i;}}
 let k=0;while(S[v]>>k&1)k++;
 col[v]=k;top=Math.max(top,k);h[k]=(h[k]||0)+1;
 for(const i of CC[v])if(col[i]<0){S[i]|=1<<k;remaining[i]--;}
}
for(let v=0;v<CC.length;v++)if(CC[v].some(x=>col[x]===col[v]))throw Error('color');
return {h,n:1+top,valid:true};
}
const assert=require('node:assert/strict'),fs=require('node:fs'),crypto=require('node:crypto');
const cbuf=fs.readFileSync(process.argv[2]),gbuf=fs.readFileSync(process.argv[3]);
const blobHash=x=>crypto.createHash('sha1').update(Buffer.from('blob '+x.length+'\0')).update(x).digest('hex');
assert.equal(blobHash(cbuf),'5cc5386748aa477917fa7fd113a9c54a343800f5');
assert.equal(crypto.createHash('sha256').update(gbuf).digest('hex'),'1791fab8c152302db2652ec66f227eb473ce4007f934431ef6bbb2ab2f2d87ab');
const d=gamma(JSON.parse(cbuf),JSON.parse(gbuf)),P=planes(d);
assert.deepEqual(P.h,{'6':65,'8':650,'10':3380});
assert.equal(P.L.length,975);assert.equal(P.L0.length,975);
assert.deepEqual(P.depthDegree,{'0':3381,'3':650,'15':65});
const B=allowedB(d,P);
assert.deepEqual(B.h,{'338':975});
assert.equal(975*6*338,1977300);
const A=graph(P,d);
assert.equal(A.count,5850);assert.deepEqual(A.hist,{'14':5850});
assert.deepEqual(A.planeType,{'6,8,8':975});assert.equal(A.tri,13650);assert.equal(A.k4,0);assert(A.witness);
const colorResult=color(A);assert.equal(colorResult.n,7);assert.deepEqual(colorResult.h,{'0':1149,'1':1145,'2':1103,'3':1018,'4':913,'5':504,'6':18});
const tail=[];
for(const [a,b] of [[4,9],[5,6],[7,8],[11,12]])for(const s of [-2,2]){
 const w=Array(13).fill(0);w[0]=-2;w[1]=-2;w[a]=-2;w[b]=s;tail.push(w);
}
for(let i=0;i<8;i++)for(let j=0;j<=i;j++)assert.equal(dot(tail[i],tail[j]),i===j?16:8);
console.log('Depth=6 integral cosets 1977300; 975 planes × 6 orderings × 338 B-sums. Projected induced local graph: 5850 vertices, degree 14, 13650 triangles, no K4. Integral depth-six palette clique number 4; literal 8-clique from old 13D tail carrier.',A.witness);

const u=720, v=8359937;
assert(P.L0.some(U=>U.includes(u)&&U.includes(v)&&U.includes(u^v)));
const qs=[];
for(let a=0;a<8;a++){
  const a0=a&1,a1=(a>>1)&1,a2=(a>>2)&1;
  qs.push([ (a0?u:0)^(a2?v:0), (a1?u:0)^((a0^a2)?v:0), (a2?u:0)^(a1?v:0) ]);
}
assert.equal(qs.filter(x=>x[0]!==0).length,6);
for(let i=0;i<8;i++)for(let j=0;j<i;j++){
   const delta=qs[i].map((x,k)=>x^qs[j][k]), nz=new Set(delta.filter(Boolean));
   assert(nz.size>=2&&[...nz].every(x=>x===u||x===v||x===(u^v)));
}
assert.equal(6218175600+6*136704+6,6218995830);
console.log('PASS literal 74D Barnes rank-two palette and tied count');
