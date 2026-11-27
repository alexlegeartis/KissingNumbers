
const enc="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
function un64(s){
 const a=new Uint8Array(s.length*3/4|0), idx=new Int16Array(128).fill(-1);
 for(let j=0;j<64;j++)idx[enc.charCodeAt(j)]=j;
 let n=0,v=0,k=0;for(let j=0;j<s.length;j++){const c=s.charCodeAt(j),r=c<128?idx[c]:-1;if(r<0)continue;v=(v<<6)|r;k+=6;if(k>=8){k-=8;a[n++]=v>>>k;}}
 return a.subarray(0,n);
}
const ascii=(a,s,t)=>Array.from(a.subarray(s,t),v=>String.fromCharCode(v)).join("");
const u16=(a,p)=>a[p]|a[p+1]<<8;
const u32=(a,p)=>(a[p]|a[p+1]<<8|a[p+2]<<16|a[p+3]<<24)>>>0;
const u64=(a,p)=>u32(a,p)+2**32*u32(a,p+4);
function unzip(a){
 const files=new Map();
 let k=0;
 while(k+30<a.length&&u32(a,k)===0x04034b50){
  const type=u16(a,k+8),flags=u16(a,k+6),preLen=u32(a,k+18),outLen=u32(a,k+22),n=u16(a,k+26),m=u16(a,k+28);
  if(flags&8)throw Error("zip data descriptor not supported");
  const f=ascii(a,k+30,k+30+n), extra=k+30+n; let c=preLen,o=outLen;
  if(c===0xffffffff||o===0xffffffff){
   let t=extra;while(t<extra+m&&u16(a,t)!==1)t+=4+u16(a,t+2);
   if(t>=extra+m)throw Error("zip64 missing");
   t+=4;if(o===0xffffffff){o=u64(a,t);t+=8;}if(c===0xffffffff)c=u64(a,t);
  }
  const begin=extra+m;
  if(type===8)files.set(f,inflate(a.subarray(begin,begin+c),o));
  else if(type===0)files.set(f,a.slice(begin,begin+c));else throw Error("unknown method");
  if(files.get(f).length!==o)throw Error("output length");
  k=begin+c;
 }
 return files;
}
function inflate(a,N){
 let off=0,C=0,n=0,O=0;
 const out=new Uint8Array(N);
 function bits(b){while(n<b){if(off>=a.length)throw Error("input underrun");C|=a[off++]<<n;n+=8;}
    const r=C&((1<<b)-1);C>>>=b;n-=b;return r;}
 function huff(L){
  const M=Math.max(...L);
  const S=new Uint16Array(M+1), D=new Uint16Array(M+1),H=new Uint32Array(1<<M);
  for(const l of L)if(l)S[l]++;
  for(let i=1;i<=M;i++)D[i]=(D[i-1]+S[i-1])<<1;
  for(let s=0;s<L.length;s++){
   const t=L[s];if(!t)continue;const u=D[t]++;let v=0;for(let i=0;i<t;i++)v=(v<<1)|((u>>>i)&1);
   for(let z=v;z<H.length;z+=1<<t)H[z]=(t<<16)|s;
  }
  return [M,H];
 }
 function get(T){
  const [M,H]=T;while(n<M&&off<a.length){C|=a[off++]<<n;n+=8;}
  const z=H[C&((1<<M)-1)],k=z>>>16;if(!k||k>n)throw Error("bad Huffman");
  n-=k;C>>>=k;return z&65535;
 }
 const perm=[16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15];
 const Lbase=[3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258];
 const Lextra=[0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0];
 const Dbase=[1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577];
 const Dextra=[0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13];
 let Lf,Df;
 for(let more=1;more;){
   more=1-bits(1);const type=bits(2);
   if(!type){bits(n%8);const len=bits(16),v=bits(16);if((len^v)!==65535)throw Error("stored len");
    for(let j=0;j<len;j++)out[O++]=bits(8);continue;}
   let L,D;
   if(type===1){Lf||=huff(Array.from({length:288},(_,i)=>i<144?8:i<256?9:i<280?7:8));Df||=huff(new Array(32).fill(5));L=Lf;D=Df;}
   else if(type===2){const nl=bits(5)+257,nd=bits(5)+1,hn=bits(4)+4,cl=new Array(19).fill(0);
     for(let j=0;j<hn;j++)cl[perm[j]]=bits(3);
     const H=huff(cl),F=[];
     while(F.length<nl+nd){const s=get(H);if(s<16){F.push(s);continue;}
      const d=s===16?3+bits(2):s===17?3+bits(3):11+bits(7),v=s===16?F[F.length-1]:0;
      for(let j=0;j<d;j++)F.push(v);}
     if(F.length!==nl+nd)throw Error("code length");
     L=huff(F.slice(0,nl));D=huff(F.slice(nl));
   }else throw Error("reserved DEFLATE type");
   for(;;){const s=get(L);if(s<256){out[O++]=s;continue;}if(s===256)break;
    if(s>285)throw Error("bad len");
    const d=s-257,k=Lbase[d]+bits(Lextra[d]),t=get(D);
    if(t>=30)throw Error("bad distance code");
    const j=Dbase[t]+bits(Dextra[t]);if(j>O||k+O>N)throw Error("bounds");
    for(let u=0;u<k;u++,O++)out[O]=out[O-j];
   }
 }
 if(O!==N)throw Error("inflate length "+O+" not "+N);return out;
}
function npy(a){
 if(a[0]!==0x93||ascii(a,1,6)!=="NUMPY")throw Error("npy");
 const m=a[6]===1?10:12,len=m===10?u16(a,8):u32(a,8),h=ascii(a,m,m+len);
 const shape=(h.match(/'shape': *\((.*?)\)/)||[0,''])[1].split(',').filter(x=>x.trim()).map(Number);
 const type=(h.match(/'descr': *'(.*?)'/)||[0,''])[1];
 if(h.includes("'fortran_order': True")||type!=="<i2")throw Error(h);
 const z=new Int16Array(a.buffer,a.byteOffset+m+len,(a.length-m-len)/2);
 if(shape.reduce((r,t)=>r*t,1)!==z.length)throw Error(h);
 return {shape,z};
}

const dot=(x,y,m=8)=>{let v=0;for(let i=0;i<m;i++)v+=x[i]*y[i];return v};
const key=x=>x.join(",");
const Tnorm=t=>t.reduce((v,w)=>v+w*w,0);
function extWall(E){
 const X=[],B=new Map();for(const [x,m,n]of E){if(!n[2]){X.push(x.slice(0,8));continue;}
 if(!n[0])continue;
 if(n[2]!==32)throw Error("unknown tail type");
 const t=x.slice(8).map(v=>v/4);let f=t.find(v=>v);if(f<0)continue;
 const k=key(t);if(!B.has(k))B.set(k,[t,[]]);B.get(k)[1].push(x.slice(0,8));
 }
 return [X,[...B.values()]];
}
function* signedShape(shape,d){
 const c=new Int8Array(d),cnt={};for(const x of shape)cnt[x]=(cnt[x]||0)+1;
 const V=Object.keys(cnt).map(Number);if(V.includes(0))throw Error("shape");
 cnt[0]=d-shape.length;
 function* gen(j){
  if(j===d){yield [...c];return;}
  if(cnt[0]){cnt[0]--;c[j]=0;yield*gen(j+1);cnt[0]++;}
  for(const x of V)if(cnt[x]){cnt[x]--;c[j]=-x;yield*gen(j+1);c[j]=x;yield*gen(j+1);cnt[x]++;}
 }
 yield*gen(0);
}
function walls(a,P){
 const [X,B]=P;
 for(const x of X)if(dot(a,x)>64)return null;
 const h=[];
 for(const [t,V]of B){
  let m=32;for(const x of V){let u=dot(a,x);if(u>m){m=u;if(m>64)return null;}}
  h.push((64-m)/4);
 }
 return h;
}
function gate(E,C,shape,T,opt={}){
 const P=extWall(E),tails=T.map(t=>P[1].map(([p])=>Math.abs(dot(t,p,4))));
 let na=0,sp=0,nh=0,nt=0,good=[],hist={},found=[];
 for(const a of signedShape(shape,8)){
  ++na;
  const h=walls(a,P);if(!h)continue;++sp;
  const poss=tails.map((q,i)=>q.every((v,j)=>v<=h[j])?i:-1).filter(i=>i!==-1);
  if(!poss.length)continue;++nh;nt+=poss.length;
  if(opt.gateOnly)continue;
  for(const i of poss){let x=a.concat(T[i]);const bad=[];for(let j=0;j<C.length;j++)if(dot(x,C[j],12)>64)bad.push(j);
   hist[bad.length]=(hist[bad.length]||0)+1; if(bad.length<=(opt.maxBad??24))found.push([x,bad]);
 }
 }
 return {shape,na,sp,nh,nt,hist,found}
}

const det3=M=>M[0]*(M[4]*M[8]-M[5]*M[7])-M[1]*(M[3]*M[8]-M[5]*M[6])+M[2]*(M[3]*M[7]-M[4]*M[6]);
function geom(P){
 const A=P[1].map(q=>q[0]),mat=[];
 for(let a=0;a<9;a++)for(let b=a+1;b<10;b++)for(let c=b+1;c<11;c++)for(let d=c+1;d<12;d++){
 const I=[a,b,c,d],R=I.map(i=>A[i]), cof=[];
 for(let i=0;i<4;i++){cof[i]=[];for(let j=0;j<4;j++){
  const mm=[];for(let k=0;k<4;k++)if(k!==i)for(let l=0;l<4;l++)if(l!==j)mm.push(R[k][l]);
  cof[i][j]=((i+j)&1?-1:1)*det3(mm);
 }}
 let D=R[0].reduce((v,x,i)=>v+x*cof[0][i],0);
 if(!D)continue;
 if(D<0){D=-D;for(let i=0;i<4;i++)for(let j=0;j<4;j++)cof[i][j]*=-1;}
 mat.push({I,cof,D});
 }
 return {A,mat};
}
function exactMax({A,mat},b){
 let max=[0,1,[0,0,0,0]],n=0;
 for(const {I,cof,D} of mat)for(let bits=0;bits<8;bits++){
  const rhs=I.map((k,i)=>((bits>>i)&1)?-b[k]:b[k]),x=[0,0,0,0];
  for(let i=0;i<4;i++)for(let j=0;j<4;j++)x[i]+=cof[j][i]*rhs[j];
  if(A.some((u,i)=>Math.abs(u.reduce((a,v,j)=>a+v*x[j],0))>D*b[i]))continue;
  n++;const w=x.reduce((s,n)=>s+n*n,0);
  if(w*max[1]>max[0]*D*D)max=[w,D*D,x,D];
 }
 return {max,n};
}

const assert=(x,m)=>{if(!x)throw Error(JSON.stringify(m));};
const Owners=[202,202,203,203,202,202,203,203,202,202,203,203,44,44,182,45,183,45,79,269,44,238,238,45,239,239,44,44,182,45,183,45,79,269,44,238,238,45,239,239,44,210,182,45,183,211,79,269,44,238,238,45,239,239,202,202,203,203,202,202,203,203,202,202,203,203,44,44,182,45,183,45,79,269,44,238,238,45,239,239,44,100,182,45,183,211,79,269,44,238,238,45,239,239,44,210,182,45,183,211,79,269,44,238,238,45,239,239,108,108,109,109,108,100,109,109,108,108,109,109,110,110,111,111,110,110,111,111,110,110,111,111,250,88,250,267,267,89,250,250,250,267,267,101,250,250,250,267,267,89,250,250,250,267,267,263,250,250,88,267,267,263,250,250,250,267,267,101,202,202,203,203,202,202,203,203,202,202,203,203,44,44,182,45,183,45,79,269,44,238,238,45,239,239,44,182,182,45,183,183,79,269,44,238,238,45,239,239,44,210,182,45,183,211,79,269,44,238,238,45,239,239,202,202,203,203,202,202,203,203,202,202,203,203,44,182,182,45,183,183,79,269,44,238,238,45,239,239,44,100,182,45,183,211,218,269,44,238,238,45,239,239,210,210,182,211,183,211,79,269,210,238,238,211,239,239,250,250,250,267,267,263,250,250,250,267,267,101,250,250,250,267,267,89,250,250,250,267,267,263,250,250,250,267,267,263,250,250,250,267,267,101,238,238,88,88,89,89,239,239,194,194,270,195,271,195,88,238,88,88,89,89,239,239,194,194,270,195,271,195,238,238,88,88,89,89,239,239,194,270,270,195,271,271,88,88,100,88,101,101,89,101,194,194,270,195,271,195,238,238,100,88,89,101,239,239,194,270,270,195,271,271,88,88,88,88,89,89,89,89,194,270,270,195,271,271,238,238,100,100,101,101,239,239,194,270,270,195,271,271,238,238,88,88,89,89,239,239,194,270,270,195,271,271,238,238,94,94,95,95,239,239,194,270,270,195,271,271,88,88,262,88,89,287,89,89,194,270,270,195,271,271,100,100,100,100,101,101,101,101,270,270,270,271,271,271,262,262,262,266,387,263,263,263,270,270,270,271,271,271,34,34,34,93,93,93,266,267,25,269,35,238,35,253,253,239,34,34,34,93,93,93,266,267,25,269,35,238,35,253,253,239,34,34,34,251,91,93,197,267,25,269,35,238,35,253,253,239,34,240,34,93,93,241,197,267,25,269,35,35,35,253,253,253,34,240,34,251,93,241,197,267,25,269,35,238,35,253,253,239,108,108,109,109,108,104,109,109,224,96,132,97,108,104,109,105,108,108,109,109,108,104,109,105,34,240,34,93,93,241,197,267,25,269,35,35,35,253,253,253,34,34,34,251,91,91,197,267,25,269,35,35,35,253,253,239,34,240,34,251,93,241,197,267,25,269,35,238,35,253,253,239,34,34,34,251,91,251,197,267,25,269,35,238,35,253,253,239,110,110,111,111,110,110,111,111,110,106,111,107,222,110,133,111,110,110,111,111,110,110,111,111,34,240,34,93,93,241,197,267,25,269,35,35,35,253,253,253,34,240,34,251,93,241,197,267,25,269,35,35,35,253,253,253,34,240,34,251,93,241,197,267,25,269,35,35,35,253,253,253,252,90,252,269,269,249,252,252,252,269,269,103,252,252,252,269,269,249,252,252,252,269,269,103,252,252,90,269,269,249,252,252,252,269,269,103,108,108,109,109,108,104,109,105,108,108,109,109,110,110,111,111,110,110,111,111,110,110,111,111,42,240,240,43,241,241,197,267,42,42,102,43,185,43,42,240,240,43,241,241,77,267,42,42,184,43,185,43,42,240,240,43,241,241,197,267,42,102,102,43,185,213,204,204,205,205,204,204,205,205,204,204,205,205,42,240,240,43,241,241,197,267,42,42,184,43,185,43,42,240,240,43,241,241,197,267,42,102,102,43,185,213,42,240,240,43,241,241,197,267,42,212,184,43,185,213,204,204,205,205,204,204,205,205,204,204,205,205,252,252,252,269,269,265,252,252,252,269,269,103,252,252,252,269,269,249,252,252,252,269,269,103,252,252,252,269,269,249,252,252,252,269,269,103,42,240,240,43,241,241,197,267,42,42,102,43,185,43,42,240,240,43,241,241,197,267,42,184,184,43,185,185,42,240,240,43,241,241,77,267,42,102,102,43,185,213,204,204,205,205,204,204,205,205,204,204,205,205,42,240,240,43,241,241,197,267,42,184,184,43,185,185,42,240,240,43,241,241,197,267,42,102,102,43,185,213,212,240,240,213,241,241,197,267,212,212,184,213,185,213,204,204,205,205,204,204,205,205,204,204,205,205,272,272,272,273,273,273,264,264,264,268,389,265,265,265,272,272,272,273,273,273,102,102,102,102,103,103,103,103,272,272,272,273,273,273,90,90,102,90,91,91,91,91,272,272,272,273,273,273,240,240,92,92,93,265,241,241,272,272,272,273,273,273,240,240,90,90,91,91,241,241,272,272,272,273,273,273,240,240,102,102,93,103,241,241,272,272,272,273,273,273,90,90,90,90,91,91,91,91,272,272,272,273,273,273,240,240,102,90,91,103,241,241,272,272,272,273,273,273,90,90,102,90,91,91,91,91,272,272,272,273,273,273,240,240,90,90,91,265,241,241,272,272,272,273,273,273,90,240,90,90,91,91,241,241,272,272,272,273,273,273,240,240,90,90,91,91,241,241];
function run(buffer) {
 const Z=npy(unzip(buffer).get("points.npy"));
 assert(Z.shape[0]===223528&&Z.shape[1]===28,"literal live 28D witness");
 const I=[0,2,4,5,6,10,11,23], J=Array.from({length:24},(_,i)=>i).filter(i=>!I.includes(i)),T=[24,25,26,27];
 const E=new Map(), C=[], P=[];
 for(let n=0;n<Z.shape[0];n++){
  const k=28*n,a=I.map(i=>Z.z[k+i]),t=T.map(i=>Z.z[k+i]), b=J.map(i=>Z.z[k+i]), B=[Tnorm(a),Tnorm(b),Tnorm(t)];
  assert(B[0]+B[1]+B[2]===128&&a.concat(b,t).every(x=>x%2===0),"exact source norm and even alphabet");
  const x=a.concat(t),k1=key(x);
  if(!B[1]){C.push(x);if(!B[0])P.push(t);continue;}
  if(E.has(k1)){E.get(k1)[1]++;continue;}
  E.set(k1,[x,1,B]);
 }
 assert(C.length===712&&E.size===3833&&new Set(C.map(key)).size===712,"old core + projected exterior");
 const poles=[...signedShape([8,8],4)], two=[...signedShape([4,4],4)],tk=new Set(two.map(key)); 
 assert(P.length===24&&new Set(P.map(key)).size===24&&poles.every(t=>P.some(x=>key(x)===key(t))),"all 24 original tail poles kept");
 for(const [x,_,r]of E.values()){
  const neg=x.slice(0,8).concat(x.slice(8).map(v=>-v));
  assert(E.has(key(neg))&&E.get(key(neg))[1]===_,"antipodal tail per projected body");
  assert(!r[2]||(r[0]&&r[2]===32&&tk.has(key(x.slice(8))))||(!r[0]&&r[2]===64),"only classified exterior tail radii");
 }
 const tails=[...signedShape([4,4,4],4),...signedShape([6,2,2,2],4)];
 const shapes=[];
 function add(bound,rem,q){if(!rem){shapes.push(q);return;}if(q.length===8)return;
  for(let i=Math.min(bound,Math.floor(Math.sqrt(rem)));i>0;i--)add(i,rem-i*i,[...q,2*i]);}
 add(4,20,[]);
 assert(shapes.length===6&&tails.length===96&&tails.every(t=>Tnorm(t)===48),"exhaustive even profiles");
 const out=shapes.map(s=>gate([...E.values()],C,s,tails,{maxBad:712}));

 const all=out.flatMap(o=>o.found).map(([z,old])=>[key(z),z,old]).sort((a,b)=>a[0]<b[0]?-1:a[0]>b[0]?1:0);
 assert(Owners.length===1200&&all.length===1200&&new Set(all.map(x=>x[0])).size===1200&&out.reduce((s,x)=>s+x.nt,0)===1200,"exhaustive donor total");
 const taken=new Map();
 for(const [i,[_,z,old]]of all.entries()){
  const c=Owners[i];assert(Number.isInteger(c)&&old.includes(c),"owner is an actual removed blocker");
  if(!taken.has(c))taken.set(c,[]);const B=taken.get(c);
  assert(B.every(q=>dot(z,q,12)>64),"same owner means conflicting new pair");
  B.push(z);
 }
 assert(taken.size===75,"owner certificate covers all available donors");
 const q=extWall([...E.values()]),A=geom(q),bmap=new Map(),caps={};
 for(const a of signedShape([8,4],8)){const h=walls(a,q);assert(h&&h.every(Number.isInteger),"integral wall normal equation");const k=key(h);if(!bmap.has(k))bmap.set(k,exactMax(A,h));const m=bmap.get(k).max;assert(m[1]>0&&m[0]<48*m[1],"arbitrary real tails too short");const n=m[0]/m[1];assert(Number.isInteger(n),"valid exact squared cap");caps[n]=(caps[n]||0)+1;}
 assert(caps[20]===148&&caps[24]===28&&caps[32]===48,"exact continuous radial obstruction");

 const B64=new Map(),E64=new Map();
 for(const x of C){const a=x.slice(0,8),t=x.slice(8);if(Tnorm(a)!==64||Tnorm(t)!==64)continue;
  const k=key(t);if(!B64.has(k))B64.set(k,[t,[]]);B64.get(k)[1].push(a);}
 for(const [x,m,r]of E.values())if(!r[0]&&r[2]===64)E64.set(key(x.slice(8)),m);
 const short=[...signedShape([8],4),...signedShape([4,4,4,4],4)].map(key);
 assert(B64.size===24&&E64.size===24&&short.every(k=>B64.has(k)&&E64.get(k)===32),"actual 384 orthoplex blocks and 768 exterior short tails");
 const F=[...B64.values()].map(([t,V])=>{
  assert(V.length===16&&new Set(V.map(key)).size===16&&V.every(a=>V.some(b=>key(b)===key(a.map(s=>-s)))),"full antipodal frame");
  const B=V.filter(a=>a.find(x=>x)!==undefined&&a.find(x=>x)>0);
  assert(B.length===8&&B.every((a,i)=>B.every((b,j)=>dot(a,b)===(i===j?64:0))),"orthogonal basis");
  return [t,B];
 });
 for(const t of signedShape([24,24],4)){
  const C6=F.filter(([y])=>dot(t,y,4)===192);assert(C6.length===6,"two axis and four sign head maximizers");
  let best=Infinity;
  for(let s=0;s<256;s++){
   const A=Array.from({length:8},(_,j)=>2*C6[0][1].reduce((u,v,i)=>u+((s>>i&1)?-1:1)*v[j],0));
   assert(Tnorm(A)===2048,"scaled corner has body radius 2048/25");
   const m=Math.max(...C6.flatMap(([_,B])=>B.map(a=>Math.abs(dot(A,a)))));
   best=Math.min(best,m);
  }
  assert(best===256,"six-frame corner enumeration");
 }
 // Rational endpoints in the real-radius margin.
 const sq=x=>x*x;
 assert(sq(479n)<23n*10000n&&4n*100n<479n,"S >= 479/100, 2*Delta/S < 1/25");
 assert(2n*sq(24n)-46n*25n===2n,"Delta <= 2/25");
 assert(3n*(479n-4n)>2n*479n&&sq(24n)+3n<46n*25n,"exactly two high head slots");
 assert(2n*(6400n-8n*479n+64n)<53n*100n,"each body-basis product < 53/2");
 assert(4n*64n*2048n-25n*7n*sq(53n)>25n*4n*sq(18n),"all first-frame absolute products >18");
 assert(8n*sq(38n)<25n*sq(22n)&&2n*(256n-22n*5n)>53n*5n,"separating six-frame corner has ample margin");
 console.log(JSON.stringify({status:"PASS",strictRealBodyFloor:82,full:Z.shape[0],exterior:E.size,oldCore:C.length,donors:1200,owners:taken.size,cap84:caps,byShape:out.map(o=>[o.shape,o.na,o.nh,o.nt])}));
}

// Incremental exact replay against the pinned complete coordinate archive.
import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
const filename=process.argv[2];
if(!filename)throw Error('Usage: node verify.mjs /path/to/coordinates28_223528.npz');
const bytes=readFileSync(filename);
if(createHash('sha256').update(bytes).digest('hex')!=='3556adb1c32b26116bd83bf0793e4eca0b007ca1c751da47c50d3e37cb66d5d2')throw Error('Incorrect pinned archive');
run(bytes);
