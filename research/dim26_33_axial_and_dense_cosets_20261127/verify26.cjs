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

function snap26(doc){
 const d=doc.added_integer.data,m=doc.midpoint_rows.data,p=doc.fraction_numerators.data,fden=doc.fraction_denominators.data;
 const X=[],M=[],R=[],F=[];for(let i=0;i<d.length;i+=26)X.push(d.slice(i,i+26));
 const mids=new Set(m),pairs=new Map();
 for(let i=0;i<X.length;i++){
  const v=X[i];if(!mids.has(i))R.push(v);else{
   const z=v.slice(0,24),k=z.join(',');
   if(!pairs.has(k))pairs.set(k,{b:z,t:[]});pairs.get(k).t.push(v.slice(24));
  }
 }
 for(let i=0;i<p.length;i+=26)F.push(p.slice(i,i+26));
 const U=[...pairs.values()].map(({b,t})=>{if(t.length!==2||t[0].some((x,j)=>x!==-t[1][j]))throw Error('pair');return {b,col:t[0][0]===t[0][1]?0:1,t};});
 if(U.length!==4416||R.length!==196||F.length!==104)throw Error('census');
 return{M:U,R,F,fden};
}
function leech(){
 const H=new Set(),gp=(1<<11)|(1<<10)|(1<<6)|(1<<5)|(1<<4)|(1<<2)|1;const gl=[0];
 for(let i=0;i<12;i++){let g=(gp<<i)|(1<<23);gl.push(...gl.slice().map(h=>h^g));}
 const out=[];const z=()=>new Int8Array(24);
 for(let i=0;i<24;i++)for(let j=i+1;j<24;j++)for(const a of [-8,8])for(const b of [-8,8]){const x=z();x[i]=a;x[j]=b;out.push(x);}
 for(const w of gl)if(pc(w)===8){
  const c=[];for(let i=0;i<24;i++)if(w>>i&1)c.push(i);
  for(let a=0;a<256;a++)if(!(pc(a)&1)){const x=z();for(let j=0;j<8;j++)x[c[j]]=a>>j&1?-4:4;out.push(x);}
 }
 for(const w of gl){
  const r=Int8Array.from({length:24},(_,i)=>w>>i&1?2:-2);
  for(let j=0;j<24;j++){const x=r.slice();x[j]*=-3;out.push(x);}
 }
 out.sort((a,b)=>{for(let i=0;i<24;i++)if(a[i]!==b[i])return a[i]-b[i];return 0;});
 if(out.length!==196560)throw Error('shell count');return out;
}

function vacant26(doc,LL,K){
 const cp=doc.midpoint_endpoints.data,x=doc.added_integer.data,mi=doc.midpoint_rows.data;
 for(let j=0;j<mi.length;j++){
  const u=LL[cp[2*j]],v=LL[cp[2*j+1]],k=mi[j]*26;
  if(!u||!v||!u.every((w,l)=>w+v[l]===2*x[k+l]))throw Error('coordinate convention '+j);
 }
 const present=new Set(K.M.map(x=>x.b.map(v=>2*v).join(',')));
 const D=doc.deleted_leech_indices.data, Z=D.map(j=>LL[j]);
 const seen=new Map();for(let i=0;i<Z.length;i++){const u=Z[i];
  for(let j=i+1;j<Z.length;j++){const v=Z[j];let d=0;for(let k=0;k<24;k++)d+=u[k]*v[k];if(d!==64)continue;
   const b=Array.from(u,(x,k)=>x+v[k]),s=b.join(',');if(!present.has(s)&&!seen.has(s))seen.set(s,[b,D[i],D[j]]);
  }
 }
 return [...seen.values()];
}

function coreGraph(K){
 const P=K.M.length,p=Int32Array.from({length:P},(_,i)=>i),s=new Uint16Array(P).fill(1);
 const find=(i)=>{let j=i;while(j!==p[j])j=p[j];while(i!==j){const t=p[i];p[i]=j;i=t;}return j;};
 const join=(i,j)=>{let a=find(i),b=find(j);if(a===b)return;if(s[a]<s[b])[a,b]=[b,a];p[b]=a;s[a]+=s[b];};
 function dot(a,b){let z=0;for(let i=0;i<24;i++)z+=a[i]*b[i];return z;}
 let E=0;for(let i=0;i<P;i++){const a=K.M[i];for(let j=i+1;j<P;j++){const b=K.M[j];if(a.col!==b.col&&dot(a.b,b.b)>32){join(i,j);E++;}}}
 const C=[],ids=new Map(),pin=[],Q=[[4,4],[4,-4]];
 for(let i=0;i<P;i++){
  const c=find(i);if(!ids.has(c)){ids.set(c,C.length);C.push([]);pin.push(false);}
  const g=ids.get(c),old=K.M[i];C[g].push(old);
  const t=Q[1-old.col];
  if(!pin[g])for(let j=0;j<K.F.length;j++){const v=K.F[j],h=dot(old.b,v);if(h+Math.abs(t[0]*v[24]+t[1]*v[25])>64*K.fden[j]){pin[g]=true;break;}}
 }
 const giant=C.map((x,i)=>[x.length,i]).sort((a,b)=>b[0]-a[0])[0][1];
 return{C,pin,E,giant,Q};
}
function near26(K,G,D){
 const old=G.C[G.giant].slice(),all=[];let w=5903;
 for(let i=old.length-1;i>0;i--){w=(Math.imul(w,1664525)+1013904223)>>>0;const j=w%(i+1);[old[i],old[j]]=[old[j],old[i]];}
 for(const [b,a,c] of D){
  let M=64,N=64,ok=true;for(const {b:t,col} of old){let d=0;for(let k=0;k<24;k++)d+=b[k]*t[k];
   if(d>128){ok=false;break;}
   if(col===0){const h=128-d;if(h<M){M=h;if(M*M+N*N<4096){ok=false;break;}}}
   else{const h=128-d;if(h<N){N=h;if(M*M+N*N<4096){ok=false;break;}}}
  }
  if(ok)all.push([b,a,c]);
 }
 return all;
}
function signLe(a,b,r){
 if(b===0n||r===0n)return a<=0n;
 if(a<=0n&&b<0n)return true;
 if(a>=0n&&b>0n)return false;
 if(b>0n)return b*b*r<=a*a;
 return a*a<=b*b*r;
}
function cuts26(b,K,G){
 const same=[];let near=true,all=[],de=[];
 for(const chunk of G.C){
  let V=[100000,100000];
  for(const {b:v,col} of chunk){let d=0;for(let j=0;j<24;j++)d+=v[j]*b[j];V[col]=Math.min(V[col],128-d);}
  same.push(V);if(Math.min(...V)<0)near=false;
 }
 if(!near)return null;
 const F=K.F.map((x,i)=>{let h=128*K.fden[i];for(let j=0;j<24;j++)h-=x[j]*b[j];return[x[24],x[25],h];});
 const ends=[[4,4,0,0,0,1]];
 function cut(x,y,h){
  const L=x*x+y*y,l=128*L-h*h;if(!L||l<0)return;
  ends.push([h*x,h*y,-y,x,l,2*L],[h*x,h*y,y,-x,l,2*L]);
 }
 for(const V of same)for(const h of V)if(h<=64)for(const [x,y] of G.Q){cut(x,y,h);cut(-x,-y,h);}
 for(const f of F)cut(...f);
 function test(ep){
  const [x0,y0,x1,y1,r,s]=ep.map(BigInt);
  const le=(x,y,t)=>{x=BigInt(x);y=BigInt(y);t=BigInt(t);
   return signLe(2n*(x*x0+y*y0)-t*s,2n*(x*x1+y*y1),r);};
  const ab=(x,y,h)=>le(x,y,h)&&le(-x,-y,h);
  for(const [x,y,h] of F)if(!le(x,y,h))return null;
  const flip=[];
  for(let i=0;i<same.length;i++){const [u,v]=same[i];
   if(ab(4,4,u)&&ab(4,-4,v))continue;
   if(!G.pin[i]&&ab(4,-4,u)&&ab(4,4,v)){flip.push(i);continue;}
   return null;
  }
  return flip;
 }
 const sols=[],seen=new Set();
 for(const ep of ends){const code=ep.join(',');if(seen.has(code))continue;seen.add(code);
  const f=test(ep);if(f!==null)sols.push([ep,f]);}
 return {checks:seen.size,sols};
}

function anglePrepare(doc){
const K=snap26(doc),A=K.R.filter(r=>r[24]===8&&r[25]===0),B=K.R.filter(r=>r[25]===8&&r[24]===0),J=[1,3,7,8,9,12,13,14,15,16,17,18,19,20,21,22];const C=[];
function dot24(x,y){return x.slice(0,24).reduce((s,x,i)=>s+x*y[i],0);}
for(const [j,F] of [[24,A],[25,B]]){
 assert.equal(F.length,48);
 const H=K.R.filter(r=>r[j]===-8&&r[49-j]===0);
 assert.deepEqual(H.map(v=>v.slice(0,24).join()).sort(),F.map(v=>v.slice(0,24).join()).sort());
 const S=new Set(F.map(v=>v.slice(0,24).join()));
 const chosen=F.filter(v=>v.slice(0,24).find(Boolean)>0);
 assert.equal(S.size,48);assert.equal(chosen.length,24);
 for(let a=0;a<chosen.length;a++){
  assert(S.has(chosen[a].slice(0,24).map(y=>-y).join()));
  for(let b=0;b<=a;b++)assert.equal(dot24(chosen[a],chosen[b]),a===b?64:0);
 }
}
let count4=0,count2=0;
function maxdot(a,rows){let m=0;for(const z of rows){let d=0;for(const [j,x] of a)d+=z[j]*x;m=Math.max(m,Math.abs(d));}return m;}
function gate(a,type){
 if(type===4)count4++;else count2++;
 const x=maxdot(a,A),y=maxdot(a,B);
 if((64-x)**2+(64-y)**2>=4096)C.push([a,x,y,type]);
}
for(let a=0;a<24;a++)for(let b=a+1;b<24;b++)for(let c=b+1;c<24;c++)for(let d=c+1;d<24;d++)for(let u=0;u<8;u++)gate([[a,4],[b,u&1?-4:4],[c,u&2?-4:4],[d,u&4?-4:4]],4);
for(let u=0;u<32768;u++)if(!(pc(u)&1))gate(J.map((j,i)=>[j,u>>i&1?-2:2]),2);
assert.equal(count4,85008);assert.equal(count2,16384);
return {C};
}
function angleTest(doc,prep){
const K=snap26(doc),G=coreGraph(K);if(!G.pin[G.giant])throw Error('giant unpinned');
const C=G.C[G.giant].slice();let w=77377;
for(let i=C.length-1;i>0;i--){w=(Math.imul(w,1664525)+1013904223)>>>0;const j=w%(i+1);[C[i],C[j]]=[C[j],C[i]];}
const H={},pass=[];
function get(a,b){let d=0;for(const [i,x] of a)d+=x*b[i];return d;}
for(const [a,x,y,type] of prep.C){
 if(x!==16||y!==16)throw Error('wrong angular gate');
 for(const s of [-1,1])for(const c of [0,1]){
  if(C.some(o=>o.col===c&&s*get(a,o.b)>=19))continue;
  const name=[type,s,c];H[name]=(H[name]||0)+1;pass.push([a,s,c,type]);
 }
}
return {giant:C.length,H,pass};
}
const assert=require('node:assert/strict');
const fs=require('node:fs'),crypto=require('node:crypto');
const b=fs.readFileSync(process.argv[2]);
assert.equal(crypto.createHash('sha256').update(b).digest('hex'),'35fd8143cf0c7e6b0146719485e00c86952cbd5e404c7881c5f6b92fa25b21ec');
const doc=decode(new Uint8Array(b)),K=snap26(doc),G=coreGraph(K);
assert.equal(G.C.length,18);assert.equal(G.E,86976);assert(G.pin[G.giant]);assert.equal(G.C[G.giant].length,4096);
// Both giant colors are centrally symmetric radial-96 exact tight frames.
for(const c of [0,1]){
 const L=G.C[G.giant].filter(x=>x.col===c).map(x=>x.b);
 assert.equal(L.length,2048);
 const S=new Set(L.map(x=>x.join()));
 assert.equal(S.size,2048);
 for(const x of L){
  assert.equal(x.reduce((a,b)=>a+b*b,0),96);
  assert(S.has(x.map(z=>-z).join()));
 }
 for(let i=0;i<24;i++)for(let j=i;j<24;j++)
  assert.equal(L.reduce((a,x)=>a+x[i]*x[j],0),i===j?8192:0);
}
const D=vacant26(doc,leech(),K),N=near26(K,G,D);let valid=0,cuts=0;
assert.equal(D.length,367270);assert.equal(N.length,2432);
for(const [body] of N){const r=cuts26(body,K,G);if(r){valid++;cuts+=r.checks;assert.equal(r.sols.length,0);}}
assert.equal(valid,256);assert.equal(cuts,6848);
const P=anglePrepare(doc);
assert.equal(P.C.filter(z=>z[3]===4).length,12544);assert.equal(P.C.filter(z=>z[3]===2).length,256);
const A=angleTest(doc,P);assert.equal(A.pass.length,0);
console.log('26D: 367270 vacant bodies; 2432 pass giant cut; 2176 rejected by orientation-independent body cuts; 256 proceed; 6848 exact arc endpoints, none feasible. Independent angular branch: 12800 fixed-axis body lines, no candidate survives the pinned giant.');